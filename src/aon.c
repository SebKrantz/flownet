#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>

/*
 * All-or-nothing assignment: one Dijkstra search per origin node, parallel over
 * origins with OpenMP.
 *
 * The search reproduces igraph's igraph_get_shortest_paths_dijkstra() (C core
 * 0.10.17, as bundled with R igraph 2.x) step for step, so it returns the same
 * path as igraph::shortest_paths() even where several paths tie on cost:
 *   - the same two-way indexed max-heap of negated distances (igraph_2wheap:
 *     shift-up continues on ties, sink prefers the left child on ties and swaps
 *     only on a strict improvement, modify sinks then shifts up at the old position);
 *   - the same relaxation: first reach when the stored distance is negative,
 *     improvement only when strictly shorter, parent edge recorded on both;
 *   - the same order of incident edges: by the other endpoint's id, and parallel
 *     edges in descending edge id (igraph's radix-sorted edge index);
 *   - the same early stop once every target of the origin has been settled.
 * Flows, path costs and edge counts are accumulated after the searches, in the
 * order of the former R loop (origins in order of first appearance, pairs in
 * input order, edges from origin to destination), so results do not depend on
 * the number of threads.
 */

#define PARENT(x)     (((x)+1)/2-1)
#define LEFTCHILD(x)  (((x)+1)*2-1)
#define RIGHTCHILD(x) (((x)+1)*2)

typedef struct {
  double *data;   /* heap values (negated distances) */
  int *index;     /* node at heap position */
  int *index2;    /* heap position + 2 of a node, 0 = not in the heap */
  int size;
} WHeap;

static inline void wh_switch(WHeap *h, int e1, int e2) {
  if(e1 != e2) {
    double tmp3 = h->data[e1];
    h->data[e1] = h->data[e2];
    h->data[e2] = tmp3;
    int tmp1 = h->index[e1], tmp2 = h->index[e2];
    h->index2[tmp1] = e2 + 2;
    h->index2[tmp2] = e1 + 2;
    h->index[e1] = tmp2;
    h->index[e2] = tmp1;
  }
}

static inline void wh_shift_up(WHeap *h, int elem) {
  while(!(elem == 0 || h->data[elem] < h->data[PARENT(elem)])) {
    int p = PARENT(elem);
    wh_switch(h, elem, p);
    elem = p;
  }
}

static inline void wh_sink(WHeap *h, int head) {
  for(;;) {
    int l = LEFTCHILD(head), r = RIGHTCHILD(head);
    if(l >= h->size) return;
    if(r == h->size || h->data[l] >= h->data[r]) {
      if(h->data[head] < h->data[l]) { wh_switch(h, head, l); head = l; } else return;
    } else {
      if(h->data[head] < h->data[r]) { wh_switch(h, head, r); head = r; } else return;
    }
  }
}

static inline void wh_push(WHeap *h, int idx, double elem) {
  int size = h->size;
  h->data[size] = elem;
  h->index[size] = idx;
  h->index2[idx] = size + 2;
  h->size = size + 1;
  wh_shift_up(h, size);
}

static inline double wh_delete_max(WHeap *h) {
  double tmp = h->data[0];
  int tmpidx = h->index[0];
  wh_switch(h, 0, h->size - 1);
  h->size--;
  h->index2[tmpidx] = 0;
  wh_sink(h, 0);
  return tmp;
}

static inline void wh_modify(WHeap *h, int idx, double elem) {
  int pos = h->index2[idx] - 2;
  h->data[pos] = elem;
  wh_sink(h, pos);
  wh_shift_up(h, pos);
}

/* other endpoint of 0-based edge e seen from 0-based node v */
#define OTHER(e, v) ((f[e] - 1 == (v)) ? t[e] - 1 : f[e] - 1)

/**
 * All-or-nothing assignment of an OD matrix
 *
 * @param from,to   Integer vectors (1-based node indices, length E)
 * @param cost      Numeric vector of non-negative edge costs (length E)
 * @param n_nodes   Number of nodes (max node index)
 * @param directed  Logical: directed graph (from -> to only) or undirected
 * @param od_from,od_to Integer vectors (1-based node indices, length P)
 * @param flow      Numeric vector of flows (length P)
 * @param retflags  Logical vector of length 3: return paths, path costs, edge counts
 * @param nthreads  Number of OpenMP threads
 * @return list(final_flows, found, paths, costs, counts); found is FALSE for pairs with
 *         an empty path (origin = destination or unreachable); paths are integer vectors of
 *         1-based edge ids from origin to destination
 */
SEXP aon_assign(SEXP from, SEXP to, SEXP cost, SEXP n_nodes, SEXP directed,
                SEXP od_from, SEXP od_to, SEXP flow, SEXP retflags, SEXP nthreads) {

  const int E = length(from), nn = asInteger(n_nodes), dir = asLogical(directed);
  const int P = length(od_from), nt = asInteger(nthreads) > 1 ? asInteger(nthreads) : 1;
  if(length(to) != E || length(cost) != E) error("from, to and cost must have the same length");
  if(length(od_to) != P || length(flow) != P) error("od_from, od_to and flow must have the same length");
  if(length(retflags) != 3) error("retflags must be a logical vector of length 3");
  const int ret_paths = LOGICAL(retflags)[0], ret_costs = LOGICAL(retflags)[1], ret_counts = LOGICAL(retflags)[2];
  const int *f = INTEGER(from), *t = INTEGER(to), *of = INTEGER(od_from), *ot = INTEGER(od_to);
  const double *w = REAL(cost), *fl = REAL(flow);

  for(int e = 0; e < E; e++) {
    if(ISNAN(w[e])) error("Weights must not contain NaN values.");
    if(w[e] < 0) error("Weights must not be negative, got %g.", w[e]);
    if(f[e] < 1 || f[e] > nn || t[e] < 1 || t[e] > nn) error("node index out of range");
  }
  for(int p = 0; p < P; p++)
    if(of[p] < 1 || of[p] > nn || ot[p] < 1 || ot[p] > nn) error("OD node index out of range");

  /* Incidence lists: one record per (node, other endpoint, edge), stably sorted by the
     other endpoint, then by node. Records are generated in descending edge id, so parallel
     edges keep that order. */
  const size_t R = dir ? (size_t) E : 2 * (size_t) E;
  int *rn = (int *) R_alloc(R, sizeof(int)), *rb = (int *) R_alloc(R, sizeof(int)), *re = (int *) R_alloc(R, sizeof(int));
  size_t k = 0;
  for(int i = E - 1; i >= 0; i--) {
    rn[k] = f[i] - 1; rb[k] = t[i] - 1; re[k++] = i;
    if(!dir) { rn[k] = t[i] - 1; rb[k] = f[i] - 1; re[k++] = i; }
  }
  int *cnt = (int *) R_alloc((size_t) nn + 1, sizeof(int)), *ord = (int *) R_alloc(R, sizeof(int));
  memset(cnt, 0, ((size_t) nn + 1) * sizeof(int));
  for(size_t i = 0; i < R; i++) cnt[rb[i] + 1]++;
  for(int i = 0; i < nn; i++) cnt[i + 1] += cnt[i];
  for(size_t i = 0; i < R; i++) ord[cnt[rb[i]]++] = (int) i;
  int *astart = (int *) R_alloc((size_t) nn + 1, sizeof(int));
  int *aedge = (int *) R_alloc(R, sizeof(int)), *aoth = (int *) R_alloc(R, sizeof(int));
  memset(astart, 0, ((size_t) nn + 1) * sizeof(int));
  for(size_t i = 0; i < R; i++) astart[rn[i] + 1]++;
  for(int i = 0; i < nn; i++) astart[i + 1] += astart[i];
  memcpy(cnt, astart, ((size_t) nn + 1) * sizeof(int));
  for(size_t j = 0; j < R; j++) {
    int i = ord[j], pos = cnt[rn[i]]++;
    aedge[pos] = re[i];
    aoth[pos] = rb[i];
  }

  /* OD pairs grouped by origin, origins in order of first appearance */
  int *gof = (int *) R_alloc((size_t) nn, sizeof(int));
  for(int i = 0; i < nn; i++) gof[i] = -1;
  int *gid = (int *) R_alloc((size_t) P + 1, sizeof(int)), *gorig = (int *) R_alloc((size_t) P + 1, sizeof(int)), ng = 0;
  for(int p = 0; p < P; p++) {
    int o = of[p] - 1;
    if(gof[o] < 0) { gof[o] = ng; gorig[ng++] = o; }
    gid[p] = gof[o];
  }
  int *gstart = (int *) R_alloc((size_t) ng + 1, sizeof(int)), *gmem = (int *) R_alloc((size_t) P + 1, sizeof(int));
  memset(gstart, 0, ((size_t) ng + 1) * sizeof(int));
  for(int p = 0; p < P; p++) gstart[gid[p] + 1]++;
  for(int g = 0; g < ng; g++) gstart[g + 1] += gstart[g];
  int *gpos = (int *) R_alloc((size_t) ng + 1, sizeof(int));
  memcpy(gpos, gstart, ((size_t) ng + 1) * sizeof(int));
  for(int p = 0; p < P; p++) gmem[gpos[gid[p]]++] = p;

  /* per pair: path length and offset into its origin's path buffer */
  int *plen = (int *) R_alloc((size_t) P + 1, sizeof(int)), **gbuf = (int **) R_alloc((size_t) ng + 1, sizeof(int *));
  size_t *poff = (size_t *) R_alloc((size_t) P + 1, sizeof(size_t));
  for(int g = 0; g < ng; g++) gbuf[g] = NULL;
  int failed = 0;

  #pragma omp parallel num_threads(nt) if(nt > 1)
  {
    double *dist = (double *) malloc((size_t) nn * sizeof(double));
    int *par = (int *) malloc((size_t) nn * sizeof(int)), *touched = (int *) malloc((size_t) nn * sizeof(int));
    char *is_target = (char *) calloc((size_t) nn, 1);
    WHeap h;
    h.data = (double *) malloc((size_t) nn * sizeof(double));
    h.index = (int *) malloc((size_t) nn * sizeof(int));
    h.index2 = (int *) calloc((size_t) nn, sizeof(int));
    h.size = 0;
    int ok = dist && par && touched && is_target && h.data && h.index && h.index2;
    if(!ok) {
      #pragma omp atomic write
      failed = 1;
    } else for(int i = 0; i < nn; i++) { dist[i] = -1.0; par[i] = 0; }

    #pragma omp for schedule(dynamic, 1)
    for(int g = 0; g < ng; g++) {
      if(!ok) continue;
      int src = gorig[g], ntouched = 0, to_reach = 0;
      for(int j = gstart[g]; j < gstart[g + 1]; j++) {
        int d = ot[gmem[j]] - 1;
        if(!is_target[d]) { is_target[d] = 1; to_reach++; }
      }
      dist[src] = 0.0; par[src] = 0; touched[ntouched++] = src;
      wh_push(&h, src, 0);
      while(h.size > 0 && to_reach > 0) {
        int u = h.index[0];
        double mindist = -wh_delete_max(&h);
        if(is_target[u]) { is_target[u] = 0; to_reach--; }
        for(int j = astart[u], jend = astart[u + 1]; j < jend; j++) {
          int e = aedge[j], v = aoth[j];
          double altdist = mindist + w[e], curdist = dist[v];
          if(curdist < 0) {
            dist[v] = altdist; par[v] = e + 1; touched[ntouched++] = v;
            wh_push(&h, v, -altdist);
          } else if(altdist < curdist) {
            dist[v] = altdist; par[v] = e + 1;
            wh_modify(&h, v, -altdist);
          }
        }
      }
      /* path lengths, then one buffer for the origin's paths */
      size_t total = 0;
      for(int j = gstart[g]; j < gstart[g + 1]; j++) {
        int p = gmem[j], d = ot[p] - 1, len = 0, act = d;
        is_target[d] = 0;
        if(d != src && par[d] != 0) {
          while(par[act]) { int e = par[act] - 1; act = OTHER(e, act); len++; }
        }
        plen[p] = len; poff[p] = total; total += (size_t) len;
      }
      int *buf = total ? (int *) malloc(total * sizeof(int)) : NULL;
      if(total && !buf) {
        #pragma omp atomic write
        failed = 1;
      } else {
        for(int j = gstart[g]; j < gstart[g + 1]; j++) {
          int p = gmem[j], act = ot[p] - 1;
          int *pp = buf + poff[p];
          for(int i = plen[p] - 1; i >= 0; i--) { int e = par[act] - 1; pp[i] = e + 1; act = OTHER(e, act); }
        }
      }
      gbuf[g] = buf;
      /* reset the workspace */
      while(h.size > 0) { h.size--; h.index2[h.index[h.size]] = 0; }
      for(int i = 0; i < ntouched; i++) { dist[touched[i]] = -1.0; par[touched[i]] = 0; }
    }
    free(dist); free(par); free(touched); free(is_target);
    free(h.data); free(h.index); free(h.index2);
  }

  if(failed) {
    for(int g = 0; g < ng; g++) free(gbuf[g]);
    error("aon_assign: out of memory");
  }

  /* results, accumulated in the order of the former R loop */
  SEXP res = PROTECT(allocVector(VECSXP, 5));
  SEXP final_flows = PROTECT(allocVector(REALSXP, E));
  SEXP found = PROTECT(allocVector(LGLSXP, P));
  SEXP paths = PROTECT(ret_paths ? allocVector(VECSXP, P) : R_NilValue);
  SEXP costs = PROTECT(ret_costs ? allocVector(REALSXP, P) : R_NilValue);
  SEXP counts = PROTECT(ret_counts ? allocVector(INTSXP, E) : R_NilValue);
  double *ff = REAL(final_flows), *pc = ret_costs ? REAL(costs) : NULL;
  int *fd = LOGICAL(found), *ec = ret_counts ? INTEGER(counts) : NULL;
  memset(ff, 0, (size_t) E * sizeof(double));
  if(ret_counts) memset(ec, 0, (size_t) E * sizeof(int));
  for(int g = 0; g < ng; g++) {
    for(int j = gstart[g]; j < gstart[g + 1]; j++) {
      int p = gmem[j], len = plen[p];
      const int *pp = gbuf[g] ? gbuf[g] + poff[p] : NULL;
      fd[p] = len > 0;
      if(len > 0) {
        double fv = fl[p];
        for(int i = 0; i < len; i++) ff[pp[i] - 1] += fv;
      }
      if(ret_costs) {
        if(len == 0) pc[p] = NA_REAL;
        else { double s = 0.0; for(int i = 0; i < len; i++) s += w[pp[i] - 1]; pc[p] = s; }
      }
      if(ret_counts) for(int i = 0; i < len; i++) ec[pp[i] - 1]++;
      if(ret_paths) {
        SEXP v = allocVector(INTSXP, len);
        SET_VECTOR_ELT(paths, p, v);
        if(len) memcpy(INTEGER(v), pp, (size_t) len * sizeof(int));
      }
    }
    free(gbuf[g]); gbuf[g] = NULL;
  }
  SET_VECTOR_ELT(res, 0, final_flows);
  SET_VECTOR_ELT(res, 1, found);
  SET_VECTOR_ELT(res, 2, paths);
  SET_VECTOR_ELT(res, 3, costs);
  SET_VECTOR_ELT(res, 4, counts);
  UNPROTECT(6);
  return res;
}
