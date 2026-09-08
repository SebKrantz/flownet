#ifdef _OPENMP
#include <omp.h>
#endif

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

typedef struct {
  int to;
  int edge;
  double cost;
  int next;
} Arc;

typedef struct {
  int node;
  int first_edge;
  double dist;
} HeapItem;

typedef struct {
  HeapItem *data;
  int size;
  int capacity;
  int overflow;
} MinHeap;

static void heap_push(MinHeap *heap, int node, int first_edge, double dist) {
  if(heap->size >= heap->capacity) {
    heap->overflow = 1;
    return;
  }

  int i = heap->size++;
  while(i > 0) {
    int parent = (i - 1) / 2;
    if(heap->data[parent].dist <= dist) break;
    heap->data[i] = heap->data[parent];
    i = parent;
  }

  heap->data[i].node = node;
  heap->data[i].first_edge = first_edge;
  heap->data[i].dist = dist;
}

static HeapItem heap_pop(MinHeap *heap) {
  HeapItem min_item = heap->data[0];
  HeapItem last = heap->data[--heap->size];

  int i = 0;
  while(1) {
    int left = 2 * i + 1;
    int right = left + 1;
    if(left >= heap->size) break;

    int child = left;
    if(right < heap->size && heap->data[right].dist < heap->data[left].dist) child = right;
    if(heap->data[child].dist >= last.dist) break;

    heap->data[i] = heap->data[child];
    i = child;
  }

  if(heap->size > 0) heap->data[i] = last;
  return min_item;
}

static int state_is_current(HeapItem item, int *first1, int *first2, double *dist1, double *dist2) {
  int node = item.node;
  return (first1[node] == item.first_edge && dist1[node] == item.dist) ||
         (first2[node] == item.first_edge && dist2[node] == item.dist);
}

static void swap_labels(int node, int *first1, int *first2, double *dist1, double *dist2) {
  int tmp_first = first1[node];
  double tmp_dist = dist1[node];
  first1[node] = first2[node];
  dist1[node] = dist2[node];
  first2[node] = tmp_first;
  dist2[node] = tmp_dist;
}

static void update_label(int node, int first_edge, double dist,
                         int *first1, int *first2, double *dist1, double *dist2,
                         MinHeap *heap) {
  if(first1[node] == first_edge) {
    if(dist < dist1[node]) {
      dist1[node] = dist;
      heap_push(heap, node, first_edge, dist);
    }
    return;
  }

  if(first2[node] == first_edge) {
    if(dist < dist2[node]) {
      dist2[node] = dist;
      if(dist2[node] < dist1[node]) swap_labels(node, first1, first2, dist1, dist2);
      heap_push(heap, node, first_edge, dist);
    }
    return;
  }

  if(dist < dist1[node]) {
    first2[node] = first1[node];
    dist2[node] = dist1[node];
    first1[node] = first_edge;
    dist1[node] = dist;
    heap_push(heap, node, first_edge, dist);
  } else if(dist < dist2[node]) {
    first2[node] = first_edge;
    dist2[node] = dist;
    heap_push(heap, node, first_edge, dist);
  }
}

static void add_arc(Arc *arcs, int *head, int idx, int from, int to, int edge, double cost) {
  arcs[idx].to = to;
  arcs[idx].edge = edge;
  arcs[idx].cost = cost;
  arcs[idx].next = head[from];
  head[from] = idx;
}

SEXP critical_link_detours(SEXP from, SEXP to, SEXP cost, SEXP n_nodes, SEXP directed,
                           SEXP sources, SEXP source_edges, SEXP nthreads) {
  if(TYPEOF(from) != INTSXP || TYPEOF(to) != INTSXP) error("from and to must be integer vectors");
  if(TYPEOF(cost) != REALSXP) error("cost must be a numeric vector");
  if(TYPEOF(n_nodes) != INTSXP || length(n_nodes) != 1) error("n_nodes must be a single integer");
  if(TYPEOF(directed) != LGLSXP || length(directed) != 1) error("directed must be TRUE or FALSE");
  if(TYPEOF(sources) != INTSXP) error("sources must be an integer vector");
  if(TYPEOF(source_edges) != VECSXP) error("source_edges must be a list");
  if(length(source_edges) != length(sources)) error("source_edges length must match sources length");
  if(TYPEOF(nthreads) != INTSXP || length(nthreads) != 1) error("nthreads must be a single integer");

  int m = length(from);
  if(length(to) != m || length(cost) != m) error("from, to, and cost must have the same length");

  int n = INTEGER(n_nodes)[0];
  if(n < 0) error("n_nodes must be non-negative");

  int is_directed = LOGICAL(directed)[0];
  int n_arcs = is_directed ? m : 2 * m;
  int nt = INTEGER(nthreads)[0];
  if(nt < 1) nt = 1;
  int n_sources = length(sources);

  SEXP result = PROTECT(allocVector(REALSXP, m));
  double *res = REAL(result);
  for(int i = 0; i < m; ++i) res[i] = R_PosInf;

  if(m == 0 || n == 0 || n_sources == 0) {
    UNPROTECT(1);
    return result;
  }

  int *from_ptr = INTEGER(from);
  int *to_ptr = INTEGER(to);
  double *cost_ptr = REAL(cost);
  int *src_ptr = INTEGER(sources);

  for(int i = 0; i < m; ++i) {
    if(from_ptr[i] < 1 || from_ptr[i] > n || to_ptr[i] < 1 || to_ptr[i] > n) {
      UNPROTECT(1);
      error("from and to values must be between 1 and n_nodes");
    }
  }
  for(int s = 0; s < n_sources; ++s) {
    if(src_ptr[s] < 1 || src_ptr[s] > n) {
      UNPROTECT(1);
      error("sources values must be between 1 and n_nodes");
    }
    SEXP es = VECTOR_ELT(source_edges, s);
    if(TYPEOF(es) != INTSXP) {
      UNPROTECT(1);
      error("source_edges must contain integer vectors");
    }
  }

  int *head = (int *) R_alloc((size_t) n + 1, sizeof(int));
  Arc *arcs = (Arc *) R_alloc((size_t) n_arcs, sizeof(Arc));
  for(int i = 0; i <= n; ++i) head[i] = -1;

  int arc_idx = 0;
  for(int i = 0; i < m; ++i) {
    add_arc(arcs, head, arc_idx++, from_ptr[i], to_ptr[i], i + 1, cost_ptr[i]);
    if(!is_directed) add_arc(arcs, head, arc_idx++, to_ptr[i], from_ptr[i], i + 1, cost_ptr[i]);
  }

  // Cap nthreads at the number of sources (no point launching more threads than work units)
  if(nt > n_sources) nt = n_sources;

  size_t per_thread = (size_t)(n + 1);
  double *dist1_buf = (double *) R_alloc((size_t) nt * per_thread, sizeof(double));
  double *dist2_buf = (double *) R_alloc((size_t) nt * per_thread, sizeof(double));
  int *first1_buf = (int *) R_alloc((size_t) nt * per_thread, sizeof(int));
  int *first2_buf = (int *) R_alloc((size_t) nt * per_thread, sizeof(int));

  // Heap upper bound. Each node is "validly" popped at most twice (once per label
  // slot, since labels only decrease and are popped in non-decreasing order), so
  // the main loop pushes at most 2 * n_arcs entries. The initial source expansion
  // pushes at most 2 per distinct neighbour, i.e. <= 2 * n. The maximum number of
  // simultaneously live entries can never exceed the total number of pushes, so
  // 2 * n_arcs + 2 * n + 16 is a provable upper bound. Sized once per thread so no
  // allocation is needed inside the parallel region.
  size_t heap_cap = (size_t) 2 * (size_t) n_arcs + (size_t) 2 * (size_t) n + 16;
  HeapItem *heap_buf = (HeapItem *) R_alloc((size_t) nt * heap_cap, sizeof(HeapItem));

  int overflow_flag = 0;

#ifdef _OPENMP
  #pragma omp parallel for num_threads(nt) if(nt > 1)
#endif
  for(int s = 0; s < n_sources; ++s) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    int source = src_ptr[s];

    double *dist1 = dist1_buf + (size_t) tid * per_thread;
    double *dist2 = dist2_buf + (size_t) tid * per_thread;
    int *first1 = first1_buf + (size_t) tid * per_thread;
    int *first2 = first2_buf + (size_t) tid * per_thread;

    MinHeap heap;
    heap.data = heap_buf + (size_t) tid * heap_cap;
    heap.size = 0;
    heap.capacity = (int) heap_cap;
    heap.overflow = 0;

    for(int i = 0; i <= n; ++i) {
      dist1[i] = R_PosInf;
      dist2[i] = R_PosInf;
      first1[i] = 0;
      first2[i] = 0;
    }

    for(int arc = head[source]; arc != -1; arc = arcs[arc].next) {
      int next_node = arcs[arc].to;
      if(next_node == source) continue;
      update_label(next_node, arcs[arc].edge, arcs[arc].cost,
                   first1, first2, dist1, dist2, &heap);
    }

    while(heap.size > 0) {
      HeapItem item = heap_pop(&heap);
      if(!state_is_current(item, first1, first2, dist1, dist2)) continue;

      for(int arc = head[item.node]; arc != -1; arc = arcs[arc].next) {
        int next_node = arcs[arc].to;
        if(next_node == source) continue;

        double next_dist = item.dist + arcs[arc].cost;
        update_label(next_node, item.first_edge, next_dist,
                     first1, first2, dist1, dist2, &heap);
      }
    }

    if(heap.overflow) {
      overflow_flag = 1;
      continue;
    }

    SEXP edges_for_source = VECTOR_ELT(source_edges, s);
    int n_es = length(edges_for_source);
    int *edge_ids = INTEGER(edges_for_source);
    for(int e = 0; e < n_es; ++e) {
      int edge = edge_ids[e];
      if(edge < 1 || edge > m) continue;
      int target = to_ptr[edge - 1];
      if(target == source) {
        res[edge - 1] = 0;
      } else if(first1[target] != edge) {
        res[edge - 1] = dist1[target];
      } else {
        res[edge - 1] = dist2[target];
      }
    }
  }

  if(overflow_flag) {
    UNPROTECT(1);
    error("Internal heap overflow in critical_link_detours; please file an issue");
  }

  UNPROTECT(1);
  return result;
}
