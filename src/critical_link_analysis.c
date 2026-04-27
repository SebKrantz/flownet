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
} MinHeap;

static void heap_init(MinHeap *heap, int capacity) {
  heap->size = 0;
  heap->capacity = capacity > 0 ? capacity : 1;
  heap->data = (HeapItem *) R_Calloc(heap->capacity, HeapItem);
}

static void heap_free(MinHeap *heap) {
  if(heap->data != NULL) R_Free(heap->data);
  heap->data = NULL;
  heap->size = 0;
  heap->capacity = 0;
}

static void heap_push(MinHeap *heap, int node, int first_edge, double dist) {
  if(heap->size >= heap->capacity) {
    heap->capacity *= 2;
    heap->data = (HeapItem *) R_Realloc(heap->data, heap->capacity, HeapItem);
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

SEXP critical_link_detours(SEXP from, SEXP to, SEXP cost, SEXP n_nodes, SEXP directed) {
  if(TYPEOF(from) != INTSXP || TYPEOF(to) != INTSXP) error("from and to must be integer vectors");
  if(TYPEOF(cost) != REALSXP) error("cost must be a numeric vector");
  if(TYPEOF(n_nodes) != INTSXP || length(n_nodes) != 1) error("n_nodes must be a single integer");
  if(TYPEOF(directed) != LGLSXP || length(directed) != 1) error("directed must be TRUE or FALSE");

  int m = length(from);
  if(length(to) != m || length(cost) != m) error("from, to, and cost must have the same length");

  int n = INTEGER(n_nodes)[0];
  if(n < 0) error("n_nodes must be non-negative");

  int is_directed = LOGICAL(directed)[0];
  int n_arcs = is_directed ? m : 2 * m;

  SEXP result = PROTECT(allocVector(REALSXP, m));
  double *res = REAL(result);
  for(int i = 0; i < m; ++i) res[i] = R_PosInf;

  if(m == 0 || n == 0) {
    UNPROTECT(1);
    return result;
  }

  int *from_ptr = INTEGER(from);
  int *to_ptr = INTEGER(to);
  double *cost_ptr = REAL(cost);

  int *head = (int *) R_alloc((size_t) n + 1, sizeof(int));
  Arc *arcs = (Arc *) R_alloc((size_t) n_arcs, sizeof(Arc));
  double *dist1 = (double *) R_alloc((size_t) n + 1, sizeof(double));
  double *dist2 = (double *) R_alloc((size_t) n + 1, sizeof(double));
  int *first1 = (int *) R_alloc((size_t) n + 1, sizeof(int));
  int *first2 = (int *) R_alloc((size_t) n + 1, sizeof(int));

  for(int i = 0; i <= n; ++i) head[i] = -1;

  int arc_idx = 0;
  for(int i = 0; i < m; ++i) {
    if(from_ptr[i] < 1 || from_ptr[i] > n || to_ptr[i] < 1 || to_ptr[i] > n) {
      error("from and to values must be between 1 and n_nodes");
    }
    add_arc(arcs, head, arc_idx++, from_ptr[i], to_ptr[i], i + 1, cost_ptr[i]);
    if(!is_directed) add_arc(arcs, head, arc_idx++, to_ptr[i], from_ptr[i], i + 1, cost_ptr[i]);
  }

  MinHeap heap;
  heap_init(&heap, n_arcs > 1024 ? n_arcs : 1024);

  for(int source = 1; source <= n; ++source) {
    R_CheckUserInterrupt();

    for(int i = 0; i <= n; ++i) {
      dist1[i] = R_PosInf;
      dist2[i] = R_PosInf;
      first1[i] = 0;
      first2[i] = 0;
    }
    heap.size = 0;

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

    for(int edge = 0; edge < m; ++edge) {
      if(from_ptr[edge] != source) continue;

      int target = to_ptr[edge];
      if(target == source) {
        res[edge] = 0;
      } else if(first1[target] != edge + 1) {
        res[edge] = dist1[target];
      } else {
        res[edge] = dist2[target];
      }
    }
  }

  heap_free(&heap);
  UNPROTECT(1);
  return result;
}
