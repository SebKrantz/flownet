#include <R.h>
#include <Rinternals.h>
#include <string.h>

#ifndef SEXPPTR_RO
#define SEXPPTR_RO(x) ((const SEXP *)DATAPTR_RO(x))  // to avoid overhead of looped VECTOR_ELT
#endif

/**
 * Check if paths have duplicated edges
 *
 * @param paths1 List of integer vectors (edge numbers for first part of paths)
 * @param paths2 List of integer vectors (edge numbers for second part of paths)
 * @param delta_ks Integer vector used as hash table (must be large enough to index all edge numbers)
 * @param uniqe_edge_id Integer vector giving undirected edge ids. Needed for checking duplicates on directed graphs. NULL for undirected graphs
 * @return Integer vector with indices of the paths without duplicate edges
 */
SEXP check_path_duplicates(SEXP paths1, SEXP paths2, SEXP delta_ks, SEXP undir_edge_id) {

  int n_paths = length(paths2);
  if (length(paths1) < n_paths) {
    error("paths1 must be at least as long as paths2");
  }
  if(n_paths == 0) return allocVector(INTSXP, 0);

  // Get pointer to delta_ks for direct indexing
  int *delta_ptr = INTEGER(delta_ks)-1;

  // Allocate buffer for results
  int *buf = (int *) R_alloc(n_paths, sizeof(int)), j = 0;

  // Iterate over each path
  const SEXP *paths1_ptr = SEXPPTR_RO(paths1);
  const SEXP *paths2_ptr = SEXPPTR_RO(paths2);

  // Check inputs
  if(isNull(undir_edge_id)) {
    for (int k = 0; k < n_paths; k++) {
      int len1 = length(paths1_ptr[k]);
      int len2 = length(paths2_ptr[k]);
      double *path1_ptr = REAL(paths1_ptr[k]);
      double *path2_ptr = REAL(paths2_ptr[k]);
      int has_duplicate = 0;
      // Check edges in path1
      for (int i = 0; i < len1; i++) delta_ptr[(int)path1_ptr[i]] = 1; // Mark edge as seen
      // check path2 for duplicates with path1
      for (int i = 0; i < len2; i++) {
          if (delta_ptr[(int)path2_ptr[i]]) {
              has_duplicate = 1;
              break; // Found duplicate
          }
      }
      // Second pass: clear the hash table
      for (int i = 0; i < len1; i++) delta_ptr[(int)path1_ptr[i]] = 0;
      // Set result: TRUE if no duplicates, FALSE if duplicates
      if(!has_duplicate) buf[j++] = k+1;
    }
  } else {
    if(length(undir_edge_id) != length(delta_ks)) error("Internal length mismatch between delta_ks and undir_edge_id. Please file an issue.");
    // Get pointer to undir_edge_id for direct indexing
    int *eid_ptr = INTEGER(undir_edge_id)-1;
    for (int k = 0; k < n_paths; k++) {
      int len1 = length(paths1_ptr[k]);
      int len2 = length(paths2_ptr[k]);
      double *path1_ptr = REAL(paths1_ptr[k]);
      double *path2_ptr = REAL(paths2_ptr[k]);
      int has_duplicate = 0;
      // Check edges in path1
      for (int i = 0; i < len1; i++) delta_ptr[eid_ptr[(int)path1_ptr[i]]] = 1; // Mark edge as seen
      // check path2 for duplicates with path1
      for (int i = 0; i < len2; i++) {
        if (delta_ptr[eid_ptr[(int)path2_ptr[i]]]) {
          has_duplicate = 1;
          break; // Found duplicate
        }
      }
      // Second pass: clear the hash table
      for (int i = 0; i < len1; i++) delta_ptr[eid_ptr[(int)path1_ptr[i]]] = 0;
      // Set result: TRUE if no duplicates, FALSE if duplicates
      if(!has_duplicate) buf[j++] = k+1;
    }
  }

  SEXP result = PROTECT(allocVector(INTSXP, j));
  if(j) memcpy(INTEGER(result), buf, sizeof(int) * j);
  UNPROTECT(1);
  return result;
}


/**
 * Mark edges as traversed by incrementing counts in edges_traversed
 *
 * @param paths List of numeric/integer vectors (edge numbers for paths)
 * @param edges_traversed Integer vector to be modified in place (must be large enough to index all edge numbers)
 * @return The modified edges_traversed vector
 */
SEXP mark_edges_traversed(SEXP paths, SEXP edges_traversed) {

  int n_paths = length(paths);

  // Get pointer to edges_traversed for direct indexing
  int *edges_ptr = INTEGER(edges_traversed);

  // Iterate over each path
  const SEXP *paths_ptr = SEXPPTR_RO(paths);

  for (int k = 0; k < n_paths; k++) {
    int path_len = length(paths_ptr[k]);
    if (path_len == 0) continue; // Skip empty paths

    double *path_ptr = REAL(paths_ptr[k]);

    // Increment count for each edge in the path
    for (int i = 0; i < path_len; i++) edges_ptr[(int)path_ptr[i] - 1]++;
  }

  return edges_traversed;
}


/**
 * Set all delta_ks values for visited edges to zero
 *
 * This function resets to zero the values in the integer vector delta_ks corresponding to
 * all edges traversed by the given set of paths. For each index in no_dups, retrieves the
 * corresponding path from paths1 and paths2, and for each edge in these paths, sets its delta_ks
 * entry to zero. Finally, for each edge in shortest_path, the delta_ks entry is also set to zero.
 *
 * Intended for use after edge count tallies (delta_ks) are no longer needed for those paths.
 *
 * @param delta_ks    Integer vector to be zeroed in-place for traversed edges
 * @param no_dups     Integer vector of indices (1-based) of non-duplicate paths
 * @param paths1      List of vectors; primary part of alternative paths (double vectors of edge IDs)
 * @param paths2      List of vectors; secondary part of alternative paths (double vectors of edge IDs)
 * @param shortest_path Numeric vector of edge IDs for the shortest path
 * @return            The modified delta_ks vector (as SEXP)
SEXP free_delta_ks(SEXP delta_ks, SEXP no_dups, SEXP paths1, SEXP paths2, SEXP shortest_path) {
  int n_no_dups = length(no_dups);
  const SEXP *paths1_ptr = SEXPPTR_RO(paths1);
  const SEXP *paths2_ptr = SEXPPTR_RO(paths2);
  int *no_dups_ptr = INTEGER(no_dups);
  double *shortest_path_ptr = REAL(shortest_path);
  int shortest_path_len = length(shortest_path);
  int *delta_ptr = INTEGER(delta_ks)-1;

  // Zero out delta_ks entries for all edges in paths1 and paths2 for non-duplicate paths
  for (int idx = 0; idx < n_no_dups; idx++) {
    int k = no_dups_ptr[idx] - 1;
    int len1 = length(paths1_ptr[k]);
    int len2 = length(paths2_ptr[k]);
    double *p1 = REAL(paths1_ptr[k]);
    double *p2 = REAL(paths2_ptr[k]);
    for (int i = 0; i < len1; i++) delta_ptr[(int)p1[i]] = 0;
    for (int i = 0; i < len2; i++) delta_ptr[(int)p2[i]] = 0;
  }

  // Zero out delta_ks entries for all edges in the shortest path
  for (int i = 0; i < shortest_path_len; i++) delta_ptr[(int)shortest_path_ptr[i]] = 0;

  return delta_ks;
}
 */

// Assign to list inside mirai daemon
SEXP set_vector_elt(SEXP x, SEXP i, SEXP elt) {
  int idx = asInteger(i) - 1;
  if(TYPEOF(x) == INTSXP) INTEGER(x)[idx] = INTEGER(elt)[0];
  else if(TYPEOF(x) == REALSXP) REAL(x)[idx] = asReal(elt);
  else SET_VECTOR_ELT(x, idx, elt);
  return R_NilValue;
}


/**
 * Assign flows to edges for multiple paths (batch AoN assignment)
 *
 * @param paths List of numeric vectors (edge indices for each path)
 * @param flows Numeric vector of flow values (one per path)
 * @param final_flows Numeric vector to accumulate flows (modified in place)
 * @param indices Integer indices of od-pairs (to) processed
 * @param od_pairs Integer vector indicating whether OD pair is valid
 * @return The modified final_flows vector
 */
SEXP assign_flows_to_paths(SEXP paths, SEXP flows, SEXP final_flows, SEXP indices, SEXP od_pairs) {
  int n_paths = length(paths);
  int n_flows = length(indices);
  if (n_paths != n_flows) {
    error("paths and flows must have the same length");
  }
  double *flows_vals = REAL(flows);
  double *final_ptr = REAL(final_flows);
  const SEXP *paths_ptr = SEXPPTR_RO(paths);
  int *idx = INTEGER(indices);
  int *odp = INTEGER(od_pairs);

  for (int k = 0; k < n_paths; k++) {
    int path_len = length(paths_ptr[k]);
    if (path_len == 0) {
      odp[idx[k]-1] = NA_INTEGER;
      continue;
    }

    double flow_val = flows_vals[idx[k]-1];
    double *path_ptr = REAL(paths_ptr[k]);

    for (int i = 0; i < path_len; i++) {
      final_ptr[(int)path_ptr[i] - 1] += flow_val;
    }
  }

  return final_flows;
}


/**
 * Compute sum of costs for each path and assign directly to indexed positions
 *
 * @param paths List of numeric vectors (edge indices for each path)
 * @param cost Numeric vector of edge costs
 * @param result Numeric vector to store results (modified in place)
 * @param indices Integer vector of 1-based indices into result where costs should be stored
 * @return The modified result vector
 */
SEXP sum_path_costs(SEXP paths, SEXP cost, SEXP result, SEXP indices) {
  int n_paths = length(paths);
  int n_indices = length(indices);
  if (n_paths != n_indices) {
    error("paths and indices must have the same length");
  }
  double *cost_ptr = REAL(cost);
  double *result_ptr = REAL(result);
  int *idx_ptr = INTEGER(indices);
  const SEXP *paths_ptr = SEXPPTR_RO(paths);

  for (int k = 0; k < n_paths; k++) {
    int path_len = length(paths_ptr[k]);
    int result_idx = idx_ptr[k] - 1;  // Convert to 0-based

    if (path_len == 0) {
      result_ptr[result_idx] = NA_REAL;
      continue;
    }

    double sum = 0.0;
    double *path_ptr = REAL(paths_ptr[k]);

    for (int i = 0; i < path_len; i++) {
      sum += cost_ptr[(int)path_ptr[i] - 1];
    }
    result_ptr[result_idx] = sum;
  }

  return result;
}

/**
 * Contract already-oriented intermediate nodes in an undirected graph.
 *
 * @param from Integer vector of origin node ids
 * @param to Integer vector of destination node ids
 * @param gid Integer vector mapping local edges to group ids
 * @param nodes Integer vector of candidate intermediate node ids
 * @return list(from = <int>, to = <int>, gid = <int>, ok = <logical>)
 */
SEXP contract_linear_nodes(SEXP from, SEXP to, SEXP gid, SEXP nodes) {
  if (TYPEOF(from) != INTSXP || TYPEOF(to) != INTSXP || TYPEOF(gid) != INTSXP || TYPEOF(nodes) != INTSXP) {
    error("C_contract_linear_nodes expects integer vectors for from, to, gid, and nodes.");
  }
  R_xlen_t n_edges = XLENGTH(from);
  R_xlen_t n_nodes = XLENGTH(nodes);
  if (XLENGTH(to) != n_edges || XLENGTH(gid) != n_edges) {
    error("Internal length mismatch: from, to, and gid must have identical length.");
  }

  SEXP out_from = PROTECT(duplicate(from));
  SEXP out_to = PROTECT(duplicate(to));
  SEXP out_gid = PROTECT(duplicate(gid));

  int *from_ptr = INTEGER(out_from);
  int *to_ptr = INTEGER(out_to);
  int *gid_ptr = INTEGER(out_gid);
  int *nodes_ptr = INTEGER(nodes);

  /*
   * `consolidate_graph_core()` remaps endpoints to dense 1..N (funique + fmatch)
   * before calling here, so node_slot[] is O(N) and avoids hashing. We still
   * scan for max_node to size that array and remain safe if inputs are not
   * remapped in a future code path.
   */
  int max_node = 0;
  for (R_xlen_t e = 0; e < n_edges; e++) {
    int v = from_ptr[e];
    if (v != NA_INTEGER && v > max_node) max_node = v;
    v = to_ptr[e];
    if (v != NA_INTEGER && v > max_node) max_node = v;
  }
  for (R_xlen_t i = 0; i < n_nodes; i++) {
    int v = nodes_ptr[i];
    if (v != NA_INTEGER && v > max_node) max_node = v;
  }
  if (max_node <= 0) {
    SEXP ok = PROTECT(ScalarLogical(FALSE));
    SEXP ans = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, out_from);
    SET_VECTOR_ELT(ans, 1, out_to);
    SET_VECTOR_ELT(ans, 2, out_gid);
    SET_VECTOR_ELT(ans, 3, ok);
    SEXP nms = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(nms, 0, mkChar("from"));
    SET_STRING_ELT(nms, 1, mkChar("to"));
    SET_STRING_ELT(nms, 2, mkChar("gid"));
    SET_STRING_ELT(nms, 3, mkChar("ok"));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(6);
    return ans;
  }

  /* node_slot[v] = i+1 if v is the (i+1)-th candidate node, else 0. */
  size_t dense_len = (size_t) max_node + 1U;
  int *node_slot = (int *) R_alloc(dense_len, sizeof(int));
  memset(node_slot, 0, dense_len * sizeof(int));
  for (R_xlen_t i = 0; i < n_nodes; i++) {
    int v = nodes_ptr[i];
    if (v == NA_INTEGER || v <= 0) continue;
    if (node_slot[v] == 0) node_slot[v] = (int) (i + 1);
  }

  /*
   * Precompute, for each candidate intermediate node (slot i in `nodes`):
   *   in_edge[i]  = smallest edge index where gft$to   == nodes[i]
   *   out_edge[i] = smallest edge index where gft$from == nodes[i]
   * Matches collapse::fmatch semantics (first match) used by the R path.
   */
  int *in_edge = (int *) R_alloc((size_t) n_nodes, sizeof(int));
  int *out_edge = (int *) R_alloc((size_t) n_nodes, sizeof(int));
  char *processed = (char *) R_alloc((size_t) n_nodes, sizeof(char));
  memset(in_edge, 0, (size_t) n_nodes * sizeof(int));
  memset(out_edge, 0, (size_t) n_nodes * sizeof(int));
  memset(processed, 0, (size_t) n_nodes * sizeof(char));

  for (R_xlen_t e = n_edges; e > 0; e--) {
    R_xlen_t idx = e - 1;
    int f = from_ptr[idx];
    int t = to_ptr[idx];
    if (f != NA_INTEGER && f > 0) {
      int pos = node_slot[f];
      if (pos) out_edge[pos - 1] = (int) e;
    }
    if (t != NA_INTEGER && t > 0) {
      int pos = node_slot[t];
      if (pos) in_edge[pos - 1] = (int) e;
    }
  }

  /*
   * Walk each linear chain exactly once. A node i is a chain head iff the
   * `from` endpoint of its in-edge is not itself an intermediate (i.e. not a
   * candidate node with both an in- and an out-edge). Each walk:
   *   - marks the outgoing edge `c_oe` as merged (from=NA, gid=gid[head_edge])
   *   - advances the head's `to` to the next node
   *   - records the in-edge of the intermediate for a final `to` rewrite
   *     so that every merged-into edge shares the same (from, to) pair
   *     after ffirst()-fill, matching the R implementation's semantics.
   * Pure cycles of intermediates have no head and are safely skipped.
   */
  int *chain_buf = (int *) R_alloc((size_t) n_nodes, sizeof(int));
  int ok_flag = 0;
  for (R_xlen_t i = 0; i < n_nodes; i++) {
    if (processed[i]) continue;
    int ie = in_edge[i];
    int oe = out_edge[i];
    if (ie == 0 || oe == 0) {
      processed[i] = 1;
      continue;
    }

    int pred = from_ptr[ie - 1];
    int pred_pos = (pred == NA_INTEGER || pred <= 0) ? 0 : node_slot[pred];
    if (pred_pos != 0 && in_edge[pred_pos - 1] != 0 && out_edge[pred_pos - 1] != 0) {
      continue; // not a chain head; will be reached by its predecessor
    }

    int head_edge = ie;
    int head_gid = gid_ptr[head_edge - 1];
    int chain_len = 0;
    R_xlen_t cur_pos = i;
    int end_node = NA_INTEGER;
    while (1) {
      int c_oe = out_edge[cur_pos];
      if (c_oe == 0) break;
      int next_node = to_ptr[c_oe - 1];
      gid_ptr[c_oe - 1] = head_gid;
      from_ptr[c_oe - 1] = NA_INTEGER;
      processed[cur_pos] = 1;
      chain_buf[chain_len++] = (int) in_edge[cur_pos];
      ok_flag = 1;
      end_node = next_node;

      if (next_node == NA_INTEGER || next_node <= 0) break;
      int next_pos = node_slot[next_node];
      if (next_pos == 0) break;
      if (processed[next_pos - 1]) break;
      if (in_edge[next_pos - 1] == 0 || out_edge[next_pos - 1] == 0) break;
      cur_pos = (R_xlen_t) (next_pos - 1);
    }

    /* Align every merged in-edge with the final chain endpoint so that
     * downstream GRP() groups these edges with the representative. */
    for (int j = 0; j < chain_len; j++) {
      to_ptr[chain_buf[j] - 1] = end_node;
    }
  }

  SEXP ok = PROTECT(ScalarLogical(ok_flag));
  SEXP ans = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(ans, 0, out_from);
  SET_VECTOR_ELT(ans, 1, out_to);
  SET_VECTOR_ELT(ans, 2, out_gid);
  SET_VECTOR_ELT(ans, 3, ok);
  SEXP nms = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(nms, 0, mkChar("from"));
  SET_STRING_ELT(nms, 1, mkChar("to"));
  SET_STRING_ELT(nms, 2, mkChar("gid"));
  SET_STRING_ELT(nms, 3, mkChar("ok"));
  setAttrib(ans, R_NamesSymbol, nms);

  UNPROTECT(6);
  return ans;
}

