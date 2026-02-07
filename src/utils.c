#include <R.h>
#include <Rinternals.h>

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

