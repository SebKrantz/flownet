#include <R.h>
#include <Rinternals.h>

#ifndef SEXPPTR_RO
#define SEXPPTR_RO(x) ((const SEXP *)DATAPTR_RO(x))  // to avoid overhead of looped VECTOR_ELT
#endif

static double POS_INF = 1.0/0.0;

/**
 * compute_path_sized_logit
 *
 * Compute path-size logit probabilities for a set of alternative and shortest paths,
 * and attributes flows to network edges for each alternative. Used in transit assignment
 * and route choice modeling with explicit path-overlap correction via the path-size logit model.
 *
 * Parameters
 * ----------
 * paths1 : SEXP (list of numeric vectors)
 *      The first segment of the path. List of alternatives, each a numeric vector (1-based edge indices).
 * paths2 : SEXP (list of numeric vectors)
 *      The second segment of the path. List paralleling paths1, each a numeric vector (1-based edge indices).
 * no_dups : SEXP (integer vector)
 *      Indices (1-based) of paths to consider (those with no duplicate edges).
 * shortest_path : SEXP (numeric vector)
 *      Edge ids (1-based) forming the shortest path for this OD pair.
 * cost : SEXP (numeric vector)
 *      Cost of each edge in the network.
 * cost_ks : SEXP (numeric vector)
 *      Total cost for each alternative path with no duplicate edges.
 * d_ij : SEXP (length 1, double)
 *      Cost of the shortest path for this OD pair.
 * beta_PSL : SEXP (length 1, double)
 *      Path-size logit parameter. Determines strength of path-size correction.
 * flow : SEXP (length 1, double)
 *      OD matrix flow value (amount of travelers/objects).
 * delta_ks : SEXP (integer vector)
 *      An integer vector, same length as number of edges, used as temporary workspace for edge overlap counts.
 * edge_probs : SEXP (double vector)
 *      A double vector, same length as number of edges, used as temporary workspace for edge probability accumulation.
 *      Only used when retvals_PSL[2] is TRUE. Reused across OD pairs (reset only for visited edges).
 * final_flows : SEXP (double vector)
 *      A vector, same length as number of edges, accumulates the assigned flow to each edge across all paths.
 * retvals_PSL : SEXP (logical vector of length 3)
 *      [0]: return edges vector?; [1] return edge counts?; [2] return edge weights?
 *
 * Returns
 * -------
 * If !any(retvals_PSL):
 *   prob_ks : SEXP (double vector of length (n_no_dups + 1))
 *      Probability vector. The first element corresponds to the shortest path, and the next n_no_dups elements
 *      correspond to the alternative paths given by no_dups.
 *
 * If any(retvals_PSL):
 *   A list with up to four elements:
 *     [[1]] prob_ks: as above
 *     [[2]] edges: edge indices (1-based)
 *     [[3]] counts: edge counts
 *     [[4]] eweights: vector with summed probabilities per edge (compact form)
 *
 * Details
 * -------
 * Computes path overlap factors (gamma) for alternatives and the shortest path.
 * Calculates logit probabilities with path-size correction. Updates a temporary vector delta_ks to count the overlap
 * of network edges among candidate paths. Assigns OD flow to edges in final_flows by the chosen probabilities.
 *
 * Edge numbering: All edge ids in the paths and shortest_path are assumed to be 1-based double values.
 */
SEXP compute_path_sized_logit(SEXP paths1, SEXP paths2, SEXP no_dups, SEXP shortest_path,
                              SEXP cost, SEXP cost_ks, SEXP d_ij, SEXP beta_PSL,
                              SEXP flow, SEXP delta_ks, SEXP edge_probs,
                              SEXP final_flows, SEXP retvals_PSL) {

  int n_no_dups = length(no_dups);
  const SEXP *paths1_ptr = SEXPPTR_RO(paths1);
  const SEXP *paths2_ptr = SEXPPTR_RO(paths2);
  int *no_dups_ptr = INTEGER(no_dups);
  double *shortest_path_ptr = REAL(shortest_path);
  int shortest_path_len = length(shortest_path);
  double *cost_ptr = REAL(cost);
  double *cost_ks_ptr = REAL(cost_ks);
  double d_ij_val = REAL(d_ij)[0];
  double beta_PSL_val = asReal(beta_PSL);
  double flow_val = asReal(flow);
  int *delta_ptr = INTEGER(delta_ks)-1; // offset for 1-based edge IDs
  double *edge_probs_ptr = REAL(edge_probs)-1; // same offset for 1-based indexing
  double *final_flows_ptr = REAL(final_flows);
  if(length(retvals_PSL) != 4) error("Internal error: retvals_PSL needs to be of length 4. Please file an issue.");
  int *ret = LOGICAL(retvals_PSL), ret_any = ret[0] + ret[1] + ret[2]; // no need to know number of edges for ret[3]


  // Step 1: Update delta_ks for paths in no_dups
  // delta_ks stores, for each edge, how many times it appears in the alternatives under consideration
  int n_edges = 0;
  if(ret_any) {
    for (int idx = 0; idx < n_no_dups; idx++) {
      int k = no_dups_ptr[idx] - 1; // Convert to 0-based
      int len1 = length(paths1_ptr[k]);
      int len2 = length(paths2_ptr[k]);
      double *p1 = REAL(paths1_ptr[k]);
      double *p2 = REAL(paths2_ptr[k]);
      for (int i = 0; i < len1; i++) {
        int ip = (int)p1[i];
        if(delta_ptr[ip] == 0) ++n_edges;
        delta_ptr[ip]++;
      }
      for (int i = 0; i < len2; i++) {
        int ip = (int)p2[i];
        if(delta_ptr[ip] == 0) ++n_edges;
        delta_ptr[ip]++;
      }
    }
    // Update delta_ks for shortest_path
    for (int i = 0; i < shortest_path_len; i++) {
      int ip = (int)shortest_path_ptr[i];
      if(delta_ptr[ip] == 0) ++n_edges;
      delta_ptr[ip]++;
    }
  } else {
    for (int idx = 0; idx < n_no_dups; idx++) {
      int k = no_dups_ptr[idx] - 1; // Convert to 0-based
      int len1 = length(paths1_ptr[k]);
      int len2 = length(paths2_ptr[k]);
      double *p1 = REAL(paths1_ptr[k]);
      double *p2 = REAL(paths2_ptr[k]);
      for (int i = 0; i < len1; i++) delta_ptr[(int)p1[i]]++;
      for (int i = 0; i < len2; i++) delta_ptr[(int)p2[i]]++;
    }
    // Update delta_ks for shortest_path
    for (int i = 0; i < shortest_path_len; i++) delta_ptr[(int)shortest_path_ptr[i]]++;
  }

  // Step 2: Compute gamma_ks and gamma_1 (path-size factors for each alternative and the shortest path)

  double *gamma_ks = NULL;
  SEXP PS = R_NilValue;
  if(ret[3]) {
    PROTECT(PS = allocVector(REALSXP, n_no_dups + 1));
    gamma_ks = REAL(PS)+1; // offset because first index needs to be shortest path
  } else {
    gamma_ks = (double *) R_alloc(n_no_dups, sizeof(double));
  }
  double gamma_1 = 0.0;

  for (int idx = 0; idx < n_no_dups; idx++) {
    int k = no_dups_ptr[idx] - 1;
    int len1 = length(paths1_ptr[k]);
    int len2 = length(paths2_ptr[k]);
    double *p1 = REAL(paths1_ptr[k]);
    double *p2 = REAL(paths2_ptr[k]);
    double sum = 0.0;
    for (int i = 0; i < len1; i++) {
      int edge = (int)p1[i];
      sum += cost_ptr[edge-1] / delta_ptr[edge];
    }
    for (int i = 0; i < len2; i++) {
      int edge = (int)p2[i];
      sum += cost_ptr[edge-1] / delta_ptr[edge];
    }
    gamma_ks[idx] = sum / cost_ks_ptr[k];
    // Ensure gamma is positive to avoid log(0) or log(negative)
    if (gamma_ks[idx] <= 0.0) gamma_ks[idx] = 1e-10;
  }

  // Compute gamma_1 for shortest_path
  for (int i = 0; i < shortest_path_len; i++) {
    int edge = (int)shortest_path_ptr[i];
    gamma_1 += cost_ptr[edge-1] / delta_ptr[edge];
  }
  gamma_1 /= d_ij_val;
  // Ensure gamma_1 is positive to avoid log(0) or log(negative)
  if (gamma_1 <= 0.0) gamma_1 = 1e-10;
  if(ret[3]) gamma_ks[0-1] = gamma_1; // store shortest path gamma in first position if requested for output

  // Step 3: Compute prob_ks using log-sum-exp for numerical stability
  SEXP prob_ks = PROTECT(allocVector(REALSXP, n_no_dups + 1));
  double *prob_ptr = REAL(prob_ks);
  memset(prob_ptr, 0, (n_no_dups + 1) * sizeof(double));

  // First pass: compute max utility for numerical stability
  double max_utility = -POS_INF;
  for (int idx = 0; idx < n_no_dups; idx++) {
    double util = -cost_ks_ptr[no_dups_ptr[idx]-1] + beta_PSL_val * log(gamma_ks[idx]);
    if(util > max_utility && util != -POS_INF) max_utility = util;
  }
  double util_shortest = -d_ij_val + beta_PSL_val * log(gamma_1);
  if(util_shortest > max_utility && util_shortest != -POS_INF) max_utility = util_shortest;

  double sum_exp = 0.0;
  if(max_utility != -POS_INF) {
    // Second pass: compute exponentials with max normalization to prevent underflow
    for (int idx = 0; idx < n_no_dups; idx++) {
      double util = -cost_ks_ptr[no_dups_ptr[idx]-1] + beta_PSL_val * log(gamma_ks[idx]);
      double util_diff = util - max_utility;
      // Clamp to reasonable range to prevent overflow/underflow
      if(util_diff < -700.0) util_diff = -700.0;  // exp(-700) ≈ 0, but won't underflow
      prob_ptr[idx+1] = exp(util_diff);
      sum_exp += prob_ptr[idx+1];
    }
    double util_diff_0 = util_shortest - max_utility;
    if(util_diff_0 < -700.0) util_diff_0 = -700.0;
    prob_ptr[0] = exp(util_diff_0);
    sum_exp += prob_ptr[0];
  }

  // Handle case where all utilities are -infinity
  // More lenient check: allow very small but positive sum_exp (can happen with extreme costs)
  // Only fail if truly invalid (NaN, <= 0, or infinity)
  if(max_utility == -POS_INF || sum_exp != sum_exp || sum_exp <= 1e-300 || sum_exp >= POS_INF) {
    // Reset delta_ks before returning
    for (int idx = 0; idx < n_no_dups; idx++) {
      int k = no_dups_ptr[idx] - 1;
      int len1 = length(paths1_ptr[k]);
      int len2 = length(paths2_ptr[k]);
      double *p1 = REAL(paths1_ptr[k]);
      double *p2 = REAL(paths2_ptr[k]);
      for (int i = 0; i < len1; i++) delta_ptr[(int)p1[i]] = 0;
      for (int i = 0; i < len2; i++) delta_ptr[(int)p2[i]] = 0;
    }
    for (int i = 0; i < shortest_path_len; i++) delta_ptr[(int)shortest_path_ptr[i]] = 0;
    UNPROTECT(1);
    return R_NilValue;
  }

  // Normalize to get probabilities
  for (int i = 0; i <= n_no_dups; i++) prob_ptr[i] /= sum_exp;

  // Step 4: Update final_flows (and edge_probs if requested)
  if (ret_any || ret[3]) {
    // Return list with path weights and edges, edge_counts, and edge_weights
    SEXP result = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(result, 0, prob_ks);
    if(ret[3]) SET_VECTOR_ELT(result, 4, PS); // include path-size factors in output if requested
    int *pe = NULL, *pec = NULL,
      k = 0, l = length(delta_ks), lp = l+1,
      ret0 = ret[0], ret1 = ret[1], ret2 = ret[2];
    // Allocate edges array if requested OR if we need it for edge weights extraction
    if(ret0) {
      SET_VECTOR_ELT(result, 1, allocVector(INTSXP, n_edges));
      pe = INTEGER(VECTOR_ELT(result, 1));
    }
    if(ret1) {
      SET_VECTOR_ELT(result, 2, allocVector(INTSXP, n_edges));
      pec = INTEGER(VECTOR_ELT(result, 2));
    }
    // Accumulate final_flows (and edge_probs if requested)
    for (int idx = 0; idx < n_no_dups; idx++) {
      double prob_val = prob_ptr[idx+1];
      if(prob_val != prob_val || prob_val <= 0) continue;
      double flow_prob = flow_val * prob_val;
      int ki = no_dups_ptr[idx] - 1;
      int len1 = length(paths1_ptr[ki]);
      int len2 = length(paths2_ptr[ki]);
      double *p1 = REAL(paths1_ptr[ki]);
      double *p2 = REAL(paths2_ptr[ki]);
      for (int i = 0; i < len1; i++) {
        int edge = (int)p1[i];
        if(ret2) edge_probs_ptr[edge] += prob_val;
        final_flows_ptr[edge-1] += flow_prob;
      }
      for (int i = 0; i < len2; i++) {
        int edge = (int)p2[i];
        if(ret2) edge_probs_ptr[edge] += prob_val;
        final_flows_ptr[edge-1] += flow_prob;
      }
    }
    // Handle shortest path
    double prob_shortest = prob_ptr[0];
    if(prob_shortest == prob_shortest && prob_shortest > 0) {
      double flow_prob_shortest = flow_val * prob_shortest;
      for (int i = 0; i < shortest_path_len; i++) {
        int edge = (int)shortest_path_ptr[i];
        if(ret2) edge_probs_ptr[edge] += prob_shortest;
        final_flows_ptr[edge-1] += flow_prob_shortest;
      }
    }
    // Loop through delta_ptr to extract edges/counts/eweights and reset
    if(ret2) {
      SET_VECTOR_ELT(result, 3, allocVector(REALSXP, n_edges));
      double *pew = REAL(VECTOR_ELT(result, 3));
      for(int i = 1; i != lp; ++i) {
        if(delta_ptr[i]) {
          if(ret1) pe[k] = i; // Store edge index if needed
          if(ret2) pec[k] = delta_ptr[i]; // Store edge count if requested
          pew[k] = edge_probs_ptr[i]; // Extract compact edge weight
          k++;
          delta_ptr[i] = 0; // Reset delta_ks for next OD pair
          edge_probs_ptr[i] = 0.0; // Reset edge_probs for next OD pair
        }
      }
    } else {
      for(int i = 1; i != lp; ++i) {
        if(delta_ptr[i]) {
          if(ret1) pe[k] = i; // Store edge index if needed
          if(ret2) pec[k] = delta_ptr[i]; // Store edge count if requested
          k++;
          delta_ptr[i] = 0; // Reset delta_ks for next OD pair
        }
      }
    }
    UNPROTECT(2+ret[3]);
    return result;
  }

  // Only update final_flows (original code path - no overhead)
  for (int idx = 0; idx < n_no_dups; idx++) {
    double prob_val = flow_val * prob_ptr[idx+1];
    if(prob_val != prob_val || prob_val <= 0) continue;
    int k = no_dups_ptr[idx] - 1;
    int len1 = length(paths1_ptr[k]);
    int len2 = length(paths2_ptr[k]);
    double *p1 = REAL(paths1_ptr[k]);
    double *p2 = REAL(paths2_ptr[k]);
    for (int i = 0; i < len1; i++) {
      int edge = (int)p1[i];
      final_flows_ptr[edge-1] += prob_val;
      delta_ptr[edge] = 0;
    }
    for (int i = 0; i < len2; i++) {
      int edge = (int)p2[i];
      final_flows_ptr[edge-1] += prob_val;
      delta_ptr[edge] = 0;
    }
  }
  // Handle shortest path
  double prob_shortest = flow_val * prob_ptr[0];
  if(prob_shortest == prob_shortest && prob_shortest > 0) {
    for (int i = 0; i < shortest_path_len; i++) {
      int edge = (int)shortest_path_ptr[i];
      final_flows_ptr[edge-1] += prob_shortest;
      delta_ptr[edge] = 0;
    }
  }

  UNPROTECT(1);
  return prob_ks;
}
