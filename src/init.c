#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP check_path_duplicates(SEXP paths1, SEXP paths2, SEXP delta_ks);
SEXP compute_path_sized_logit(SEXP paths1, SEXP paths2, SEXP no_dups, SEXP shortest_path,
                               SEXP cost, SEXP cost_ks, SEXP d_ij, SEXP beta_PSL, SEXP flow,
                               SEXP delta_ks, SEXP final_flows);
SEXP mark_edges_traversed(SEXP paths, SEXP edges_traversed);

static const R_CallMethodDef CallEntries[] = {
  {"C_check_path_duplicates", (DL_FUNC) &check_path_duplicates, 3},
  {"C_compute_path_sized_logit", (DL_FUNC) &compute_path_sized_logit, 11},
  {"C_mark_edges_traversed", (DL_FUNC) &mark_edges_traversed, 2},
  {NULL, NULL, 0}
};

void R_init_mmflowr(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
