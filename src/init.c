#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP C_cpp_gmm_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_cpp_gmm_fit_many(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"C_cpp_gmm_fit", (DL_FUNC) &C_cpp_gmm_fit, 7},
  {"C_cpp_gmm_fit_many", (DL_FUNC) &C_cpp_gmm_fit_many, 8},
  {NULL, NULL, 0}
};

void R_init_CellCharterR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
