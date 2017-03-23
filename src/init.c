#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "normpsy.h"

static R_FortranMethodDef FortRout[] = {
  {"backtransformation", (DL_FUNC) &F77_SUB(backtransformation), 10},
  {NULL, NULL, 0}
};


void R_init_NormPsy(DllInfo * dll)
{
  R_registerRoutines(dll, NULL, NULL, FortRout, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}

