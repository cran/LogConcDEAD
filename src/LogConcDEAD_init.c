#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void logconestw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP convhull(SEXP, SEXP);
extern SEXP delaunayn(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"logconestw", (DL_FUNC) &logconestw, 10},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"convhull", (DL_FUNC) &convhull, 2},
    {"delaunayn",   (DL_FUNC) &delaunayn,   2},
    {NULL, NULL, 0}
};

void R_init_LogConcDEAD(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

