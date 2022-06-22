#include <stdlib.h>
#include "R-hdpMultinomial_iterate.h"

SEXP hdpMultinomial_iterate(SEXP hdpin, SEXP numiter, SEXP doconparam, SEXP dolik, SEXP numallocs, SEXP dodebug)
{
  int ni, docp, dl, nallocs;

  GetRNGstate();

  ni = asInteger(numiter);
  docp = asInteger(doconparam);
  dl = asInteger(dolik);
  nallocs = asInteger(numallocs);

  DEBUG = asInteger(dodebug);

  HDP *hdp = rReadHDP(hdpin);

  SEXP LIK = PROTECT(allocVector(REALSXP, ni));

  SEXP ALLOCATIONS;

  if (nallocs > 0) {
    ALLOCATIONS = PROTECT(allocVector(INTSXP, nallocs));
  } else {
    ALLOCATIONS = R_NilValue;
  }
  rdebug0(1,"Running hdpMultinomial_iterate.\n");
  if (nallocs > 0) {
    hdp_iterate(hdp, REAL(LIK), ni, docp, dl, INTEGER(ALLOCATIONS));
  } else {
    hdp_iterate(hdp, REAL(LIK), ni, docp, dl, NULL);
  }
  rdebug0(1,"Finished hdpMultinomial_iterate.\n");

  SEXP hdpout = PROTECT(duplicate(hdpin));
  rWriteHDP(hdpout,hdp);

  SEXP result = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(result, 0, hdpout);
  SET_VECTOR_ELT(result, 1, LIK);
  SET_VECTOR_ELT(result, 2, ALLOCATIONS);

  if (nallocs > 0) {
    UNPROTECT(4);
  } else {
    UNPROTECT(3);
  }

  PutRNGstate();

  return result;
}


