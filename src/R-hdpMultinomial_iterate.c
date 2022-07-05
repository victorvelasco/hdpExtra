#include <stdlib.h>
#include "R-hdpMultinomial_iterate.h"

SEXP hdpMultinomial_iterate(SEXP hdpin, SEXP numiter, SEXP doconparam, SEXP dolik, SEXP numallocs, SEXP numcp, SEXP dodebug)
{
  int ni, docp, dl, nallocs, ncp;

  GetRNGstate();

  ni = asInteger(numiter);
  docp = asInteger(doconparam);
  dl = asInteger(dolik);
  nallocs = asInteger(numallocs);
  ncp = asInteger(numcp);

  DEBUG = asInteger(dodebug);

  HDP *hdp = rReadHDP(hdpin);

  SEXP LIK = PROTECT(allocVector(REALSXP, ni));
  SEXP CONPARAMETER = PROTECT(allocVector(REALSXP, ni*ncp));

  SEXP ALLOCATIONS;

  if (nallocs > 0) {
    ALLOCATIONS = PROTECT(allocVector(INTSXP, nallocs));
  } else {
    ALLOCATIONS = R_NilValue;
  }
  rdebug0(1,"Running hdpMultinomial_iterate.\n");
  if (nallocs > 0) {
    hdp_iterate(hdp, REAL(LIK), ni, ncp, docp, dl, INTEGER(ALLOCATIONS), REAL(CONPARAMETER));
  } else {
    hdp_iterate(hdp, REAL(LIK), ni, ncp, docp, dl, NULL, REAL(CONPARAMETER));
  }
  rdebug0(1,"Finished hdpMultinomial_iterate.\n");

  SEXP hdpout = PROTECT(duplicate(hdpin));
  rWriteHDP(hdpout,hdp);

  SEXP result = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(result, 0, hdpout);
  SET_VECTOR_ELT(result, 1, LIK);
  SET_VECTOR_ELT(result, 2, CONPARAMETER);
  SET_VECTOR_ELT(result, 3, ALLOCATIONS);


  if (nallocs > 0) {
    UNPROTECT(5);
  } else {
    UNPROTECT(4);
  }

  PutRNGstate();

  return result;
}


