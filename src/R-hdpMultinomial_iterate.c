#include <stdlib.h>
#include "R-hdpMultinomial_iterate.h"

SEXP hdpMultinomial_iterate(SEXP hdpin, SEXP numiter, SEXP doconparam, SEXP dolik, SEXP numallocs, SEXP numcp, SEXP dodebug)
{
  int ni, docp, dl, nallocs, ncp, nclust;

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
  SEXP NUMCLUST = PROTECT(allocVector(INTSXP, ni));

  SEXP ALLOCATIONS;

  if (nallocs > 0) {
    ALLOCATIONS = PROTECT(allocVector(INTSXP, nallocs));
  } else {
    ALLOCATIONS = R_NilValue;
  }
  rdebug0(1,"Running hdpMultinomial_iterate.\n");
  if (nallocs > 0) {
    hdp_iterate(hdp, REAL(LIK), ni, ncp, docp, dl, INTEGER(ALLOCATIONS), REAL(CONPARAMETER), INTEGER(NUMCLUST));
  } else {
    hdp_iterate(hdp, REAL(LIK), ni, ncp, docp, dl, NULL, REAL(CONPARAMETER), INTEGER(NUMCLUST));
  }
  rdebug0(1,"Finished hdpMultinomial_iterate.\n");

  SEXP hdpout = PROTECT(duplicate(hdpin));
  rWriteHDP(hdpout,hdp);

  SEXP result = PROTECT(allocVector(VECSXP, 5));
  SET_VECTOR_ELT(result, 0, hdpout);
  SET_VECTOR_ELT(result, 1, LIK);
  SET_VECTOR_ELT(result, 2, CONPARAMETER);
  SET_VECTOR_ELT(result, 3, ALLOCATIONS);
  SET_VECTOR_ELT(result, 4, NUMCLUST);


  if (nallocs > 0) {
    UNPROTECT(6);
  } else {
    UNPROTECT(5);
  }

  PutRNGstate();

  return result;
}


