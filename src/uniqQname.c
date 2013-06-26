#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "header.h"
#include "hash.h"

SEXP uniqQname(SEXP qname, SEXP totReadsR, SEXP pos, SEXP mpos, SEXP names){

  int totReads, hashSize, i, l, count=0, *qname_p;
  hash_t *hashP, myhash;
  PROTECT(totReadsR);
  totReads=INTEGER(totReadsR)[0];
  hashP = &myhash;
  hashSize=(int)(totReads/2);
  hash_init(hashP, hashSize);
  PROTECT(qname);
  PROTECT(pos);
  int *p_pos=INTEGER(pos);
  qname_p=INTEGER(qname);
  PROTECT(mpos);
  int *p_mpos=INTEGER(mpos);   

//vector to return new read ids
  PROTECT(names = coerceVector(names, STRSXP));

  char **tmpres;
  tmpres = malloc(floor(totReads/2) * sizeof(char*));
  for (i=0; i<floor(totReads/2); i++) tmpres[i] = malloc(200 * sizeof(char));

  char *tmp, *idtmp;
  tmp = malloc(200 * sizeof(int));
  idtmp = malloc(30 * sizeof(int));
  for (i=0; i<totReads; i++) {
    //strcpy(tmp,CHAR(STRING_ELT(qname, i)));
    sprintf(tmp, "%d", qname_p[i]);
    strcat(tmp, ".");
    if(p_pos[i]<p_mpos[i]) sprintf(idtmp, "%d",  p_pos[i]);
    else sprintf(idtmp, "%d",  p_mpos[i]);
    strcat(tmp, idtmp);
    strcat(tmp, ".");
    if(p_pos[i]<p_mpos[i]) sprintf(idtmp, "%d",  p_mpos[i]);
    else sprintf(idtmp, "%d",  p_pos[i]);
    strcat(tmp, idtmp);
    SET_STRING_ELT(names, i, mkChar(tmp));
    l=hash_lookup(hashP, tmp);
    if(l!=HASH_FAIL) {
      hash_update(hashP, tmp, l+1);
      if((l+1)==3) {
	strcpy(tmpres[count], CHAR(STRING_ELT(names, i)));
	count++;
      }
    } else hash_insert(hashP, tmp, 1);
  }
  free(tmp);
  free(idtmp);
  hash_destroy(hashP);

  SEXP res;
  PROTECT(res = allocVector(STRSXP, count));
  for (i=0; i<count; i++) {
    SET_STRING_ELT(res, i, mkChar(tmpres[i]));
  }

  SEXP ans;
  PROTECT(ans = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ans, 0, names);
  SET_VECTOR_ELT(ans, 1, res);
  
  UNPROTECT(7);
  return(ans);
}
