#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "hash.h"

int compare (const void * a, const void * b);
int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}


SEXP joinExons(SEXP sexons, SEXP sreads, SEXP stot){

  int hashSize, i, *exons, *reads, len, l, **links, tot;
  hash_t *myhashP, myhash;

  PROTECT(sexons = coerceVector(sexons, INTSXP));
  PROTECT(sreads = coerceVector(sreads, INTSXP));
  PROTECT(stot = coerceVector(stot, INTSXP));
  len=length(sexons);

  exons = INTEGER(sexons);
  reads = INTEGER(sreads);
  tot = INTEGER(stot)[0];

  myhashP = &myhash;
  hashSize = tot;
  hash_init(myhashP, hashSize);


  links=malloc((tot+1) * sizeof(int *));
  for(i=0; i<tot; i++) links[i] = malloc(50 * sizeof(int));

  int counter=0;

  char id[100];
  for(i=0; i<len; i++){
    sprintf(id, "%d", reads[i]);
    l=hash_lookup(myhashP, id);
    if(l!=HASH_FAIL) {
      links[l][0]++;
      if(links[l][0] % 49 == 0) links[l] = realloc(links[l], (links[l][0]+50) * sizeof(int)); 
      links[l][links[l][0]] = exons[i];
   }
    else {
      hash_insert(myhashP, id, counter);
      links[counter][0] = 1;
      links[counter][1] = exons[i];
      counter++;
      if(counter>=tot) break;
    }
  }

  char **ans;
  int j=0, k, chk=0, finalSize;
  ans = malloc((counter+1) * sizeof(char *));
  
  hash_destroy(myhashP);  
  hash_init(myhashP, counter);
  for(i=0; i<counter; i++) {
    if(links[i][0]>1){
      qsort(&links[i][1], links[i][0], sizeof(int), compare);
      chk=0;
      for(k=2; k<links[i][0]+1; k++) if(links[i][k-1] != links[i][k]) chk++;
      if(chk>0){
	ans[j] = malloc(15 * (links[i][0]+1) * sizeof(char));
	sprintf(id, "%d", links[i][1]);
	strcpy(ans[j], id);
	strcat(ans[j], ".");
	for(k=2; k<links[i][0]+1; k++) {
	  if(links[i][k] !=links[i][k-1]){
	    sprintf(id, "%d", links[i][k]);
	    strcat(ans[j], id);
	    strcat(ans[j], ".");
	  }
	}
	l=hash_lookup(myhashP, ans[j]);
	if(l!=HASH_FAIL) hash_update(myhashP, ans[j], l+1);
	else hash_insert(myhashP, ans[j], 1);
	j++;
      }
    }
  }
  finalSize=j;
 
  int *tmpcounts, *pcounts;
  tmpcounts = malloc(finalSize * sizeof(int));
  char **tmpkey;
  tmpkey = malloc(finalSize * sizeof(char *));
  for(i=0; i<finalSize; i++) tmpkey[i] = malloc(200 * sizeof(char));
  hash_node_t *bucket;
  j=0;
  for(i=0; i<myhash.size; i++) {
    if(myhash.bucket[i]!=NULL) {
      bucket=myhash.bucket[i];
      while(bucket) {
	//Malloc here space for tmpkey instead of fixed 200 chars
	tmpkey[j] = malloc((strlen(bucket->key)+1) * sizeof(char));
	strcpy(tmpkey[j], bucket->key);
	tmpcounts[j] = bucket->data;
	bucket = bucket->next;
	j++;
      }
    }
  }
  int ksize=j;
 
  SEXP res;
  SEXP keys;
  SEXP counts;
  PROTECT(keys = allocVector(STRSXP, j));
  PROTECT(counts = allocVector(INTSXP, j));
  PROTECT(res = allocVector(VECSXP, 2));
  pcounts = INTEGER(counts);

  for(i=0; i<j; i++) {
    SET_STRING_ELT(keys, i, Rf_mkChar(tmpkey[i]));
    pcounts[i] = tmpcounts[i];
  }
  SET_VECTOR_ELT(res, 0, keys);
  SET_VECTOR_ELT(res, 1, counts);

  for(i=0; i<finalSize; i++) free(ans[i]);
  for(i=0; i<ksize; i++) free(tmpkey[i]);
  for(i=0; i<tot; i++) free(links[i]);
  
  free(ans);
  free(links);
  free(tmpcounts);
  UNPROTECT(6);
  return(res);
}

