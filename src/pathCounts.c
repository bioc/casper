#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "header.h"
#include "hash.h"
#include "functions.h"
#include "fragFunc.h"

void countPaths(int pos, path_t *frags, hash_t *pathsHashPtr);
int sort(const void *x, const void *y); //int sort(const void *x, const void *y);
int sort2(const void *x, const void *y);
void addPath(int *unex, int *unrid, hash_t *hash, int totEx);

SEXP pathCounts(SEXP reid, SEXP rid, SEXP exst, SEXP exid){
    //void pathCounts(char **reid, int *rid, int *exst, int *exid, int nreads){
    
    int totF, hashSize, i, *p_rid, *p_exst, *p_exid, nreads;
    hash_node_t *bucket;
    path_t *frags;
    hash_t  *pathsHashPtr, pathsHash, *fragsHashPtr, fragsHash;
    pathsHashPtr=&pathsHash;
    fragsHashPtr=&fragsHash;
    hashSize=pow(2,25);	
    verbose=0;
	
    PROTECT(rid = coerceVector(rid, INTSXP));
    PROTECT(exst = coerceVector(exst, INTSXP));
    PROTECT(exid = coerceVector(exid, INTSXP));
    PROTECT(reid = coerceVector(reid, STRSXP));
    nreads=length(rid);
    
    hash_init(fragsHashPtr, hashSize);
    hash_init(pathsHashPtr, hashSize);
    
    p_rid=INTEGER(rid);
    p_exst=INTEGER(exst);
    p_exid=INTEGER(exid);
    
    path_t **pfrags=&frags;
    totF = buildFrags(fragsHashPtr, reid, p_rid, p_exst, p_exid, nreads, pfrags);
        
    for(i=0; i<fragsHash.size; i++) {
        if(fragsHash.bucket[i]!=NULL)  {
            bucket=fragsHash.bucket[i];
            while(bucket) {
                countPaths(bucket->data, frags, pathsHashPtr);
	        free(frags[bucket->data].exons);
                free(frags[bucket->data].starts);
                free(frags[bucket->data].rids);
                bucket=bucket->next;
            }
        }
    }
    
    SEXP key, pathc, tot;
    PROTECT(key = allocVector(STRSXP, nreads));
    PROTECT(pathc = allocVector(INTSXP, nreads));
    PROTECT(tot = allocVector(INTSXP, 1));
    
    int *p_pathc, count, *p_tot;
    count=0;
    
    p_pathc=INTEGER(pathc);
    p_tot=INTEGER(tot);
    
    for(i=0; i<pathsHash.size; i++) {
        if(pathsHash.bucket[i]!=NULL)  {
            bucket=pathsHash.bucket[i];
            while(bucket) {
                SET_STRING_ELT(key, count, mkChar(bucket->key));
                p_pathc[count] = bucket->data;
                bucket=bucket->next;
                count++;
            }
        }
    }
    
    p_tot[0] = count;
    
    SEXP ans;
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, key);
    SET_VECTOR_ELT(ans, 1, pathc);
    SET_VECTOR_ELT(ans, 2, tot);
    
    hash_destroy(pathsHashPtr);
    hash_destroy(fragsHashPtr);
    UNPROTECT(8);
    free(frags);
    return(ans);
}


void countPaths(int pos, path_t *frags, hash_t *pathsHashPtr){
    int i, *ord, *unex, *unrid, **standex, totEx;
    
    ord=malloc((frags[pos].nexon+1) * sizeof(int));
    unex=malloc((frags[pos].nexon+1) * sizeof(int));
    unrid=malloc((frags[pos].nexon+1) * sizeof(int));
    standex=malloc((frags[pos].nexon+1) * sizeof(int *));
    for(i=0; i<frags[pos].nexon; i++) standex[i]=malloc(4 * sizeof(int));
    
    for (i=0; i<frags[pos].nexon; i++) {
        standex[i][0]=frags[pos].starts[i];
        standex[i][1]=frags[pos].exons[i];
        standex[i][2]=frags[pos].rids[i];
    }
    
    qsort(standex, frags[pos].nexon, sizeof(int **), sort);
    
    totEx=0;
    unex[totEx]=standex[0][1];
    unrid[totEx]=standex[0][2];
    totEx++;
    for(i=1; i<frags[pos].nexon; i++){
      if((standex[i][0] != standex[i-1][0])||(standex[i][2] != standex[i-1][2])) {
	unex[totEx]=standex[i][1];
	unrid[totEx]=standex[i][2];
	totEx++;
      }
    }
    
    addPath(unex, unrid, pathsHashPtr, totEx);
    free(ord);
    free(unex);
    free(unrid);
    for(i=0; i<frags[pos].nexon; i++) free(standex[i]);
    free(standex);
    
}

void addPath(int *unex, int *unrid, hash_t *hash, int totEx){
    int l;
    char *pastr, *tmp;
    tmp=malloc(50 * sizeof(char));
    pastr=malloc((totEx+1)*50 * sizeof(char));
    verbose=0;
    //Check for overlapping ends
    int nleft=0, nright=0, *lread, *rread;
    lread=malloc((totEx+1) * sizeof(int));
    rread=malloc((totEx+1) * sizeof(int));
    for(l=0; l<totEx; l++){
      if(unrid[l]==1) {
	lread[nleft]=unex[l];
	nleft++;
      } else {
	rread[nright]=unex[l];
	nright++;
      }
    }

  
    strcpy(pastr, ".");
    sprintf(tmp, "%d", unex[0]);
    strcat(pastr, tmp);

    if(totEx>1){
      for(l=1; l<nleft; l++){
	strcat(pastr, ".");
	sprintf(tmp, "%d", lread[l]);
	strcat(pastr, tmp);
      }
      strcat(pastr, "-");
      for(l=0; l<nright; l++){
	sprintf(tmp, "%d", rread[l]);
	strcat(pastr, tmp);
	strcat(pastr, ".");
      } 
    }

   
    l=hash_lookup(hash, pastr);
    if(l!=HASH_FAIL) {
        hash_update(hash, pastr, l+1);
    }
    else {
        hash_insert(hash, pastr, 1); 
    }
    free(lread);
    free(rread);
    free(tmp);
    free(pastr);  
}

int sort(const void *x, const void *y) {
    int *a = *(int**)x, *b= *(int**)y;
    int res;
    if((a[2]-b[2])!=0) res = a[2]-b[2];
    else res = a[0]-b[0];
    return( res );
}

int sort2(const void *x, const void *y) {
    int *a = *(int**)x, *b= *(int**)y;
    int res;
    res = (a[0] - b[0]) * (a[1] - b[1]);
    return( a[2] - b[2] );
}

