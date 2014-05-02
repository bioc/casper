#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "header.h"
#include "hash.h"
#include "functions.h"


SEXP procBam(SEXP qname, SEXP chr, SEXP start, SEXP mpos, SEXP cigar, SEXP strand, SEXP totFrags, SEXP totReads, SEXP flag, SEXP totJunx, SEXP len, SEXP strs, SEXP key, SEXP chrom, SEXP rid, SEXP rstrand, SEXP jchrom, SEXP jstrs, SEXP jlen, SEXP rflag){
	read_t *frags;
	int totF, j, l, hashSize, i, frags_size, *qname_p;
	hash_t *fragsHashPtr, fragsHash;
	hash_node_t *bucket;
	char *echrom="\0";
	verbose=0;
	
	// Alloc memory for fragments
	PROTECT(totFrags);
	frags_size = INTEGER(totFrags)[0];
	frags = malloc((frags_size + 1) * sizeof(read_t));
	
    	fragsHashPtr = &fragsHash;		
	hashSize=frags_size+100;hash_init(fragsHashPtr, hashSize);
    
	PROTECT(qname);
	qname_p = INTEGER(qname);
	PROTECT(chr);
	PROTECT(start);
	int *p_start=INTEGER(start);
	PROTECT(mpos);
	int *p_mpos=INTEGER(mpos);
	PROTECT(cigar);
	PROTECT(totReads);
	PROTECT(strand);
	int *p_strand=INTEGER(strand);
	int *p_flag=INTEGER(flag);
	//int reads_size = INTEGER(totReads)[0];
	char *str;
	str=malloc(100 * sizeof(char));
	totF=0;

      	// Find paired reads and fill fragment's object 
	for (i=0; i<frags_size; i++) {
	  if(p_start[i]<p_mpos[i]) sprintf(str, "%d.%d.%d", qname_p[i], p_start[i], p_mpos[i]);
	  else sprintf(str, "%d.%d.%d", qname_p[i], p_mpos[i], p_start[i]);
	  l=hash_lookup(fragsHashPtr, (char *) str);//CHAR(STRING_ELT(qname, i)));
	  if(l!=HASH_FAIL) {
	    if(length(chr)>1) addRead2Frag((char *) str, CHAR(STRING_ELT(chr, i)), p_start[i], p_strand[i], i, l, frags, 2);
	    else addRead2Frag((char *) str, echrom, p_start[i], p_strand[i], i, l, frags, 2);
	  }
	  else {
            hash_insert(fragsHashPtr, (char *) str, totF);
	    if(length(chr)>1) addRead2Frag((char *) str, CHAR(STRING_ELT(chr, i)), p_start[i], p_strand[i], i, totF, frags, 1);
	    else addRead2Frag((char *) str, echrom, p_start[i], p_strand[i], i, totF, frags, 1);
            totF++;
	  }
	}
	free(str);
	int tmp, *cigs, counter=0, jcounter=0;
	
	PROTECT(len);
	int *p_len=INTEGER(len);
	PROTECT(jlen);
        int *p_jlen=INTEGER(jlen);

	PROTECT(chrom);
	PROTECT(jchrom);

	PROTECT(strs);
	int *p_strs=INTEGER(strs);
        PROTECT(jstrs);
        int *p_jstrs=INTEGER(jstrs);

	PROTECT(key);
	PROTECT(rid);
	int *p_rid=INTEGER(rid);
        PROTECT(rstrand);
	int *p_rstrand = INTEGER(rstrand), ini;
	int *p_rflag = INTEGER(rflag);

	for(i=0; i<fragsHash.size; i++) {
	  if(fragsHash.bucket[i]!=NULL)  {
            bucket=fragsHash.bucket[i];
            while(bucket) {
	      tmp=bucket->data;
	      cigs=malloc(100 * sizeof(int));
	      if(frags[tmp].nreads==2) {
		cigs = procCigar(m_strdup(CHAR(STRING_ELT(cigar, frags[tmp].strand_1))), cigs);
		frags[tmp].len_1 = p_start[frags[tmp].strand_1];
		ini=1;
		if(cigs[1]<0) ini=2;
		for(j=ini; j<cigs[0]+1;j++) {
		  if(cigs[j]>0){
		    SET_STRING_ELT(key, counter, mkChar(bucket->key));
		    if(length(chr)>1) SET_STRING_ELT(chrom, counter, mkChar(CHAR(STRING_ELT(chr, frags[tmp].strand_1))));
		    if(length(rflag)>1) p_rflag[counter] = p_flag[frags[tmp].strand_1];
		    p_strs[counter] = frags[tmp].len_1;
		    p_len[counter] = frags[tmp].len_1+(cigs[j]-1);
		    p_rid[counter] = 1;
		    p_rstrand[counter] = p_strand[frags[tmp].strand_1];
		    frags[tmp].len_1 += cigs[j];
		    counter++;
		  } else {
		    if(INTEGER(totJunx)[0]>1){
		      if((cigs[0]>1)&&(j<cigs[0])) {
			if(length(chr)>1) SET_STRING_ELT(jchrom, jcounter, mkChar(CHAR(STRING_ELT(chr, frags[tmp].strand_1))));
			p_jstrs[jcounter] = frags[tmp].len_1;
			frags[tmp].len_1 -= cigs[j];		    
			p_jlen[jcounter] = frags[tmp].len_1-1;
			jcounter++;
		      }
		    }
		  }
		}
		cigs = procCigar(m_strdup(CHAR(STRING_ELT(cigar, frags[tmp].strand_2))), cigs);
		ini=1;
                if(cigs[1]<0) {
		  frags[tmp].len_2=p_start[frags[tmp].strand_2] + cigs[1]*-1;
		  ini=2;
		} else frags[tmp].len_2=p_start[frags[tmp].strand_2];
		for(j=ini; j<cigs[0]+1;j++) {        
		  if(cigs[j]>0){
		    SET_STRING_ELT(key, counter, mkChar(bucket->key));
		    if(length(chr)>1) SET_STRING_ELT(chrom, counter, mkChar(CHAR(STRING_ELT(chr, frags[tmp].strand_2))));
		    if(length(rflag)>1) p_rflag[counter] = p_flag[frags[tmp].strand_2];
		    p_strs[counter] = frags[tmp].len_2;
		    p_len[counter] = frags[tmp].len_2+(cigs[j]-1);
		    p_rid[counter] = 2;
		    p_rstrand[counter] = p_strand[frags[tmp].strand_2];
		    frags[tmp].len_2+=cigs[j];
		    //printf("I'm here %d %d %d %d %d %d\n", frags[tmp].len_2, counter, p_strs[counter], p_len[counter], cigs[j], j);
		    counter++;
		  } else {
		    if(INTEGER(totJunx)[0]>1){
		      if((cigs[0]>1)&&(j<cigs[0])) {
			if(length(chr)>1) SET_STRING_ELT(jchrom, jcounter, mkChar(CHAR(STRING_ELT(chr, frags[tmp].strand_2))));
			p_jstrs[jcounter] = frags[tmp].len_2;
			frags[tmp].len_2 -= cigs[j];
			p_jlen[jcounter] = frags[tmp].len_2-1;
			jcounter++;
		      }
		    }
		  }
		}
  	      }
	      free(cigs);
	      bucket=bucket->next;
            }
	  }
	}

	SEXP reads;
	PROTECT(reads = allocVector(VECSXP, 10));
	SET_VECTOR_ELT(reads, 0, len);
	SET_VECTOR_ELT(reads, 1, strs);
	SET_VECTOR_ELT(reads, 2, key);
	SET_VECTOR_ELT(reads, 3, chrom);
	SET_VECTOR_ELT(reads, 4, rid);
	SET_VECTOR_ELT(reads, 5, rstrand);
	SET_VECTOR_ELT(reads, 6, jchrom);
        SET_VECTOR_ELT(reads, 7, jstrs);
        SET_VECTOR_ELT(reads, 8, jlen);
        SET_VECTOR_ELT(reads, 9, rflag);
	free(frags);
	hash_destroy(fragsHashPtr);
	UNPROTECT(18);
	return(reads);
}
