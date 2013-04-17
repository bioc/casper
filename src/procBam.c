#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "header.h"
#include "hash.h"
#include "functions.h"


SEXP procBam(SEXP qname, SEXP chr, SEXP start, SEXP cigar, SEXP strand, SEXP totFrags, SEXP totReads, SEXP totJunx, SEXP len, SEXP strs, SEXP key, SEXP chrom, SEXP rid, SEXP rstrand, SEXP jchrom, SEXP jstrs, SEXP jlen){
	read_t *frags;
	int totF, j, l, hashSize, i, *frags_size, *reads_size;
	hash_t *fragsHashPtr, fragsHash;
	hash_node_t *bucket;
	verbose=0;
    
	
	// Alloc memory for fragments
	PROTECT(totFrags = coerceVector(totFrags, INTSXP));
	frags_size = INTEGER(totFrags);
	frags = malloc((frags_size[0] + 1) * sizeof(read_t));
	
    	fragsHashPtr = &fragsHash;		
	hashSize=frags_size[0]+100;
	hash_init(fragsHashPtr, hashSize);
    
	PROTECT(qname = coerceVector(qname, STRSXP));
	PROTECT(chr = coerceVector(chr, STRSXP));
	PROTECT(start = coerceVector(start, INTSXP));
	int *p_start=INTEGER(start);
	PROTECT(cigar = coerceVector(cigar, STRSXP));
	PROTECT(totReads = coerceVector(totReads, INTSXP));
	PROTECT(strand = coerceVector(strand, INTSXP));
	int *p_strand=INTEGER(strand);
	reads_size = INTEGER(totReads);
	
	totF=0;
		// Find paired reads and fill fragment's object 
	for (i=0; i<frags_size[0]; i++) {
   	  l=hash_lookup(fragsHashPtr, CHAR(STRING_ELT(qname, i)));
	  if(l!=HASH_FAIL) {
            addRead2Frag(CHAR(STRING_ELT(qname, i)), CHAR(STRING_ELT(chr, i)), p_start[i], p_strand[i], CHAR(STRING_ELT(cigar, i)), l, frags, 2);
	  }
	  else {
            hash_insert(fragsHashPtr, CHAR(STRING_ELT(qname, i)), totF);
            addRead2Frag(CHAR(STRING_ELT(qname, i)), CHAR(STRING_ELT(chr, i)), p_start[i], p_strand[i], CHAR(STRING_ELT(cigar, i)), totF, frags, 1);
            totF++;
	  }
	}
	
	//SEXP len, strs, flag, key, chrom, rid;	
	int tmp, *cigs, counter=0, jcounter=0;
	
	PROTECT(len = coerceVector(len, INTSXP));
	int *p_len=INTEGER(len);
	PROTECT(jlen = coerceVector(jlen, INTSXP));
        int *p_jlen=INTEGER(jlen);

	PROTECT(chrom = coerceVector(chrom, STRSXP));
	PROTECT(jchrom = coerceVector(jchrom, STRSXP));

	PROTECT(strs = coerceVector(strs, INTSXP));
	int *p_strs=INTEGER(strs);
        PROTECT(jstrs = coerceVector(jstrs, INTSXP));
        int *p_jstrs=INTEGER(jstrs);

	PROTECT(key = coerceVector(key, STRSXP));
	PROTECT(rid = coerceVector(rid, INTSXP));
	int *p_rid=INTEGER(rid);
        PROTECT(rstrand = coerceVector(rstrand, INTSXP));
	int *p_rstrand = INTEGER(rstrand), ini;
	for(i=0; i<fragsHash.size; i++) {
	  if(fragsHash.bucket[i]!=NULL)  {
            bucket=fragsHash.bucket[i];
            while(bucket) {
	      tmp=bucket->data;
	      cigs=malloc(50 * sizeof(int));
	      if(frags[tmp].nreads==2) {
		cigs=procCigar(m_strdup(frags[tmp].cigar_1), cigs);
		frags[tmp].len_1=frags[tmp].st_1;
		ini=1;
		if(cigs[1]<0) ini=2;
		  //frags[tmp].len_1-=cigs[1];
		for(j=ini; j<cigs[0]+1;j++) {
		  if(cigs[j]>0){
		    SET_STRING_ELT(key, counter, mkChar(bucket->key));
		    SET_STRING_ELT(chrom, counter, mkChar(frags[tmp].chr_1));
		    p_strs[counter] = frags[tmp].len_1;
		    p_len[counter] = frags[tmp].len_1+(cigs[j]-1);
		    p_rid[counter] = 1;
		    p_rstrand[counter] = frags[tmp].strand_1;
		    frags[tmp].len_1 += cigs[j];
		    counter++;
		  } else {
		    if((cigs[0]>1)&&(j<cigs[0])) {
		      SET_STRING_ELT(jchrom, jcounter, mkChar(frags[tmp].chr_1));
		      p_jstrs[jcounter] = frags[tmp].len_1;
		      frags[tmp].len_1 -= cigs[j];		    
		      p_jlen[jcounter] = frags[tmp].len_1-1;
		      jcounter++;
		    }
		  }
		}
		cigs=procCigar(m_strdup(frags[tmp].cigar_2), cigs);
		frags[tmp].len_2=frags[tmp].st_2;
		ini=1;
                if(cigs[1]<0) ini=2;
		  //                  frags[tmp].len_2-=cigs[1];
		for(j=ini; j<cigs[0]+1;j++) {        
		  if(cigs[j]>0){
		    SET_STRING_ELT(key, counter, mkChar(bucket->key));
		    SET_STRING_ELT(chrom, counter, mkChar(frags[tmp].chr_2));
		    p_strs[counter] = frags[tmp].len_2;
		    p_len[counter] = frags[tmp].len_2+(cigs[j]-1);
		    p_rid[counter] = 2;
		    p_rstrand[counter] = frags[tmp].strand_2;
		    frags[tmp].len_2+=cigs[j];
		    counter++;
		  } else {
		    if((cigs[0]>1)&&(j<cigs[0])) {
		      SET_STRING_ELT(jchrom, jcounter, mkChar(frags[tmp].chr_2));
		      p_jstrs[jcounter] = frags[tmp].len_2;
		      frags[tmp].len_2 -= cigs[j];
		      p_jlen[jcounter] = frags[tmp].len_2-1;
		      jcounter++;
		    }
		  }
		}
		free(frags[tmp].chr_2);
		free(frags[tmp].cigar_2);
  	      }
	      free(frags[tmp].qname);
	      free(frags[tmp].chr_1);
	      free(frags[tmp].cigar_1);
	      free(cigs);
	      bucket=bucket->next;
            }
	  }
	}
	SEXP reads;
	PROTECT(reads = allocVector(VECSXP, 9));
	SET_VECTOR_ELT(reads, 0, len);
	SET_VECTOR_ELT(reads, 1, strs);
	SET_VECTOR_ELT(reads, 2, key);
	SET_VECTOR_ELT(reads, 3, chrom);
	SET_VECTOR_ELT(reads, 4, rid);
	SET_VECTOR_ELT(reads, 5, rstrand);
	SET_VECTOR_ELT(reads, 6, jchrom);
        SET_VECTOR_ELT(reads, 7, jstrs);
        SET_VECTOR_ELT(reads, 8, jlen);
	free(frags);
	hash_destroy(fragsHashPtr);
	UNPROTECT(17);
	return(reads);
}
