#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "header.h"
#include "hash.h"
#include "simReadsfunc.h"
#include <R_ext/Utils.h>
     


SEXP casperSimC(SEXP gene_exp, SEXP var_exp, SEXP var_num, SEXP var_len, SEXP exon_num, SEXP exon_st, SEXP exon_end, SEXP exon_id, SEXP len_distrV, SEXP len_distrD, SEXP st_distrV, SEXP st_distrD, SEXP read_len, SEXP nn, SEXP tx_strand, SEXP lr_fileR, SEXP chr, SEXP rseed, SEXP rbam, SEXP rinsideBam, SEXP verbose){

  int i=0, *ge, *vn, *vl, *en, *es, *ee, *ei, *txstr, ngenes, rl, sdlen, n, bam, insideBam;
  double *ve, *sdv, *sdd, *ldd, *ldv;
  FILE *LRFILE=NULL;
  SEXP startsTmp;

  PROTECT(gene_exp);// = coerceVector(gene_exp, INTSXP));
  PROTECT(var_exp);// = coerceVector(var_exp, REALSXP));
  PROTECT(var_num);// = coerceVector(var_num, INTSXP));
  PROTECT(var_len);// = coerceVector(var_len, INTSXP));
  PROTECT(exon_num);// = coerceVector(exon_num, INTSXP));
  PROTECT(exon_st);// = coerceVector(exon_st, INTSXP));
  PROTECT(exon_end);// = coerceVector(exon_end, INTSXP));
  PROTECT(exon_id);// = coerceVector(exon_id, INTSXP));
  PROTECT(tx_strand);// = coerceVector(tx_strand, INTSXP));
  PROTECT(len_distrV);// = coerceVector(len_distrV, INTSXP));
  PROTECT(len_distrD);// = coerceVector(len_distrD, REALSXP));
  PROTECT(st_distrV);// = coerceVector(st_distrV, REALSXP));
  PROTECT(st_distrD);// = coerceVector(st_distrD, REALSXP));
  PROTECT(read_len);// = coerceVector(read_len, INTSXP));
  PROTECT(nn);// = coerceVector(nn, INTSXP));
  PROTECT(lr_fileR);// = AS_CHARACTER(lr_fileR));
  PROTECT(chr);// = coerceVector(chr, STRSXP));
  PROTECT(rbam);// = coerceVector(rbam, INTSXP));
  PROTECT(rinsideBam);// = coerceVector(rinsideBam, INTSXP));
  PROTECT(verbose);



  ge = INTEGER(gene_exp);
  ve = REAL(var_exp);
  vn = INTEGER(var_num);
  vl = INTEGER(var_len);
  en = INTEGER(exon_num);
  es = INTEGER(exon_st);
  ee = INTEGER(exon_end);
  ei = INTEGER(exon_id);
  txstr = INTEGER(tx_strand);
  ldv = REAL(len_distrV);
  ldd = REAL(len_distrD);
  sdv = REAL(st_distrV);
  sdd = REAL(st_distrD);
  rl = INTEGER(read_len)[0];
  n = INTEGER(nn)[0];
  bam = INTEGER(rbam)[0];
  insideBam = INTEGER(rinsideBam)[0];
  ngenes = length(var_num);
  sdlen = length(st_distrD);
  int ldlen = length(len_distrV);
  gene_t *genes;
  genes = malloc((ngenes+1) * sizeof(gene_t));
  build_genes(genes, ve, vn, vl, en, es, ee, ei, txstr, ngenes, chr);
  
  int gene, var, len, *gansS, *vansS, *lansS, *strS, seed=INTEGER(rseed)[0], *pos, st;
  double *sansS;


  SEXP gans, vans, lans, sans, ans, strs, qname, rname, strand, posr, cigar;
  PROTECT(gans = allocVector(INTSXP, n));
  PROTECT(vans = allocVector(INTSXP, n));
  PROTECT(lans = allocVector(INTSXP, n));
  PROTECT(sans = allocVector(REALSXP, n));
  PROTECT(strs = allocVector(INTSXP, n));
  PROTECT(ans = allocVector(VECSXP, 12));
  if(insideBam==1){
    PROTECT(qname = allocVector(STRSXP, n*2)); 
    PROTECT(rname = allocVector(STRSXP, n*2));
    PROTECT(strand = allocVector(STRSXP, n*2));
    PROTECT(posr = allocVector(INTSXP, n*2));
    PROTECT(cigar = allocVector(STRSXP, n*2));
  } else {
    PROTECT(qname = allocVector(STRSXP, 1));
    PROTECT(rname = allocVector(STRSXP, 1));
    PROTECT(strand = allocVector(STRSXP, 1));
    PROTECT(posr = allocVector(INTSXP, 1));
    PROTECT(cigar = allocVector(STRSXP, 1));
  }
  gansS = INTEGER(gans);
  vansS = INTEGER(vans);
  lansS = INTEGER(lans);
  sansS = REAL(sans);
  strS = INTEGER(strs);
  pos = INTEGER(posr);
  
   if(bam==1) LRFILE = fopen(CHAR(STRING_ELT(lr_fileR, 0)), "a");


  int j, *starts, gap, totp=0;
  char **cigars = NULL;
  char seqstr[rl+1], seqnuc[2]="=", tmpchar[100];
  hash_t *paths, paths_pted;
  PROTECT(startsTmp=allocVector(INTSXP, 3));
  starts = INTEGER(startsTmp);
  
  if(bam==1){
    cigars = malloc(3 * sizeof(char *));
    for(i=0; i<3; i++) cigars[i] = malloc(100 * sizeof(char));
    strcpy(seqstr, seqnuc);
    for (i=0; i<rl-1; i++) strcat(seqstr, seqnuc);
  }
  
  paths = &paths_pted;
  hash_init(paths, NextPow2(length(exon_st)*20)); 

  int cnt=0;
  char geStr[2], *last;
  
  i=0;

  while(i<n) {
    j=0;
    cnt=0;
    gene = ge[i];
    var = choose_var(genes[gene]);
    //printf("%d %d\n", var, i);
    st=-1;
    //len = ldv[i];
    len = choose_len(genes[gene].vars[var].len, ldv, ldd, ldlen);
    if(len>genes[gene].vars[var].len) { len=genes[gene].vars[var].len; }
    while(j==0){
      cnt++;
      if(len==genes[gene].vars[var].len) st=1; 
      else st = choose_st(len, genes[gene].vars[var].len, sdv, sdd, sdlen, genes[gene].vars[var].strand);
      if(st>=0) {        
	//printf("%d\n", st);
	if(bam==1){
          starts=build_cigar(genes[gene].vars[var], len, st, rl, cigars, genes[gene].vars[var].strand);
          if(genes[gene].vars[var].strand==1) vansS[i] = starts[0];
	  else vansS[i] = starts[1];
	  if(genes[gene].vars[var].strand == 1) gap = starts[1]-(starts[0]+rl);
	  else gap = starts[0]-(starts[1]+rl);
	  if(genes[gene].vars[var].strand==1) strcpy(geStr, "+");
	  else strcpy(geStr, "-");
	  last = cigars[0] + (int) (strlen((const char *)cigars[0]))-1;
	  if(strcmp(last, "N")!=0){
	    last = cigars[1] + (int) (strlen((const char *)cigars[1]))-1;
	    if(strcmp(last, "N")!=0){
	      if((strcmp(cigars[0], "\0")!=0)&&(strcmp(cigars[1], "\0")!=0)) fprintf(LRFILE, "%s.%d.%d\t99\t%s\t%d\t64\t%s\t=\t%d\t%d\t%s\t%s\tXS:A:%s\tRG:Z:%d\n", genes[gene].chr, i, var+1, genes[gene].chr, starts[0], cigars[0], starts[1], gap, seqstr, seqstr, geStr, var+1);

	      if(strcmp("", CHAR(STRING_ELT(lr_fileR, 0)))!=0) {
		if((strcmp(cigars[0], "\0")!=0)&&(strcmp(cigars[1], "\0")!=0)) fprintf(LRFILE, "%s.%d.%d\t147\t%s\t%d\t64\t%s\t=\t%d\t%d\t%s\t%s\tXS:A:%s\tRG:Z:%d\n", genes[gene].chr, i, var+1, genes[gene].chr, starts[1], cigars[1], starts[0], gap, seqstr, seqstr, geStr, var+1);
	      }
	    }
	  }
	}
	
	build_path(genes[gene].vars[var], len, st, rl, paths, genes[gene].vars[var].strand, starts);
	if(genes[gene].vars[var].strand==1) vansS[i] = starts[0];
        else vansS[i] = starts[1];
	
	if(insideBam==1){
	  sprintf(tmpchar, "%d.%d", i, var+1); SET_STRING_ELT(qname, i*2, mkChar(tmpchar)); SET_STRING_ELT(qname, i*2+1, mkChar(tmpchar));
	  SET_STRING_ELT(rname, i*2, mkChar(genes[gene].chr)); SET_STRING_ELT(rname, i*2+1, mkChar(genes[gene].chr));
	  if(genes[gene].vars[var].strand==0) { SET_STRING_ELT(strand, i*2, mkChar("-")); SET_STRING_ELT(strand, i*2+1, mkChar("-")); }
	  else { SET_STRING_ELT(strand, i*2, mkChar("+")); SET_STRING_ELT(strand, i*2+1, mkChar("+")); }
	  pos[i*2] = starts[3]; pos[i*2+1] = starts[4];
	  SET_STRING_ELT(cigar, i*2, mkChar(cigars[0])); SET_STRING_ELT(cigar, i*2+1, mkChar(cigars[1]));
	}
	totp+=starts[2];
	//free(starts);
	strS[i] = genes[gene].vars[var].strand;
	gansS[i] = genes[gene].vars[var].len;
	lansS[i] = len;
	sansS[i] = st;
	i++;
	j=1;
	if(i % (int) floor(n*0.1) == 0) {
          if(INTEGER(verbose)[0]==1) Rprintf("%d %% of fragments simulated\n", (int) ceil(((double)i/(double)n)*100));
	  R_CheckUserInterrupt();
	}
      }
      if(cnt==1000) {
	Rprintf("Fragment lost %d %d\n", i, n);
	break;
      }
    }
  }

  hash_node_t *bucket;
  int count=0, *p_pathc;
  SEXP key, pathc;
  PROTECT(key = allocVector(STRSXP, totp));  
  PROTECT(pathc = allocVector(INTSXP, totp));
  p_pathc = INTEGER(pathc);

  
  for(i=0; i<paths_pted.size; i++) {
    if(paths_pted.bucket[i]!=NULL) {
      bucket=paths_pted.bucket[i];
      while(bucket) {
	SET_STRING_ELT(key, count, mkChar(bucket->key));
	p_pathc[count] = bucket->data;
	bucket=bucket->next;
	count++;
      }
    }
  }
  hash_destroy(paths);
         
  
  SET_VECTOR_ELT(ans, 0, gans);
  SET_VECTOR_ELT(ans, 1, vans);
  SET_VECTOR_ELT(ans, 2, lans);
  SET_VECTOR_ELT(ans, 3, sans);
  SET_VECTOR_ELT(ans, 4, strs);
  SET_VECTOR_ELT(ans, 5, key);
  SET_VECTOR_ELT(ans, 6, pathc);
  SET_VECTOR_ELT(ans, 7, qname);
  SET_VECTOR_ELT(ans, 8, rname);
  SET_VECTOR_ELT(ans, 9, strand);
  SET_VECTOR_ELT(ans, 10, posr);
  SET_VECTOR_ELT(ans, 11, cigar);

  UNPROTECT(34);

  for (i=0; i<ngenes; i++){
    free(genes[i].chr);
    for(j=0; j<vn[i]; j++) {
      free(genes[i].vars[j].exst);
      free(genes[i].vars[j].exen);
      free(genes[i].vars[j].exid);
    }
    free(genes[i].vars);
  }
  free(genes);   
  if(bam==1){
    for(i=0; i<3; i++) free(cigars[i]);
    free(cigars);
    fclose(LRFILE);
  }
  return(ans);

}
