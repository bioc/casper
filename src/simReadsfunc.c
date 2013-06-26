#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include "hash.h"
#include "header.h"

void build_genes(gene_t *genes, double *ve, int *vn, int *vl, int *en, int *es, int *ee, int *ei, int *txstr, int ngenes, SEXP chr){
  int i, j, k, varpos=0, expos=0, tmp;
  
  for (i=0; i<ngenes; i++){
    genes[i].nvar = vn[i];
    genes[i].chr = malloc((strlen(CHAR(STRING_ELT(chr,i)))+1) * sizeof(char));
    strcpy(genes[i].chr, CHAR(STRING_ELT(chr, i)));
    genes[i].vars = malloc((vn[i]+1) * sizeof(var_t));
    for(j=0; j<vn[i]; j++){
      genes[i].vars[j].len = vl[varpos];
      genes[i].vars[j].nex = en[varpos];
      genes[i].vars[j].exp = ve[varpos];
      genes[i].vars[j].strand = txstr[varpos];
      genes[i].vars[j].exst = malloc((genes[i].vars[j].nex+1) * sizeof(int));
      genes[i].vars[j].exen = malloc((genes[i].vars[j].nex+1) * sizeof(int));
      genes[i].vars[j].exid = malloc((genes[i].vars[j].nex+1) * sizeof(int));
      for(k=0; k<genes[i].vars[j].nex; k++){
	genes[i].vars[j].exst[k] = es[expos];
	genes[i].vars[j].exen[k] = ee[expos];
	genes[i].vars[j].exid[k] = ei[expos];
	expos++;
      }
      if(genes[i].vars[j].strand==-1) {
	for(k=0; k < genes[i].vars[j].nex/2; k++){
	  tmp = genes[i].vars[j].exst[k];
	  genes[i].vars[j].exst[k] = genes[i].vars[j].exst[genes[i].vars[j].nex - k -1];
	  genes[i].vars[j].exst[genes[i].vars[j].nex - k -1] = tmp;
	  tmp = genes[i].vars[j].exen[k];
          genes[i].vars[j].exen[k] = genes[i].vars[j].exen[genes[i].vars[j].nex - k -1];
          genes[i].vars[j].exen[genes[i].vars[j].nex - k -1] = tmp;
	  tmp = genes[i].vars[j].exid[k]; 
	  genes[i].vars[j].exid[k] = genes[i].vars[j].exid[genes[i].vars[j].nex - k -1]; 
	  genes[i].vars[j].exid[genes[i].vars[j].nex - k -1] = tmp; 

	}
	}
      varpos++;
    }
  }
}
  

//##########################
//---------- Choose Gene
//###########################

int choose_gene(double *exp, int ngenes){
  int i;
  double ran, tmp=0;

  ran = rand() / ( RAND_MAX + 1.0 );
   for(i=0; i<ngenes; i++){
     if((tmp<=ran) && (ran<tmp+exp[i])) return(i);
    tmp+=exp[i];
  }
  Rprintf("Error: no gene chosen\n");
  return(0);
}

//#############------------
//-------- Choose variant
//############-----------

int choose_var(gene_t gene){
  int i;
  double ran, tmp=0;
  ran = (double)rand() / (double)( RAND_MAX - 1);
  for(i=0; i<gene.nvar; i++){
    if((tmp<=ran) && (ran<tmp+gene.vars[i].exp)) return(i);
    tmp+=gene.vars[i].exp;
  }
  tmp=0;
  Rprintf("Error: no variant chosen: %d\n", gene.nvar);
  for(i=0; i<gene.nvar; i++){
    tmp+=gene.vars[i].exp;    
    Rprintf("%f %f\n", gene.vars[i].exp, tmp);
  }
  //  exit(0);
  return(0);
}

//#############------------
//-------- Choose length
//############-----------

int choose_len(int *ldv, double *ldd, int ldlen) {
  int i;
  double tmp=0, ran;
  ran = rand() / ( RAND_MAX + 1.0 );
  for(i=0; i<ldlen; i++){
    if((tmp<=ran) && (ran < tmp + ldd[i])) return(ldv[i]);
    tmp += ldd[i];
  }  
  Rprintf("Error: no length chosen\n");
  return(0);
}

//#############------------
//-------- Choose start
//############-----------

double cumu_fragsta(double x, double *startcdf, double lencdf)
{
  if (x<=0) return 0;
  if (x>=1) return 1;
  int idx= (int) (x*lencdf);
  double y1= startcdf[idx], x1= (double) idx / (lencdf-1);
  idx++;
  double y2= startcdf[idx], x2= (double) idx / (lencdf-1);
  return y1 + (x-x1) * (y2-y1)/(x2-x1);
}

int choose_st(int fraglen, int varlen, double *sdv, double *sdd, int sdlen, int strand){
  int stdlen;
  stdlen = varlen - fraglen + 1;
  if(stdlen < 0) return(-1);
  if(stdlen==0) return(1);
  double maxp=cumu_fragsta((double)stdlen/(double)varlen, sdd, sdlen);  
  double ran = ((double)rand() / (double) RAND_MAX)*maxp;
  return(((int)(cumu_fragsta(ran, sdv, sdlen)*varlen))+1);
}

//#############------------
//------- Build Cigar
//#############------------

void add_match(char *str, int match){
  char tmp[100];
  sprintf(tmp, "%dM", match);
  strcat(str, tmp);
}

void add_gap(char *str, int gap){
  char tmp[100];
  sprintf(tmp, "%dN", gap);
  strcat(str,tmp);
} 

int *build_path(var_t var, int len, int st, int rl, hash_t *path, int strand, int *starts){

  int en, rst, ren, wis, sum;
  char *pa, id[100];
  int i, pos=0, here, skip, l;

  pa = malloc((40 * var.nex) * sizeof(char));
  strcpy(pa, ".");
  if(strand==1) {
    rst = st + len - rl;
    en = st + rl - 1;
    ren = rst + rl;
  }
  else {
    if((var.len - st - len + 2)<0) Rprintf("%d %d %d %d %d %d\n", st, en, rst, ren, var.len, len);    
    st = var.len - st - len + 2;
    rst = st + len - rl; 
    en = st + rl - 1;
    ren = rst + rl - 1;
    if(st<0) Rprintf("%d %d %d %d %d %d\n", st, en, rst, ren, var.len, len);    
}

  here=0;
  sum=1;
  int chk=0;
  for(i=0; i<var.nex; i++) {
    skip=0;
    chk=0;
    wis = abs(var.exen[i] - var.exst[i])+1;
    if((sum<=st) && (st<sum+wis)) {
      sprintf(id, "%d", var.exid[i]);
      strcat(pa, id);
      pos=i;
      here=1;
      skip=1;
    }
    if((sum<=en) && (en<sum+wis)) {
      if(pos!=i){
	strcat(pa, ".");
	sprintf(id, "%d", var.exid[i]);
	strcat(pa, id);
      }
      break;
    }
    if((skip==0) && (here>0)) {
      strcat(pa, ".");
      sprintf(id, "%d", var.exid[i]);
      strcat(pa, id);
    }
    sum+=wis;
  }

  strcat(pa, "-");
  sum=1;
  here=0;
  chk=0;
  for(i=0; i<var.nex; i++) {
    skip=0;
    chk=0;
    wis = abs(var.exen[i] - var.exst[i]) + 1;
    if((sum<=rst) && (rst<sum+wis)){
      sprintf(id, "%d", var.exid[i]);
      strcat(pa, id);
      pos=i;
      here=1;
      skip=1;
    }
    if((sum<=ren) && (ren<sum+wis)) {
      if(pos!=i){
	strcat(pa, ".");
	sprintf(id, "%d", var.exid[i]);
	strcat(pa, id);
      }
      strcat(pa, ".");
      break;
    }
    if((skip==0) && (here>0)) {
      strcat(pa, ".");
      sprintf(id, "%d", var.exid[i]);
      strcat(pa, id);
    }
    sum+=wis;
  }
  starts[2]=0;
  
  l=hash_lookup(path, pa);
  if(l!=HASH_FAIL) hash_update(path, pa, l+1); 
  else {hash_insert(path, pa, 1); starts[2]=1;}
  starts[0] = st;
  starts[1] = rst;

  free(pa);
  return(0);
}

unsigned NextPow2( unsigned x ) {
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return ++x;
}


int *build_cigar(var_t var, int len, int st, int rl, char **cigars, int strand){
  int i, rst, rltmp, *ans, sum=1, wis, done=0;
  ans = malloc(3 * sizeof(int));

  
  if(strand==1){ rst = st + len - rl;} else {st = var.len - st - len + 2; rst = st + len - rl;}
  
  for(i=0; i<var.nex; i++) {
    wis = abs(var.exen[i] - var.exst[i]);
    if((sum<=st) && (st<sum+wis)) {
      st = var.exst[i] + st - sum;
      done++;
      if(done==2) break;
    }
    if((sum<=rst) && (rst < sum+wis)) {
      rst = var.exst[i] + rst - sum;
      done++;
      if(done==2) break;
    }
    sum += wis;
  }

  ans[0] = st;
  ans[1] = rst;
  strcpy(cigars[0], "\0");
  strcpy(cigars[1], "\0");
  
  //Build left read
  sum=1;
  rltmp = rl;
  if(var.len>=rl){
    for(i=0; i<var.nex; i++) {
      if((var.exst[i] <= st) && (st <= var.exen[i])) {
	if(st+rltmp <= var.exen[i]) { add_match(cigars[0], rltmp); rltmp=0; break;} 
	else {
	  if((abs(var.exst[i+1] - var.exen[i]) -1)>0){
	    add_match(cigars[0], abs(var.exen[i] - st ) + 1);// +1 ?
	    rltmp -= abs(var.exen[i] - st) + 1;
	    if(rltmp<=0) break;
	    add_gap(cigars[0], abs(var.exst[i+1] - var.exen[i]) -1); // +1 ?
	    st = var.exst[i+1];
	  } else {
	    if(st+rltmp <= var.exen[i+1]) {
	      add_match(cigars[0], rltmp); rltmp=0; break;}
	    else{
	      i++;
	      if((abs(var.exst[i+1] - var.exen[i]) -1)>0){
		add_match(cigars[0], abs(var.exen[i] - st) + 1);// +
		rltmp -= abs(var.exen[i] - st) + 1;
		if(rltmp<=0) break;
		add_gap(cigars[0], abs(var.exst[i+1] - var.exen[i]) -1); // +1 ?     
		st = var.exst[i+1];
	      } else {
		if(st+rltmp <= var.exen[i+1]) {
		  add_match(cigars[0], rltmp); rltmp=0; break;
		}
		else {
		  i++;
		  if((abs(var.exst[i+1] - var.exen[i]) -1)>0){
		    add_match(cigars[0], abs(var.exen[i] -st) + 1);
		    rltmp -= abs(var.exen[i] - st) + 1;
		    if(rltmp<=0) break;
		    add_gap(cigars[0], abs(var.exst[i+1] - var.exen[i]) -1);
		    st = var.exst[i+1];
		  } else {
		    if(st+rltmp <= var.exen[i+1]) {
		      add_match(cigars[0], rltmp); rltmp=0; break;
		    }
		  }
		}  
	      }
	    }
	  }
	}
    }
    }

    //Build right read
  rltmp = rl;
  sum=1;
  if(rst+rl<=var.exen[var.nex-1] + 1){
  for(i=0; i<var.nex; i++) {
    if((var.exst[i] <= rst) && (rst <= var.exen[i])) {
      if(rst+rltmp <= var.exen[i]) {add_match(cigars[1], rltmp); rltmp=0; break; }
      else {
	if((abs(var.exst[i+1] - var.exen[i]) -1)>0){
	  add_match(cigars[1], abs(var.exen[i] - rst) + 1); 
	  rltmp -= abs(var.exen[i] - rst)+1;
	  if(rltmp<=0) break;
	  add_gap(cigars[1], abs(var.exst[i+1] - var.exen[i]) -1);
	  rst = var.exst[i+1];
	} else {
          if(rst+rltmp <= var.exen[i+1]) {
            add_match(cigars[1], rltmp); rltmp=0; break;
	  }
          else{
	    i++;
	    if((abs(var.exst[i+1] - var.exen[i]) -1)>0){
	      add_match(cigars[1], abs(var.exen[i] - rst) + 1);// +                                                                                                    
	      rltmp -= abs(var.exen[i] - rst)+1;
	      if(rltmp<=0) break;
	      add_gap(cigars[1], abs(var.exst[i+1] - var.exen[i]) -1); // +1 ?                                                                       
	      rst = var.exst[i+1];
	    } else{
	      if(rst+rltmp <= var.exen[i+1]) {
		add_match(cigars[1], rltmp); rltmp=0; break;
	      } else{
		i++;
		if((abs(var.exst[i+1] - var.exen[i]) -1)>0){
		  add_match(cigars[1], abs(var.exen[i] - rst) + 1);
		  rltmp -= abs(var.exen[i] - rst)+1;
		  if(rltmp<=0) break;
		  add_gap(cigars[1], abs(var.exst[i+1] - var.exen[i]) -1);                                                                                           
		  rst = var.exst[i+1];
		} else{
		  if(rst+rltmp <= var.exen[i+1]) {
		    add_match(cigars[1], rltmp); rltmp=0; break;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  } else {
    add_match(cigars[1], rl);
    ans[1]=var.exen[var.nex-1]-rl;    
  }
  }
  else {
    add_match(cigars[0], rl);
    add_match(cigars[1], rl);
    ans[0]=var.exst[0];
    ans[1]=var.exst[0];
  }
  return(ans);
}




