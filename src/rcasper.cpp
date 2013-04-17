#include <sstream>

#include <time.h>

#include <R.h>

#include <Rinternals.h>

#include "stdafx.h"

#include "casper.h"

#include "rcasper.h"

#include "seppel.h"

#include "cppmemory.h"



int verbose = 1;



double cumu_fragsta(double x) 

{

	if (x<=0) return 0;

	if (x>=1) return 1;

	int idx= (int) (x*lencdf);

	double y1= startcdf[idx], x1= (double) idx / (lencdf-1);

	idx++;

	double y2= startcdf[idx], x2= (double) idx / (lencdf-1);

	return y1 + (x-x1) * (y2-y1)/(x2-x1);

}



void importFragments(int np, SEXP pnames, int *pathCounts, int strand, int inv, DataFrame *df){
  
  for (int i = 0; i < np; i++) 
    
    {
      

      //Set counts

      int count = pathCounts[i];
      
      
      //Set nb of left & right visited exons
      
      const char* pname = CHAR(STRING_ELT(pnames, i));

      char* varname = new char[strlen(pname) + 1];
      
      char* left= varname;
      
      strcpy(left, pname);
      
      if (left[0] != '.') continue;
      
      if (left[strlen(pname)-1] != '.') continue;
      
      char* mid = strchr(left, '-');
      
      if (mid == NULL) continue;
      
      mid[0] = '\0';
      
      char* right = mid+1;
      
      int leftc = 0, rightc = 0;
		
      if(strand==-1) {
	
	char *tmp= left;

	left= right;
	
	right= tmp;
	
      }
      
      for (int l = strlen(left)-1; l >= 0; l--) { if (left[l] == '.') leftc++; }  
      
      for (int r = strlen(right)-1; r >= 0; r--) { if (right[r] == '.') rightc++; } 
      
      
      
      Fragment* f = new Fragment(leftc, rightc, count, i);
      

      //Set sequence of visited exons
      
      if ((leftc>0) && (rightc>0)) {

	if(strand==1) {
	  
	  left = left+1;
	  
	  right[strlen(right)-1] = '\0';
	  
	} else {
	  
	  right = right+1;
	  
	  left[strlen(left)-1] = '\0';
	  
	}
	
	char* item;
	
	

	item = strtok(left, ".");  //split string

	for (int j = 0; item != NULL; j++) {
	  
	  int eid = atoi(item);
	  
	  if(strand==1) f->left[j] = eid;
	  
	  else f->left[leftc-j-1] = eid;
	  
	  item = strtok(NULL, ".");
	  
	}

	item = strtok(right, ".");
	
	for (int j = 0; item != NULL; j++) {
	  
	  int eid = atoi(item);
			  
	  if(strand==1) f->right[j] = eid;
	  
	  else f->right[rightc-j-1] = eid;
	  
	  item = strtok(NULL, ".");

	}


	
	bool c1= (strand==1) && (f->left[0] <= f->right[0]) && (f->left[f->leftc -1] <= f->right[f->rightc -1]);

	bool c2= (strand== -1) && (f->left[0] >= f->right[0]) && (f->left[f->leftc -1] >= f->right[f->rightc -1]);
	
	if (c1 || c2) {
	  
	  if(inv==0) {

	    df->addData(f);

	  } else {

	    df->addDataM(f);

	  }
	
	} else {
	  
	  delete f;
	  
	}
      } 

      delete [] varname;
    }



}


DataFrame* importDataFrame(SEXP exonsR, SEXP exonwidthR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP strandR)

{

  int nexons=Rf_length(exonsR), nfraglen=Rf_length(fraglenR), readLength=INTEGER(readLengthR)[0], strand=INTEGER(strandR)[0];

  int *exons=INTEGER(exonsR), *exonwidth=INTEGER(exonwidthR), *lenvals=INTEGER(lenvalsR);

  double *fraglen= REAL(fraglenR);



  //Define fragment length/start distributions and create DataFrame                                                                                                                                                                            

  DiscreteDF* fraglen_dist = new DiscreteDF(fraglen, lenvals, nfraglen);



  lencdf= Rf_length(fragstaR);

  startcdf= REAL(fragstaR);

  //::fun_fragsta = fragstaR;                                                                                                                                                                                                                  



  DataFrame* df = new DataFrame(fraglen_dist, cumu_fragsta);

  df->frag_readlen= readLength;



  //Add exons to DataFrame                                                                                                                                                                                                                     

  for (int i = 0; i < nexons; i++)

    {

      Exon *ex= new Exon(exons[i], exonwidth[i]);

      df->addExon(ex);

    }



  //Add exon path counts                                                                                                                                                                                                                       

  int np = LENGTH(pathCountsR);

  SEXP pnames = getAttrib(pathCountsR, R_NamesSymbol);

  int* pathCounts = INTEGER(pathCountsR);


  if(strand!=0) {

    importFragments(np, pnames, pathCounts, strand, 0, df);

  } else {

    importFragments(np, pnames, pathCounts, 1, 0, df);

    importFragments(np, pnames, pathCounts, -1, 1, df);

  }



  return df;


}



void importTranscripts(set<Variant*, VariantCmp> *initvars, DataFrame* df, SEXP transcriptsR, SEXP strandR)

{

	int nt = LENGTH(transcriptsR);

	SEXP tnames = getAttrib(transcriptsR, R_NamesSymbol);

	Variant *v;

	SEXP trow;

	int ntsub, *tvals;

	int strand=INTEGER(strandR)[0];

	for (int i = 0; i < nt; i++) {

		trow = VECTOR_ELT(transcriptsR, i);

		ntsub = LENGTH(trow);

		tvals = INTEGER(trow);

		

		vector<Exon*>* el = new vector<Exon*>();

		for (int s = 0; s < ntsub; s++) {

                  int eid = tvals[s];

		  Exon* ex = df->id2exon[eid];

		  el->push_back(ex);

                }



		v = new Variant(el);

		v->id= i;


		if((tvals[0] > tvals[1]) && strand==0) v->antisense=TRUE;

		int nbchar= Rf_length(STRING_ELT(tnames,i));

		v->name= string(CHAR(STRING_ELT(tnames, i)), nbchar);

		initvars->insert(v);

		delete el;

	}



}



extern "C"

{



  //Evaluate log(likelihood(pi)) + log(prior(pi)) for a grid of pi values in the known variants case

  //Also returns posterior mode, Hessian at posterior mode & path probabilities for each variant

  // Input

  // - gridR: matrix with nb rows = nb of variants, nb cols = nb grid points

  // (other args as in calcKnownSingle)

  SEXP lhoodGrid(SEXP gridR, SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP priorqR, SEXP strandR) {



    //CREATE CASPER INSTANCE

    DataFrame* df = importDataFrame(exonsR, exonwidthR, pathCountsR, fragstaR, fraglenR, lenvalsR, readLengthR, strandR);



    set<Variant*, VariantCmp> *initvars = new set<Variant*, VariantCmp>();

    importTranscripts(initvars, df, transcriptsR, strandR);

    df->fixUnexplFrags(initvars, 0); // Discard fragments that are unexplained by know variants

    double priorq = REAL(priorqR)[0];

    int totC=0;

    list<Fragment*>::const_iterator fi;

    for (fi = df->data.begin(); fi != df->data.end(); fi++)

      {

	Fragment* f = *fi;

	totC += f->count;

      }



    Model* model = new Model(initvars);

    Casper* casp = new Casper(model, df);

    Casper::priorq = priorq;

    Casper::em_maxruns = 1000;

    Casper::em_tol= 0.00001;



    int vc = model->count();

    double* em = casp->calculateMode();



    int nrow= Rf_nrows(gridR), ncol = Rf_ncols(gridR);



    SEXP ans;

    PROTECT(ans= allocVector(VECSXP, 7));



    //RETURN LOG-LIKELIHOOD + LOG-PRIOR

    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP,ncol));

    double *logpos= REAL(VECTOR_ELT(ans,0));

    gridR = coerceVector(gridR, REALSXP); //coerce matrix to vector

    double *pivals= REAL(gridR);

    for (int i=0; i<ncol; i++) logpos[i]= casp->priorLikelihoodLn(pivals + nrow*i);



    //RETURN POSTERIOR MODE & HESSIAN

    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP,vc));  //stores estimated expression

    SET_VECTOR_ELT(ans, 2, allocVector(STRSXP,vc)); //stores variant names

    SET_VECTOR_ELT(ans, 3, allocVector(REALSXP,1)); //stores logpos at the mode

                

    double *expr= REAL(VECTOR_ELT(ans,1));

    SEXP varnamesR= VECTOR_ELT(ans,2);              		

    double *fopt= REAL(VECTOR_ELT(ans,3));



    for (int j=0; j< vc; j++) {

      Variant* v = model->get(j);

      int varidx= model->indexOf(v);

      if(totC>0) expr[j] = em[varidx]; //estimated expression

      else expr[j] = 0;

      if (initvars->count(v)>0) v->name= (*initvars->find(v))->name;  //respect initial variant names

      const char *cname= (v->name).c_str();

      SET_STRING_ELT(varnamesR,j,mkChar(cname));  //variant name

    }

    fopt[0]= casp->priorLikelihoodLn(em);



    SET_VECTOR_ELT(ans, 4, allocMatrix(REALSXP,vc-1,vc-1)); //stores variance of estimated expression (logit scale)

    double *Svec= REAL(VECTOR_ELT(ans,4));

    double **Sinv= dmatrix(1,vc,1,vc), **S= dmatrix(1,vc,1,vc);

    if(totC>0){

      casp->normapprox(Sinv, em, vc, 1);

      bool posdef;

      inv_posdef(Sinv,vc-1,S,&posdef);

      for (int i=1; i<vc; i++) Svec[i-1 + (i-1)*(vc-1)]= S[i][i];

      for (int i=1; i<vc; i++) for (int j=i+1; j<vc; j++) Svec[j-1 + (i-1)*(vc-1)]= Svec[i-1 + (j-1)*(vc-1)]= S[i][j];

    }

    free_dmatrix(Sinv,1,vc,1,vc); free_dmatrix(S,1,vc,1,vc);



    //RETURN PATH PROBABILITIES FOR EACH VARIANT

    int np = (casp->frame->data).size();

    SET_VECTOR_ELT(ans, 5, allocMatrix(REALSXP,vc,np));  //stores estimated expression

    double *probmatrix = REAL(VECTOR_ELT(ans, 5));

    for (int i=0; i<vc; i++) {

      Variant *v= casp->model->get(i);

      map<Fragment*, double> vprobs= df->probabilities(v);

      int j=0;

      for (list<Fragment*>::iterator fi= df->data.begin(); fi != df->data.end(); fi++) {

	Fragment *f= *fi;

	probmatrix[i+vc*j]= vprobs[f];

	j++;

      }

    }



    SEXP dimnames, rownames;

    PROTECT(dimnames = allocVector(VECSXP, 2));

    PROTECT(rownames = allocVector(STRSXP, np));

    int j=0;

    for (list<Fragment*>::iterator fi= df->data.begin(); fi != df->data.end(); fi++) {

      Fragment *f= *fi;

      std::stringstream sstr;

      if (INTEGER(strandR)[0]==1) {

        for (int i=0; i < f->leftc; i++) { sstr << "."; sstr << f->left[i]; }

        sstr << "-";

        for (int i=0; i < f->rightc; i++) { sstr << f->right[i]; sstr << "."; }

      } else {

        for (int i=(f->rightc-1); i >= 0; i--) { sstr << "."; sstr << f->right[i]; }

        sstr << "-";

        for (int i=(f->leftc-1); i >= 0; i--) { sstr << f->left[i]; sstr << "."; }

      }

      std::string str1= sstr.str();

      const char *cname= str1.c_str();

      SET_STRING_ELT(rownames,j,mkChar(cname));  //path name



      j++;

    }

    SET_VECTOR_ELT(dimnames, 1, rownames);

    setAttrib(VECTOR_ELT(ans,5), R_DimNamesSymbol, dimnames);



    //RETURN INTEGRATED LIKELIHOOD

    SET_VECTOR_ELT(ans, 6, allocVector(REALSXP,2));

    Casper::is_runs= 10000;

    REAL(VECTOR_ELT(ans,6))[0]= casp->calculateIntegral(1);

    REAL(VECTOR_ELT(ans,6))[1]= casp->calculateIntegral(2);



    UNPROTECT(3);

    return ans;

  }



  //Estimate isoform expression for multiple genes. Calls calcKnownSingle repeadtedly

  //Input args same as for calcDenovoMultiple

  SEXP calcKnownMultiple(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP geneidR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP priorqR, SEXP strandR, SEXP citypeR, SEXP niterR, SEXP burninR, SEXP verboseR) 

	{

		int i, ngenes=LENGTH(geneidR);

		SEXP ansMultiple, ansSingle;

		double paccept=0;



		PROTECT(ansMultiple= allocVector(VECSXP, ngenes));



		for (i=0; i<ngenes; i++) {

		  ansSingle= calcKnownSingle(&paccept, VECTOR_ELT(exonsR,i), VECTOR_ELT(exonwidthR,i), VECTOR_ELT(transcriptsR,i), VECTOR_ELT(pathCountsR,i), fragstaR, fraglenR, lenvalsR, readLengthR, priorqR, VECTOR_ELT(strandR, i), citypeR, niterR, burninR);

		  SET_VECTOR_ELT(ansMultiple,i,ansSingle);

		}

		if (INTEGER(verboseR)[0]==1 && INTEGER(citypeR)[0]==2) Rprintf("Average MH acceptance rate %f\n",paccept/(ngenes+.0));

		

		UNPROTECT(1);

		return ansMultiple;

	}



  //Estimate isoform expression (known variants case)

  //Input args same as for calcDenovoSingle

  // - citypeR: set to 0 to return only estimated expression. 1: additionally, return asymptotic variances. 2: additionally, return posterior samples

  // - niterR: number of iterations to obtain posterior samples

  // - burninR: number of burnin iterations

  SEXP calcKnownSingle(double *paccept, SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP priorqR, SEXP strandR, SEXP citypeR, SEXP niterR, SEXP burninR)

	{

		DataFrame* df = importDataFrame(exonsR, exonwidthR, pathCountsR, fragstaR, fraglenR, lenvalsR, readLengthR, strandR);


		set<Variant*, VariantCmp> *initvars = new set<Variant*, VariantCmp>();

		importTranscripts(initvars, df, transcriptsR, strandR);

		df->fixUnexplFrags(initvars, 0); // Discard fragments that are unexplained by know variants

		double priorq = REAL(priorqR)[0];

		int totC=0;

		list<Fragment*>::const_iterator fi;

		Fragment* f;

		for (fi = df->data.begin(); fi != df->data.end(); fi++)

		  {

		    f = *fi;

		    totC++;

		    break;

		  }

		if(INTEGER(strandR)[0]==0) {

		  for (fi = df->dataM.begin(); fi != df->dataM.end(); fi++)
  
		    {
  
		      f = *fi;

		      totC++;

		      break;

		    }

		}


		//Model* model = new Model(new vector<Variant*>(initvars->begin(), initvars->end()));

		Model* model = new Model(initvars);

		Casper* casp = new Casper(model, df);

		Casper::priorq = priorq; //		casp->priorq = priorq;

		Casper::em_maxruns = 1000;

		Casper::em_tol= 0.00001;

		totC = casp->totCounts();

		int vc = model->count();

		double* em = casp->calculateMode();



         	SEXP ans;

                PROTECT(ans= allocVector(VECSXP, 5));



  		SET_VECTOR_ELT(ans, 0, allocVector(REALSXP,vc));  //stores estimated expression

                SET_VECTOR_ELT(ans, 1, allocVector(STRSXP,vc)); //stores variant names	   

		

                

		double *expr= REAL(VECTOR_ELT(ans,0));

                SEXP varnamesR= VECTOR_ELT(ans,1);              		

		

	 	for (int j=0; j< vc; j++) {

		  Variant* v = model->get(j);

		  int varidx= model->indexOf(v);

		  if(totC>0) expr[j] = em[varidx]; //estimated expression

		  else expr[j] = 0;

		  if (initvars->count(v)>0) v->name= (*initvars->find(v))->name;  //respect initial variant names

		  const char *cname= (v->name).c_str();

		  SET_STRING_ELT(varnamesR,j,mkChar(cname));  //variant name

                }



		if (INTEGER(citypeR)[0]>0) {

		  SET_VECTOR_ELT(ans, 2, allocVector(REALSXP,vc)); //stores variance of estimated expression (logit scale)

		  double *vexpr= REAL(VECTOR_ELT(ans,2));



		  if (INTEGER(citypeR)[0]==1) {

		    if (totC>0) casp->asymptoticSE(vexpr, em, vc); else { for (int j=0; j<vc; j++) vexpr[j] = 0; }

		  } else if (INTEGER(citypeR)[0]==2) {

		    if (totC>0) {

		      double **S= dmatrix(1,vc,1,vc), **Sinv= dmatrix(1,vc,1,vc);

		      casp->normapprox(Sinv, em, vc, 1);

		      bool posdef;

		      inv_posdef(Sinv,vc-1,S,&posdef);

		      int niter= INTEGER(niterR)[0], burnin= INTEGER(burninR)[0];

		      SET_VECTOR_ELT(ans, 3, allocVector(REALSXP,vc*(niter-burnin))); //stores posterior samples

		      double *pi = REAL(VECTOR_ELT(ans,3));

		      double pacc, integralIS;

		      casp->IPMH(pi, &pacc, &integralIS, niter, burnin, em, Sinv);

		      (*paccept) += pacc;

		      free_dmatrix(S,1,vc,1,vc); free_dmatrix(Sinv,1,vc,1,vc);

		    } else {

		      SET_VECTOR_ELT(ans, 3, allocVector(REALSXP,vc));

		      double *pi = REAL(VECTOR_ELT(ans,3));

		      for(int j=0; j<vc; j++) pi[j]=0;

		    }

		  }

		}


		SET_VECTOR_ELT(ans, 4, allocVector(INTSXP, 1));
		
		int *totCp = INTEGER(VECTOR_ELT(ans, 4));

		totCp[0] = totC;
		
		UNPROTECT(1);

		delete df;

		delete initvars;

		delete model;

		delete casp;

		zaparray(em);



		return ans;

	}



  SEXP calcDenovoMultiple(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP geneidR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP modelUnifPriorR, SEXP nvarPriorR, SEXP nexonPriorR, SEXP priorqR, SEXP minppR, SEXP selectBest, SEXP methodR, SEXP niterR, SEXP exactMarginalR, SEXP verboseR, SEXP strandR) 

	{

	  //De novo isoform discovery and estimate expression for multiple genes. Calls calcDenovoSingle repeadtedly

	  int i, ngenes=LENGTH(geneidR), ngenes10, verbose= INTEGER(verboseR)[0];

	  SEXP ansMultiple, ansSingle;



		PROTECT(ansMultiple= allocVector(VECSXP, ngenes));



		if (ngenes>10) ngenes10= ngenes/10; else ngenes10=1;		

		for (i=0; i<ngenes; i++) {

		  int nexons = min(LENGTH(VECTOR_ELT(exonsR,i)), LENGTH(nvarPriorR));

		  ansSingle= calcDenovoSingle(VECTOR_ELT(exonsR,i), VECTOR_ELT(exonwidthR,i), VECTOR_ELT(transcriptsR,i), VECTOR_ELT(pathCountsR,i), fragstaR, fraglenR, lenvalsR, readLengthR, modelUnifPriorR, VECTOR_ELT(nvarPriorR,nexons-1), VECTOR_ELT(nexonPriorR,nexons-1), priorqR, minppR, selectBest, methodR, VECTOR_ELT(niterR,i), exactMarginalR, VECTOR_ELT(strandR, i));

		  SET_VECTOR_ELT(ansMultiple,i,ansSingle);

		  if (verbose && (i%ngenes10)==0) Rprintf(".");

		}

		

		UNPROTECT(1);

		return ansMultiple;

	}



  SEXP calcDenovoSingle(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP modelUnifPriorR, SEXP nvarPriorR, SEXP nexonPriorR, SEXP priorqR, SEXP minppR, SEXP selectBest, SEXP methodR, SEXP niterR, SEXP exactMarginalR, SEXP strandR)

	{

	//De novo isoform discovery and estimate expression for a single gene

	//Input

	// - exons: vector with exon ids

	// - exonwidth: vector exon widths

	// - transcripts: list of transcripts. Names indicate transcript id. Each element contains list of exon ids. Used to initialize model

	// - pathCounts: vector with path counts. Names indicate series of visited exons e.g. e1.e2-e5

	// - fragsta: function that returns start distrib cdf

	// - fraglen: vector with fragment length distrib, i.e. P(length=lenvals[0]),P(length=lenvals[1]),... up to max length

	// - lenvals: vector with possibles values for length

	// - readLength: read length

        // - modelUnifPrior: when set to FALSE use a uniform prior on the model space (nvarPrior & nexonPrior are ignored)

	// - nvarPrior: vector with NegBinom parameters for prior prob of nb expressed variants

	// - nexonPrior: vector with Beta-Binomial param for prior prob of nb exons in a variant

	// - priorq: prior on model-specific parameters pi ~ Dirichlet(priorq)

	// - minpp: only models with post prob >= minpp are reported

	// - selectBest: set to !=0 to return only results for model with highest posterior probability

	// - methodR: 1 is exhaustive enumeration; 2 random walk MCMC; 3 prior proposal MCMC; 0 auto (exhaustive for <=4 exons, RW-MCMC otherwise)

	// - niter: number of mcmc iterations

	// - exactMarginal: set to 0 to estimate post prob from proportion of MCMC visits; set to 1 to use marginal likelihoods of MCMC visited models
 
   

	  int  selBest = INTEGER(selectBest)[0], method= INTEGER(methodR)[0], niter= INTEGER(niterR)[0];

	  int exactMarginal= INTEGER(exactMarginalR)[0];

	  double minpp = REAL(minppR)[0], priorq = REAL(priorqR)[0];



	  DataFrame* df = importDataFrame(exonsR, exonwidthR, pathCountsR, fragstaR, fraglenR, lenvalsR, readLengthR, strandR); //read input



	  set<Variant*, VariantCmp> *initvars = new set<Variant*, VariantCmp>();

	  importTranscripts(initvars, df, transcriptsR, strandR);



	  //df->debugprint();  //debug

	  //Model* tmpm = new Model(initvars);  //debug

	  //tmpm->debugprint(); //debug



	  int discarded = df->fixUnexplFrags(initvars, 1); //add variants to initvars (initial model)

	  //Model* tmp2 = new Model(initvars);  //debug

	  //tmp2->debugprint();  //debug



	  Casper::priorq = priorq;

	  Casper::em_maxruns = 100;

	  Casper::em_tol= 0.001;



	  map<Model*, double, ModelCmp> resProbs;

	  map<Model*, double*, ModelCmp> resModes;



	  Seppel* seppl;

	  if (INTEGER(modelUnifPriorR)[0]) {

	    seppl= new Seppel(df);

	  } else {

	    seppl = new Seppel(df, REAL(nvarPriorR), REAL(nexonPriorR));

	  }



		if (method == 1 || method == 0 && df->exons.size() <= 4) 

		{

			seppl->exploreExact();
			
			resProbs = seppl->resultPPIntegral();

		} 

		else if (method == 3)

		{

			seppl->exploreUnif(niter);

			if (exactMarginal) {

			  resProbs= seppl->resultPPIntegral();

			} else {

			  resProbs = seppl->resultPPMCMC();

			}

		}

		else

		{

			Model* startmodel = new Model(initvars);

			seppl->exploreSmart(startmodel, niter);

			if (exactMarginal) {

			  resProbs= seppl->resultPPIntegral();

			} else {

			  resProbs = seppl->resultPPMCMC();

			}

		}

		resModes = seppl->resultModes();

		

		// END OF CALCULATIONS



		

		Model* bestModel;

		double bestModelProb = -1, bestModelPrior = -1;



		map<Model*, double, ModelCmp>::iterator mvi = resProbs.begin();

		while (mvi != resProbs.end())

		{

			if (mvi->second > bestModelProb)

			{

				bestModelProb = mvi->second;

				bestModel = mvi->first;

				bestModelPrior = seppl->calculatePrior(bestModel);

			}



			if (mvi->second < minpp) 

			{

				resModes.erase(mvi->first);

				resProbs.erase(mvi++);

			}

			else

			{

				++mvi;

			}

		}



		// END OF FILTERING

		

		SEXP ans;

		PROTECT(ans= allocVector(VECSXP, 8));



		//Report posterior probabilities

		int i=0;

		int nrowpi=0, nx;

		if (selBest==0) { nx= resProbs.size(); } else { nx=1; }

		SET_VECTOR_ELT(ans, 0, allocMatrix(REALSXP,nx,3));

		SET_VECTOR_ELT(ans, 5, allocVector(STRSXP,nx));

		double *resProbsR= REAL(VECTOR_ELT(ans,0));

		

		set<Variant*, VariantCmp>::const_iterator vi;

		map<Model*, double, ModelCmp>::const_iterator mi;

		if (selBest==0) {

			for (mi = resProbs.begin(), i=0; mi != resProbs.end(); mi++, i++) {

				Model* m = mi->first;

				resProbsR[i]= i;

				resProbsR[i+nx]= resProbs[m];

				resProbsR[i+2*nx]= exp(seppl->calculatePrior(m));

				nrowpi+= m->count();

				//SET_STRING_ELT(VECTOR_ELT(ans,5),i,mkChar(m->getCodeStr(allpossvariants)));  //model name

			}

		} else {

			resProbsR[0]= 0;

			resProbsR[nx]= resProbs[bestModel];

			resProbsR[2*nx]= exp(seppl->calculatePrior(bestModel));

			nrowpi+= bestModel->count();

			//SET_STRING_ELT(VECTOR_ELT(ans,5),0,mkChar(bestModel->getCodeStr(allpossvariants)));  //model name

		}



		//Report estimated expression

		SET_VECTOR_ELT(ans, 1, allocMatrix(REALSXP,nrowpi,2));  //stores model id and estimated expression

		SET_VECTOR_ELT(ans, 2, allocVector(STRSXP,nrowpi)); //stores variant names

		double *expr= REAL(VECTOR_ELT(ans,1));

		SEXP varnamesR= VECTOR_ELT(ans,2);



		int rowid=0, nrowexons=0;

		set<Variant*, VariantCmp> allvariants;

		if (selBest==0) {

			for (mi = resProbs.begin(), i=0; mi != resProbs.end(); mi++, i++) {

				Model* m = mi->first;

				for (int j=0; j< m->count(); j++) {

					expr[rowid]= i;  //model id

					Variant* v = m->get(j);

					int varidx= m->indexOf(v);

					expr[rowid+nrowpi]= resModes[m][varidx]; //estimated expression

					if (initvars->count(v)>0) v->name= (*initvars->find(v))->name;  //respect initial variant names

					const char *cname= (v->name).c_str();;

					SET_STRING_ELT(varnamesR,rowid,mkChar(cname));  //variant name

					if (allvariants.count(v)==0) {

						allvariants.insert(v);

						nrowexons+= v->exonCount;

					}

					rowid++;

				}

			}

		} else {

			Model* m = bestModel;

			for (int j=0; j< m->count(); j++) {

				expr[rowid]= i;  //model id

				Variant* v = m->get(j);

				int varidx= m->indexOf(v);

			        expr[rowid+nrowpi]= resModes[m][varidx]; //estimated expression

				if (initvars->count(v)>0) v->name= (*initvars->find(v))->name;  //respect initial variant names

				const char *cname= (v->name).c_str();

				SET_STRING_ELT(varnamesR,rowid,mkChar(cname));  //variant name

				if (allvariants.count(v)==0) {

					allvariants.insert(v);

					nrowexons+= v->exonCount;

				}

				rowid++;

			}

		}



		//Report exons in each variant

		SET_VECTOR_ELT(ans, 3, allocVector(INTSXP,nrowexons));  //stores exon id

		SET_VECTOR_ELT(ans, 4, allocVector(STRSXP,nrowexons));   //stores variant name

		int *ranges= INTEGER(VECTOR_ELT(ans,3));

		SEXP varnames2R= VECTOR_ELT(ans,4);



		rowid=0;

		for (vi= allvariants.begin(); vi != allvariants.end(); vi++) {

			Variant *v= *vi;

			for (int j=0; j< v->exonCount; j++) {

				ranges[rowid]= v->exons[j]->id ; //exon id

				const char *cname= (v->name).c_str();

				SET_STRING_ELT(varnames2R,rowid,mkChar(cname));  //variant name

				rowid++;

			}

		}



		//Report sum (integrated likelihood * prior)/exp(integralMax) across all models

		SET_VECTOR_ELT(ans,6, allocVector(REALSXP,2));

		REAL(VECTOR_ELT(ans,6))[0]= seppl->integralSum;

		REAL(VECTOR_ELT(ans,6))[1]= seppl->integralMax;



		SET_VECTOR_ELT(ans,7, allocVector(INTSXP,1));

		INTEGER(VECTOR_ELT(ans,7))[0]= discarded;



		UNPROTECT(1);

		delete seppl;

		delete df;

		delete initvars;



		return(ans);

	}



}

