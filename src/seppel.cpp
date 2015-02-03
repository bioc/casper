#include <algorithm>
#include <iostream>
#include <limits>
#include "cppmemory.h"
#include "dropVariant.h"
#include "seppel.h"
using namespace std;


Seppel::Seppel(DataFrame* frame, set<Variant*>* knownVars, int integrateMethod)

{

  this->frame = frame;

  this->modelUnifPrior= 1;

  this->knownVars= knownVars;

  this->integrateMethod= integrateMethod;


  varis = new vector<Variant*>();

  varisSet = new set<Variant*>();

  models = new vector<Model*>();

  modelsSet = new set<Model*>();

}



Seppel::Seppel(DataFrame* frame, set<Variant*>* knownVars, double* nvarPrior, double* nexonPrior, double* prioradj, int integrateMethod) {

  int E= frame->exons.size();
  double a=nexonPrior[0], b=nexonPrior[1], tmp;

  this->frame = frame;
  this->modelUnifPrior= 0;
  this->knownVars= knownVars;
  this->integrateMethod= integrateMethod;

  varis = new vector<Variant*>();
  varisSet = new set<Variant*>();
  models = new vector<Model*>();
  modelsSet = new set<Model*>();


  //Prior on number of exons in a variant

//  double offset= 0.01, dbar= prioradj[1];
//  if (dbar < E-1.0) { //if known vars discard >1 exon (on the average), adjust prior mode to dbar/E
// 
//    if ((a+b)<(2.0+offset)) {
//      double pa <- a/(a+b);
//      double c1 <- 2.0+offset-a-b;
//      a= a+c1*pa;
//      b= b+c1*(1-pa);
//    }
//    c2= (dbar/E)*(a+b-2) -a +1;
//    if (c2 > b-offset) { c2= b-offset; } //ensure parameters >0  
//    a <- a+c2; b <- b-c2;
//  }

  tmp = 0;
  for (int i=1; i<= E; i++) {
    tmp = lnbeta(a +i-1, b +E-1 -(i-1)) - lnbeta(a,b);
    priorpNbExons.push_back( exp(tmp) * (i+.0)/(E+.0) );
    nvarsPoibin.push_back((int) (choose(E,i)+.1));
  }


  //Prior on number of variants in the model (truncated to <=1000)
  int imax= (int) min(pow(2, E) -1, 1000.0);
  double psum=0, prob= 1-nvarPrior[0]; //success probability (from R we pass failure prob)

  tmp = 0;

  for (int i=1; i<= imax; i++) {
    tmp = dnegbinomial(i, nvarPrior[1], prob, 0);
    priorpNbVars.push_back( log(tmp) );
    psum += tmp;
  }

  psum = log(psum);

  for (int i=0; i< imax; i++) priorpNbVars[i] -= psum;

}



Seppel::~Seppel() {



  //Free all considered variants & models

  std::for_each(varis->begin(), varis->end(), DeleteFromVector());

  varis->clear();



  std::for_each(models->begin(), models->end(), DeleteFromVector());

  models->clear();



  std::for_each(varisSet->begin(), varisSet->end(), DeleteFromVector());

  varisSet->clear();



  std::for_each(modelsSet->begin(), modelsSet->end(), DeleteFromVector());

  modelsSet->clear();



  //Free modes for each model

  map<Model*, double*, ModelCmp>::const_iterator it;

  for (it = modes.begin(); it != modes.end(); it++) { delete [] it->second; }



}





double* Seppel::initMode(Model* model, Model* similarModel) {



    int n= model->count(), nSimilar= similarModel->count(), ncommon=0;

    double sumcommon=0;

    double* pi = new double[n];

    double* piSimilar= modes[similarModel];



    for (int i=0; i<n; i++) pi[i]= 0;



    for (int i=0; i<nSimilar; i++) {

      Variant* v= similarModel->get(i);

      if (model->contains(v)) {

	double elem= piSimilar[i];

	//double elem= piSimilar[similarModel->indexOf(v)];

	pi[model->indexOf(v)]= elem;

	sumcommon += elem;

	ncommon++;

      }

    }

    double pcommon= ((double) ncommon)/((double) n);

    double norm= pcommon/sumcommon;

    if (ncommon==n) {

      for (int i=0; i<n; i++) pi[i] *= norm;

    } else {

      double newpi= (1.0-pcommon)/((double) (n-ncommon));

      for (int i=0; i<n; i++) {

        if (pi[i]> 0) pi[i] *= norm; else pi[i]= newpi;

      }

    }

    return pi;

}


double Seppel::calcIntegral(Model* model, Model* similarModel) {

  return calcIntegral(model, similarModel, true);

}


double Seppel::calcIntegral(Model* model, Model* similarModel, bool knownVarsCheck) 

{

  double *mode;

  if (model == NULL) return 1;

  if (integrals.count(model) > 0) return integrals[model];

  if (modes.count(similarModel)==0) return this->calcIntegral(model);

  unsigned int nknownVars= knownVars->size();

  if (knownVarsCheck && (nknownVars > 0)) {

    unsigned int knownInModel=0;
    vector<Variant *>::const_iterator itvarvec;

    itvarvec= (model->items).begin();

    while ((knownInModel<nknownVars) && (itvarvec != (model->items).end())) {

      knownInModel += knownVars->count(*itvarvec);

      itvarvec++;

    }

    if (knownInModel < nknownVars) {

      integrals[model]= 1;
      return 1;

    }

  }

  double like = 1, prior = 1;

  Casper* casp = new Casper(model, frame);



  if (casp->isValid()) {



    mode = this->initMode(model,similarModel);

    casp->calculateMode(mode);

    modes[model] = mode;

    like = casp->calculateIntegral(mode, model->count(), this->integrateMethod);

    prior = calculatePrior(model);

    like += prior;



  }

  integrals[model] = like;

  priorprobs[model] = prior;



  delete casp;

  return like;

}


double Seppel::calcIntegral(Model* model) {

  return calcIntegral(model, true);

}


double Seppel::calcIntegral(Model* model, bool knownVarsCheck) 

{

  if (model == NULL) return 1;

  if (integrals.count(model) > 0) return integrals[model];


  unsigned int nknownVars= knownVars->size();

  if (knownVarsCheck && (nknownVars > 0)) {

    unsigned int knownInModel=0;
    vector<Variant *>::const_iterator itvarvec;

    itvarvec= (model->items).begin();

    while ((knownInModel<nknownVars) && (itvarvec != (model->items).end())) {

      knownInModel += knownVars->count(*itvarvec);

      itvarvec++;

    }

    if (knownInModel < nknownVars) {

      integrals[model]= 1;
      return 1;

    }

  }



  double like = 1, prior = 1;

  Casper* casp = new Casper(model, frame);



  if (casp->isValid()) {

    double* mode = casp->calculateMode();

    modes[model] = mode;

    like = casp->calculateIntegral(mode, model->count(), this->integrateMethod);

    prior = calculatePrior(model);

    like += prior;

  }

  integrals[model] = like;

  priorprobs[model] = prior;



  delete casp;

  return like;

}


void Seppel::exploreExactFast(set<Variant*, VariantCmp> *initvaris) {

  int maxdrops, k, nmodels;

  frame->allVariants(varis, initvaris);

  Model* m= new Model(varis);

  modelsSet->insert(m);

  //int maxdropit=1;
  //if (m->items.size() > 1) maxdropit = (int) ceil(log2((double) m->items.size()) + 1);

  //Choose maxdrops such that max models <10,000,000, i.e. sum_{k= nvar-maxdrops}^{nvar} lchoose(nvar, k) < 10,000,000
  nmodels= 0; k= (varis->size())+1;
  while ((k>=1) && (nmodels < 10000000)) {
    k--;
    nmodels += choose((double) varis->size(), (double) k);
  }
  maxdrops= (varis->size() - k);

  exploreSubmodels(m, maxdrops); //exhaustively consider submodels of m by dropping up to maxdrops variables

}


void Seppel::exploreExact(set<Variant*, VariantCmp> *initvaris) {

  frame->allModels(varis, models, initvaris);  //store all possible variants & models



	vector<Model*>::const_iterator mi;

	for (mi = models->begin(); mi != models->end(); mi++)

	{

		Model* model = *mi;

		calcIntegral(model);

	}

}

void Seppel::exploreUnif(int runs, set<Variant*, VariantCmp> *initvaris) {



  vector<Variant*>* varis = new vector<Variant*>();

  vector<Model*>* models = new vector<Model*>();

  frame->allModels(varis, models, initvaris);  //store all possible variants & models



  vector<Model*>* possiblemodels = new vector<Model*>();

	vector<Model*>::const_iterator ami;

	for (ami = models->begin(); ami != models->end(); ami++)

	{

		Casper* ncasp = new Casper(*ami, frame);

		if (ncasp->isValid())

		{

			possiblemodels->push_back(ncasp->model);

			counts[ncasp->model] = 0;

		}

		delete ncasp;

	}

	if (possiblemodels->size() == 0)

	{

		return;

	}



	int onum = runifdisc(0, possiblemodels->size() - 1);

	Model* omodl = possiblemodels->at(onum);

	double olike = calcIntegral(omodl);
	


	int accepted = 0;

	

	for (int r = 0; r < runs; r++)

	{

		int nnum = runifdisc(0, possiblemodels->size() - 1);

		Model* nmodl = possiblemodels->at(nnum);



		double nlike = calcIntegral(nmodl);

		if (nlike != 1)

		{

			double l = nlike - olike;

			double p = exp(l);

			double x = runif();

			if (x <= p)

			{

				olike = nlike;

				omodl = nmodl;

				accepted++;

			}

		}

		counts[omodl]++;

	}

	delete [] possiblemodels;

}

void Seppel::exploreSmart(Model* startmodel, int runs)

{



  modelsSet->insert(startmodel);

  //models->push_back(startmodel);

	Model *omodl = startmodel;

	double olike = calcIntegral(omodl);
	
	SmartModelDist* odist = new SmartModelDist(this, frame, omodl, 0.8, modelsSet);

	

	int accepted = 0;
	double nlike, nprob, oprob, l, lp, x;



	for (int r = 0; r < runs; r++)

	{

	  Model* nmodl = odist->sample(varisSet);  //propose new model, add variants to varisSet

	  modelsSet->insert(nmodl);

	  //models->push_back(nmodl);

	  nlike = calcIntegral(nmodl,omodl);

	  //double nlike = calcIntegral(nmodl);



		if (nlike != 1)

		{

		  SmartModelDist* ndist = new SmartModelDist(this, frame, nmodl, 0.8, modelsSet);  //create proposal, add considered models to models



			nprob = odist->densityLn(nmodl);

			oprob = ndist->densityLn(omodl);

			l = nlike - olike + oprob - nprob;

			lp = exp(l);

			x = runif();

			if (x < lp) {

			  omodl = nmodl;

			  delete odist;

			  odist = ndist;

			  olike = nlike;

			  accepted++;

			} else {

			  delete ndist;

			}

		}



		double count = 0;

		if (counts.count(omodl) > 0)

		{

			count = counts[omodl];

		}

		counts[omodl] = count + 1;	

	}

	delete odist;



}


void Seppel::exploreSubmodels(Model* model, int maxdrops) {
 //Exhaustively consider all submodels of a given model obtained by dropping up to maxdrops variants from the model
 //Uses a recursive scheme so that for any submodel with 0 prob, no further submodels are considered.
  double nlike;
  vector<Variant*> dropvars (model->items);
  vector<Variant*> keepvars;

  modelsSet->insert(model);
  nlike= calcIntegral(model);
  if (nlike != 1) {
    exploreSubmodels(model, &keepvars, &dropvars, false, 0, maxdrops);
  }

}

void Seppel::exploreSubmodels(Model* model, vector<Variant*>* keepvars, vector<Variant*>* dropvars, bool eval_pp, int dropcount, int maxdrops) {
  //- model: current model, used to speed up calculations by giving starting value for posterior mode of submodels
  //- keepvars: these variants are always kept in the model
  //- dropvars: variants to be recursively dropped from the model
  //- eval_pp: if set to true the posterior probability of the input model is evaluated
  //- dropcount: number of variants that have been dropped in previous recursive calls
  //- maxdrops: when dropcount >= maxdrops the algorithm stops
  bool stopnow= false;

  if (eval_pp) {
    double pp;
    Model* newmodel;
    //Create model with (keepvars,dropvars)
    vector<Variant*> newitems;
    newitems.reserve(keepvars->size() + dropvars->size());
    newitems.insert(newitems.end(), keepvars->begin(), keepvars->end());
    newitems.insert(newitems.end(), dropvars->begin(), dropvars->end());
    newmodel= new Model(&newitems);
    modelsSet->insert(newmodel);
    pp= calcIntegral(newmodel, model, false);
    model= newmodel;
    if (pp==1) { stopnow= true; }
  }

  if ((dropvars->size() == 0) || (dropvars->size()==1 && keepvars->size()==0) || (dropcount >= maxdrops)) {
    stopnow= true;
  }

  if (!stopnow) {
    //dnew= dropvars[-1]; knew= c(keepvars,dropvars[1])
    vector<Variant*> dnew, knew;
    dnew.reserve(dropvars->size() - 1);
    dnew.insert(dnew.end(), (dropvars->begin())+1, dropvars->end());
    knew.reserve(keepvars->size()+1);
    knew.insert(knew.end(), keepvars->begin(), keepvars->end());
    knew.insert(knew.end(), dropvars->begin(), dropvars->begin()+1);
    exploreSubmodels(model, keepvars, &dnew, true, dropcount+1, maxdrops); //models without Variant 1
    exploreSubmodels(model, &knew, &dnew, false, dropcount, maxdrops);  //models with Variant 1
  }

}

/*
OLD VERSION OF exploreSubmodels. INEFFICIENT AS IT REQUIRED ENUMERATING ONE SAME MODEL MANY TIMES

void Seppel::exploreSubmodelsOld(Model* model, int maxdropit, int maxmodels) {
 //Exhaustively consider all submodels of a given model, up to a certain depth controlled by maxdropit.
 //Uses an iterative scheme so that for any submodel with 0 prob, no further submodels are considered.
 //Input
 // - model: model for which we want to consider all possible submodels
 // - maxdropit: max number of iterations in the scheme. For maxdropit=1 submodels dropping 1 variant are considered, else up to 2^(maxdropit-1) variants are dropped
 // - maxmodels: if at any iterations more than maxmodels would be considered, the iterative scheme stops. Defaults to 2^20

  int i, j, nvars= model->count();
  double nlike;
  Variant *v;
  vector<Variant*> newitems (model->items);
  Model *newmodel;
  dropVariant *dropvars= new dropVariant(nvars), *tmp;

  modelsSet->insert(model);
  nlike= calcIntegral(model);

  //Initialize: check dropping individual variants
  for (i=0; i<nvars; i++) {

    v= newitems[0];                    //pointer to variant i
    newitems.erase(newitems.begin());  //remove variant i

    if (knownVars->count(v) == 0) {

      newmodel= new Model(&newitems);  //new model without variant i
      modelsSet->insert(newmodel);
   
      nlike = calcIntegral(newmodel,model,false);
      if (nlike != 1) {
        int *varsin;
        varsin= ivector(0,nvars-1);
        for (j=0; j<i; j++) varsin[j]= 1;
        varsin[i]= 0;
        for (j=i+1; j<nvars; j++) varsin[j]= 1;
        dropvars->add(varsin);
      }

    }

    newitems.push_back(v); //add back variant i

  }

  //Iterate: drop combinations of variables
  i=0; int ncomb= dropvars->size() * (dropvars->size() - 1) / 2;

  while ((i<maxdropit) && (dropvars->size()>1) && (ncomb <= maxmodels) ) {
    tmp= dropvars->combinations();
    delete dropvars;
    dropvars= tmp;
    map<string, int*>::const_iterator mi;
    mi= dropvars->submodels.begin();

    while(mi != dropvars->submodels.end()) {
      //Define new model
      newitems.clear();
      for (j=0; j<nvars; j++) {
        if (mi->second[j] == 1) newitems.push_back(model->get(j));
      }
      newmodel= new Model(&newitems); //new model
      modelsSet->insert(newmodel); //store model

      //If new model has 0 prob, erase from dropvars
      nlike = calcIntegral(newmodel,model); //post prob of newmodel
      if (nlike == 1) {
	string s= mi->first;
	mi++;
	dropvars->erase(s);
      } else {
	mi++;
      }
    }

    i++;
    ncomb= dropvars->size() * (dropvars->size() - 1) / 2;

  }

  delete dropvars;
}
*/


map<Model*, double*, ModelCmp> Seppel::resultModes()

{

	return modes;

}

map<Model*, double, ModelCmp> Seppel::resultPPIntegral()

{

	map<Model*, double, ModelCmp> probs;



	integralMax = -DBL_MAX;

	

	map<Model*, double, ModelCmp>::const_iterator mi;

	for (mi = integrals.begin(); mi != integrals.end(); mi++)

	{

		if (mi->second == 1)

		{

			continue;

		}


		integralMax = max(mi->second, integralMax);

	}



	integralSum = 0;

	for (mi = integrals.begin(); mi != integrals.end(); mi++)

	{

		if (mi->second == 1)

		{

			continue;

		}

		integralSum += exp(mi->second - integralMax);

	}

	double lsum = integralMax + log(integralSum);



	for (mi = integrals.begin(); mi != integrals.end(); mi++) 

	{

		if (mi->second == 1)

		{

			continue;

		}

		probs[mi->first] = exp(mi->second - lsum);

	}



	return probs;

}

map<Model*, double, ModelCmp> Seppel::resultPPMCMC()

{

	map<Model*, double, ModelCmp> probs;

	map<Model*, double, ModelCmp>::const_iterator mi;

	double total = 0;

	for (mi = counts.begin(); mi != counts.end(); mi++)

	{

		total = total + mi->second;

	}



	for (mi = counts.begin(); mi != counts.end(); mi++)

	{

		if (mi->second > 0)

		{

			probs[mi->first] = mi->second / total;

		}

	}

	return probs;

}



void Seppel::normalizeIntegrals(double *probs, double *values, int n)

{

	double imax = -DBL_MAX;

	

	for (int i = 0; i < n; i++)

	{

		imax = max(values[i], imax);

	}

	

	double asum = 0;

	for (int i = 0; i < n; i++)

	{

		asum += exp(values[i] - imax);

	}



	double lsum = imax + log(asum);

	

	//double* probs = new double[n];

	double psum = 0;

	for (int i = 0; i < n; i++)

	{

		probs[i] = exp(values[i] - lsum);

		psum += probs[i];

	}

	for (int i = 0; i < n; i++)

	{

		probs[i] = probs[i] / psum;

	}

	//return probs;

}





double Seppel::calculatePrior(Model* model) {

  if (priorprobs.count(model) > 0) return priorprobs[model];

  if (modelUnifPrior==1) {

    return 0.0;

  } else {

    int E= frame->exons.size(), nbVars= model->count();

    if (nbVars > (int) priorpNbVars.size()) {  //nb variants > min(2^E -1, 1000) has 0 prob

      return -std::numeric_limits<double>::infinity();

    } else {

      double ans;

      //Prior on nb variants in the model

      ans= priorpNbVars[nbVars-1];



      //Prior on nb exons per variant

      if (nbVars < pow(2,E) -1) {  //when all variants were selected, prob of selected exons is 1 so this computation is skipped

        vector<int> Sk (E,0);

        for (int i=0; i< nbVars; i++) {  

          Variant* v= model->get(i); 
	  
          Sk[v->exonCount -1]++;

        }

        vector<int> Fk (E);

        for (int i=0; i< E; i++) Fk[i]= nvarsPoibin[i] - Sk[i];



        for (int i=0; i< E; i++) {

          ans += Sk[i] * log(priorpNbExons[i]) + Fk[i] * log(1-priorpNbExons[i]);
          //printf("i=%d, %f, %f\n", i, Sk[i] * log(priorpNbExons[i]), Fk[i] * log(1-priorpNbExons[i]));

        }



        if (E <= 20) {

	  ans -= dpoissonbin(model->count(), &priorpNbExons, &nvarsPoibin, 1, &Tvector, &poibinProbs); //Poisson-Binomial

        } else {

	  ans -= dpoisson(model->count(), 1.0, 1); //Poisson approx (law of small numbers)

        }

      }



      return ans;



    }



  }

}



