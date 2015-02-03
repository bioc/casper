#include <algorithm>
#include "dataframe.h"
#include "cppmemory.h"
#include <list>
#include <Rinternals.h>
using namespace std;


DataFrame::DataFrame(DiscreteDF* fraglen_dist, double (*fragsta_cumu)(double))

{

	this->fraglen_dist = fraglen_dist;

	this->fragsta_cumu = fragsta_cumu;

	this->frag_readlen = 75;



        fraglen_minx= fraglen_dist->value(0);

        fraglen_maxx= fraglen_dist->value((fraglen_dist->size)-1);

	//	for (fraglen_minx = 0; fraglen_dist->probability(fraglen_minx) == 0; fraglen_minx++) ;

	//	for (fraglen_maxx = fraglen_dist->size - 1; fraglen_dist->probability(fraglen_maxx) == 0; fraglen_maxx--) ;

}



DataFrame::~DataFrame() {

  std::for_each(exons.begin(), exons.end(), DeleteFromVector());

  std::for_each(data.begin(), data.end(), DeleteFromVector());

  std::for_each(dataM.begin(), dataM.end(), DeleteFromVector());

  //FreeClear( &exons );

  //FreeClear( &data );

  delete fraglen_dist;

}



void DataFrame::addData(Fragment* f)

{

	this->data.push_back(f);

}

void DataFrame::addDataM(Fragment* f)

{

	this->dataM.push_back(f);

}

void DataFrame::addExon(Exon* e)

{

	e->num = this->exons.size();

	this->id2exon[e->id] = e;

	this->exons.push_back(e);

}



map<Fragment*, double> DataFrame::probabilities(Variant* v)

{

	if (this->cache.count(v) > 0)

	{

		return this->cache[v];

	}



	list<Fragment*>::const_iterator fi;


	if(v->antisense){

	  
	  for (fi = dataM.begin(); fi != dataM.end(); fi++)
	    
	    {
	      
	      Fragment* f = *fi;

	      double p = probability(v, f);

	      if (p > 0.0)
  
		{
		 
		  cache[v][f] = p;

                }

	    }
	
	} else {


	  for (fi = data.begin(); fi != data.end(); fi++)
	    
	    {
	      
	      Fragment* f = *fi;
	      
	      double p = probability(v, f);

	      if (p > 0.0)
		
		{

		  cache[v][f] = p;
		  
		}

	    }
	}

	return cache[v];

}

double DataFrame::probability(Variant* v, Fragment* f)

{

	double p = 0.0;


	if (v->contains(f))
	  
	  {
	    //printf("Not discarded %d %d %d %d %d %d\n", f->left[0], f->left[f->leftc - 1], f->right[0], f->right[f->rightc - 1], v->exons[0]->id, v->exons[v->exonCount-1]->id);	    
	    int fs = v->indexOf(f->left[0]);
	    
	    int fe = v->indexOf(f->left[f->leftc - 1]);
	    
	    int bs = v->indexOf(f->right[0]);
	    
	    int be = v->indexOf(f->right[f->rightc - 1]);

	    p = prob(fs, fe, bs, be, v->positions, v->length);
	    //printf("next %d %d %d %d %d %d %0.12f\n", f->left[0], f->left[f->leftc - 1], f->right[0], f->right[f->rightc - 1], v->exons[0]->id, v->exons[v->exonCount-1]->id, p);	    
	      
	  } //else printf("discarded %d %d %d %d %d %d\n", f->left[0], f->left[f->leftc - 1], f->right[0], f->right[f->rightc - 1], v->exons[0]->id, v->exons[v->exonCount-1]->id);

	return p;

}



double DataFrame::prob(int fs, int fe, int bs, int be, int* pos, double T) {

  int vals, readLength; 
  double probs;
  DiscreteDF* lengthDist;

  if (T <= fraglen_minx) {  //if shortest fragment longer than current transcript
    vals= (int) T; 
    probs=1;
    lengthDist= new DiscreteDF(&probs,&vals,1);
  } else {
    lengthDist= fraglen_dist;
  }

  if ((frag_readlen) > ((int) T)) { //if left+right read don't fit in current transcript
    readLength= (int) (T-1);
  } else {
    readLength= frag_readlen;
  }

  // lower bound for start of left transcript
  double a1 = max(pos[fs], pos[fe] - readLength + 1);
  // upper bound for start of left transcript
  double b1 = min(pos[fs + 1] - 1, pos[fe + 1] - readLength);
  // lower bound for end of right transcript
  double a2 = max(pos[bs] + readLength, pos[be] + 1);
  // upper bound for end of right transcript
  double b2 = min(pos[bs + 1] + readLength - 1, pos[be + 1]);

  double psum = 0;

  for (int i=0; i< lengthDist->size; i++) { //stop before T

    double l= lengthDist->value(i);
    double mb= (T-l+1.0)/T;
    double rb = min(min(b1, b2 - l) / T, mb);
    double lb = min((max(a1, a2 - l) - 1.0) / T, mb);

    if (lb >= rb) { continue; } //if( i== (lengthDist->size -1)) {printf("%f %f %f %f %f %f %f %f %f %d\n", lb, rb, a1, b1, a2, b2, l, T, mb, readLength); }continue; }
      
    double punc = (fragsta_cumu(rb) - fragsta_cumu(lb)) / fragsta_cumu(mb);
    double factor = 0;

    if (l <= T && punc > 0) {
      factor = lengthDist->probability(i);
      if ((T < fraglen_maxx) && (T > fraglen_minx)) { factor /= lengthDist->cumulativeProbability((int)(T-fraglen_minx)); }
    }
    //printf(" %0.12f %0.12f %d\n", factor, psum, i);

    psum += punc * factor;

  }

  if (T <= fraglen_minx) { 
    delete lengthDist; 
  }

  return psum;

}



int DataFrame::fixUnexplFrags(set<Variant*, VariantCmp>* initvars, std::map<Variant*,std::string>* varshortnames, int* geneid, int denovo) {

	// copy all fragments
	set<Fragment*>* queue = new set<Fragment*>(data.begin(), data.end());
	set<Fragment*>::iterator itqueue;

	// remove the fragments from the queue we can explain with our variants
	set<Variant*, VariantCmp>::iterator vi;

	for (vi = initvars->begin(); vi != initvars->end(); vi++) {

		map<Fragment*, double> probs = probabilities(*vi);

		map<Fragment*, double>::iterator si;

		for (si = probs.begin(); si != probs.end(); si++) {

		  set<Fragment*>::iterator ri = queue->find(si->first);

		  if (ri != queue->end()) queue->erase(ri);

		}

	}


	int discarded= 0;

	if (denovo) {

	  //Propose new variants
	  set <string> allvarnames;
	  set <Variant*, VariantCmp> newvaris, bestvaris, *varisptr;
	  set <Variant*, VariantCmp>::const_iterator itvarset;
	   
	  for (itvarset = initvars->begin(); itvarset != initvars->end(); itvarset++) allvarnames.insert((*itvarset)->exoncomb);
	   
	  for (itqueue = queue->begin(); itqueue != queue->end(); itqueue++) {
	   
	    bool explained = false;
	    Fragment* frag = (*itqueue);
	   
	    path2Variants(&newvaris, &bestvaris, &allvarnames, &explained, initvars, frag); //propose new variants, add their names to allvarnames

	    if (!explained) {
	      discarded++;
	      data.remove(frag);
	    }
 
	  }

	  if ((initvars->size() + newvaris.size()) > 20) {  //if too many new variants were proposed, used only best ones

	    varisptr= &bestvaris;

	    for (itvarset = newvaris.begin(); itvarset != newvaris.end(); itvarset++) {

	      if (bestvaris.count(*itvarset) == 0) delete (*itvarset);  //delete variants that will not be used

	    }

	  } else {

	    varisptr= &newvaris;

	  }

	  //Add new variants to initvars, short names to varshortnames
	  for (itvarset = varisptr->begin(); itvarset != varisptr->end(); itvarset++) {
     
	    Variant* nv = (*itvarset);

	    initvars->insert(nv);

	    std::ostringstream out;
	    out << "CASP.";
	    out << geneid[0];
	    out << ".";
	    out << (varshortnames->size()+1);
	    (*varshortnames)[nv] += out.str();

	  }

	} else {  //not denovo

	    // discard all unexplained fragments by known variants in known case

	    while (queue->size() > 0) {

	      Fragment* frag = *queue->begin();

              queue->erase(queue->begin());

	      data.remove(frag);

	      discarded++;

	    }

	}

	delete queue;

	return discarded;

}


/* Old version: it only proposed one variant per fragment
int DataFrame::fixUnexplFrags(set<Variant*, VariantCmp>* initvars, std::map<Variant*,std::string>* varshortnames, int* geneid, int denovo) {

	// copy all fragments

	set<Fragment*>* queue = new set<Fragment*>(data.begin(), data.end());


	// remove the fragments from the queue we can explain with our variants

	set<Variant*, VariantCmp>::iterator vi;

	for (vi = initvars->begin(); vi != initvars->end(); vi++) {

		map<Fragment*, double> probs = probabilities(*vi);

		map<Fragment*, double>::iterator si;

		for (si = probs.begin(); si != probs.end(); si++) {

		  set<Fragment*>::iterator ri = queue->find(si->first);

		  if (ri != queue->end()) queue->erase(ri);

		}

	}



	int discarded = 0;

	if(denovo){

	  while (queue->size() > 0) {  // while we still have unexplained fragments

	      // Pop the first fragment and propose variants

	      Fragment* frag = *queue->begin();

	      queue->erase(queue->begin());
     
	      Variant* nv = path2Variant(frag);


	      // Check if the new variant can explain the fragment

	      map<Fragment*, double> probs = probabilities(nv);

	      if (probs.count(frag) > 0) {

		  initvars->insert(nv);

		  std::ostringstream out;
		  out << "CASP.";
		  out << geneid[0];
		  out << ".";
		  out << (varshortnames->size()+1);
		  (*varshortnames)[nv] += out.str();


			// delete all fragments that this variant can explain

		  map<Fragment*, double>::iterator si;

		  for (si = probs.begin(); si != probs.end(); si++) {

		      set<Fragment*>::iterator ri = queue->find(si->first);

		      if (ri != queue->end()) queue->erase(ri);

		    }

		} else {

			// this fragment cant be explained

		  discarded++;

		  data.remove(frag);

		}

	  }

	} else {  //not denovo

	    // discard all unexplained fragments by know variants in known case

	    while (queue->size() > 0) {

	      Fragment* frag = *queue->begin();

              queue->erase(queue->begin());

	      data.remove(frag);

	      discarded++;

	    }

	  }

	delete queue;

	return discarded;

}
*/



Variant* DataFrame::path2Variant(Fragment* f) 

{

	int eid=0; Exon *ex;

	vector<Exon*>::iterator itexon;

	vector<Exon*>* el = new vector<Exon*>();

	for (itexon= exons.begin(); (*itexon)->id != f->left[0]; itexon++) {
	  
		ex= (*itexon);
		
		el->push_back(ex);

	}

	for (int i=0; i< f->leftc; i++) {

		eid = f->left[i];

		ex = id2exon[eid];

		el->push_back(ex);

	}

	if (eid != f->right[0]) {

		eid = f->right[0];

		ex = id2exon[eid];

		el->push_back(ex);

	}

	for (int i=1; i< f->rightc; i++) {

		eid = f->right[i];

		ex = id2exon[eid];

		el->push_back(ex);

	}
	
	while ((*itexon)->id != eid) { itexon++; }

	itexon++;

	while (itexon != exons.end()) {

		ex= (*itexon);

		el->push_back(ex);

		itexon++;

	}

	Variant* v = new Variant(el);

	delete el;

	return v;

}



//For each variant in initvaris containing all exons in the path set by f, add a new variant to newvaris
//
// Input
// - initvaris: initial set of variants that will act as templates for the new variants
// - f: fragment we wish to explain
//
// Ouput
// - newvaris: set where elements will be added (can be non-empty at entry)
// - bestvaris: variant assigning highest prob to is added to bestvar
// - allvarnames: new variants with exoncomb in allvarnames are discarded. allvarnames is updated with new variants in newvaris
// - explained: returns true if a new variant assigns prob >0 to frag, else returns false
//
// Note: as a side effect, path probabilities for new variants are computed & stored in this->cache
void DataFrame::path2Variants(set<Variant*, VariantCmp> *newvaris, set<Variant*, VariantCmp> *bestvaris, set <string> *allvarnames, bool *explained, set <Variant*, VariantCmp> *initvaris, Fragment* f) {

  (*explained) = false;

  double maxprob=0;
  Variant *bestvar;

  set <Variant*, VariantCmp>::iterator itvarset;

  for (itvarset = initvaris->begin(); itvarset != initvaris->end(); itvarset++) {

    bool extendvar=true;
    int eid=0, curvarex, strand=1; Exon *ex; Variant *curvar= (*itvarset);

    //Check that variant contains all exons in the path
    for (int i=0; (i<(f->leftc)) && extendvar; i++) {
      ex= id2exon[f->left[i]];
      if (!curvar->contains(ex)) { extendvar= false; }
    }

    for (int i=0; (i<(f->rightc)) && extendvar; i++) {
      ex= id2exon[f->right[i]];
      if (!curvar->contains(ex)) { extendvar= false; }
    }

    if (extendvar) {

      vector<Exon*>::iterator itexon;
      vector<Exon*>* el = new vector<Exon*>();
      curvarex = 0; ex= (curvar->exons)[0];
       
      if(((f->leftc)>1) & (f->left[0]>f->left[1])) strand=-1;
      if(((f->rightc)>1) & (f->right[0]>f->right[1])) strand=-1;
      if((f->left[0] > f->right[0]) | (ex->id > (curvar->exons)[curvar->exonCount-1]->id)) strand=-1;
       
      if(strand==1){
        while ((curvarex < curvar->exonCount) && (ex->id < (f->left[0]))) {
       	el->push_back(ex);
       	curvarex++;
       	ex= (curvar->exons)[curvarex];
        }
       
        for (int i=0; i< f->leftc; i++) {
       	eid = f->left[i];
       	if (id2exon.count(eid)>0) { ex = id2exon[eid]; } else { Rf_error("Exon %d in path counts not found in genomeDB!\n",eid); }
       	el->push_back(ex);
        }
        
        for (int i=0; i< f->rightc; i++) {
       	eid = f->right[i];
       	if (id2exon.count(eid)>0) { ex = id2exon[eid]; } else { Rf_error("Exon %d in path counts not found in genomeDB!\n",eid); }
       	if (eid > f->left[f->leftc -1]) el->push_back(ex);
        }
        
        if (eid < ((curvar->exons)[curvar->exonCount -1])->id) {
       	ex= (curvar->exons)[curvarex];
       	while (ex->id <= eid) { curvarex++; ex= (curvar->exons)[curvarex]; }
       	while ( curvarex < curvar->exonCount ) {
       	  el->push_back(ex);
       	  curvarex++;
       	  ex= (curvar->exons)[curvarex];
       	}
        } 
      } else {
       
        while ((curvarex < curvar->exonCount) && (ex->id > (f->left[0]))) {
       	el->push_back(ex);
       	curvarex++;
       	ex= (curvar->exons)[curvarex];
        }
       
        for (int i=0; i< f->leftc; i++) {
       	eid = f->left[i];
       	if (id2exon.count(eid)>0) { ex = id2exon[eid]; } else { Rf_error("Exon %d in path counts not found in genomeDB!\n",eid); }
       	el->push_back(ex);
        }
       
       
        for (int i=0; i< f->rightc; i++) {
       	eid = f->right[i];
       	if (id2exon.count(eid)>0) { ex = id2exon[eid]; } else { Rf_error("Exon %d in path counts not found in genomeDB!\n",eid); }
          if (eid > f->left[f->leftc -1]) el->push_back(ex);
        }
       
        if (eid > ((curvar->exons)[curvar->exonCount -1])->id) {
          ex= (curvar->exons)[curvarex];
          while (ex->id >= eid) { curvarex++; ex= (curvar->exons)[curvarex]; }
          while ( curvarex < curvar->exonCount ) {
            el->push_back(ex);
            curvarex++;
            ex= (curvar->exons)[curvarex];
          }
        }
       
      }
       	
      if (el->size() > 0) {
       
        Variant* v = new Variant(el);
       
        if (allvarnames->count(v->exoncomb) == 0) {
	  allvarnames->insert(v->exoncomb);
	  double fragprob = probability(v, f);
	  if (fragprob > 0) {
	    *explained = true;
	    newvaris->insert(v);
	    map<Fragment*, double> probs = probabilities(v);  //update cache with all path prob for new variant
	    if (fragprob > maxprob) { maxprob= fragprob; bestvar= v; }
	  } else {
	    delete v;
	  }
        } else {
	  delete v;
        }
       
      }
       
      delete el;
    }  //end if (extendvar)
  }

  if (*explained) bestvaris->insert(bestvar);

}





void DataFrame::allModelsRec(vector<Variant*>* stack, unsigned int level, vector<Variant*>* vars, vector<Model*>* models)

{

	if (vars->size() == level)

	{

		if (stack->size() > 0)

		{

			Model* m = new Model(stack);

			models->push_back(m);

		}

		return;

	}



	stack->push_back(vars->at(level));

	allModelsRec(stack, level + 1, vars, models);

	stack->pop_back();	

	allModelsRec(stack, level + 1, vars, models);

}

void DataFrame::allModels(vector<Variant*> *varis, vector<Model*> *models, set<Variant*, VariantCmp> *initvaris) {

  set<string> inithash;
  set<Variant*, VariantCmp>::iterator itvec;

  for (itvec = initvaris->begin(); itvec != initvaris->end(); itvec++) {

    inithash.insert((*itvec)->exoncomb);

  }


  vector<Exon*>* estack = new vector<Exon*>();

  allVariantsRec(estack, 0, varis, &inithash);

  for (itvec = initvaris->begin(); itvec != initvaris->end(); itvec++) varis->push_back(*itvec);
	

  vector<Variant*>* vstack = new vector<Variant*>();

  allModelsRec(vstack, 0, varis, models);



  delete estack; 

  delete vstack; 

}


void DataFrame::allModels(vector<Variant*> *varis, vector<Model*> *models, vector<Variant*> *initvaris) {

  set<string> inithash;
  vector<Variant*>::iterator itvec;

  for (itvec = initvaris->begin(); itvec != initvaris->end(); itvec++) {

    inithash.insert((*itvec)->exoncomb);

  }



  vector<Exon*>* estack = new vector<Exon*>();

  allVariantsRec(estack, 0, varis, &inithash);

  for (itvec = initvaris->begin(); itvec != initvaris->end(); itvec++) varis->push_back(*itvec);
	

  vector<Variant*>* vstack = new vector<Variant*>();

  allModelsRec(vstack, 0, varis, models);



  delete estack; 

  delete vstack; 

}


//Old version not considering initial variant names
//void DataFrame::allModels(vector<Variant*> *varis, vector<Model*> *models) {
// 
//  vector<Exon*>* estack = new vector<Exon*>();
// 
//  allVariantsRec(estack, 0, varis);
// 
// 	
// 
//  vector<Variant*>* vstack = new vector<Variant*>();
// 
//  allModelsRec(vstack, 0, varis, models);
// 
// 
// 
//  delete estack; 
// 
//  delete vstack; 
// 
//}


void DataFrame::allVariantsRec(vector<Exon*>* stack, unsigned int level, vector<Variant*>* varis, set<string>* inithash) {

  if (exons.size() == level) {

    if (stack->size() > 0)	{

      Variant* v = new Variant(stack);

      if (inithash->count(v->exoncomb) == 0) {

	varis->push_back(v);

      } else {

	delete v;

      }

    }

    return;

  }

  stack->push_back(exons.at(level));

  allVariantsRec(stack, level + 1, varis, inithash);

  stack->pop_back();

  allVariantsRec(stack, level + 1, varis, inithash);

}



void DataFrame::allVariants(vector<Variant*> *varis, set<Variant*, VariantCmp> *initvaris) {

  set<string> inithash;
  set<Variant*, VariantCmp>::iterator itvec;

  for (itvec = initvaris->begin(); itvec != initvaris->end(); itvec++) {

    inithash.insert((*itvec)->exoncomb);

  }

  vector<Exon*>* estack = new vector<Exon*>();

  allVariantsRec(estack, 0, varis, &inithash);

  for (itvec = initvaris->begin(); itvec != initvaris->end(); itvec++) varis->push_back(*itvec);

  delete estack;

}

void DataFrame::allVariants(vector<Variant*> *varis, vector<Variant*> *initvaris) {

  set<string> inithash;
  vector<Variant*>::iterator itvec;

  for (itvec = initvaris->begin(); itvec != initvaris->end(); itvec++) {

    inithash.insert((*itvec)->exoncomb);

  }

  vector<Exon*>* estack = new vector<Exon*>();

  allVariantsRec(estack, 0, varis, &inithash);

  for (itvec = initvaris->begin(); itvec != initvaris->end(); itvec++) varis->push_back(*itvec);

  delete estack;

}

//Old version, does not respect already existing variants

//void DataFrame::allVariantsRec(vector<Exon*>* stack, unsigned int level, vector<Variant*>* varis) {
// 
//  if (exons.size() == level) {
// 
//    if (stack->size() > 0)	{
// 
//      Variant* v = new Variant(stack);
// 
//      varis->push_back(v);
// 
//    }
// 
//    return;
// 
//  }
// 
//  stack->push_back(exons.at(level));
// 
//  allVariantsRec(stack, level + 1, varis);
// 
//  stack->pop_back();
// 
//  allVariantsRec(stack, level + 1, varis);
// 
//}

//void DataFrame::allVariants(vector<Variant*> *varis) {
// 
//  vector<Exon*>* estack = new vector<Exon*>();
// 
//  allVariantsRec(estack, 0, varis);
// 
//  delete estack;
// 
//}


/*
void DataFrame::debugprint() {		
	// Exons
	Rprintf("Exons:\n");
	vector<Exon*>::const_iterator ei;
	for (ei = exons.begin(); ei != exons.end(); ei++) {
		Exon* e = *ei;
		Rprintf("%i\t%i\n", e->id, e->length);
	}
	Rprintf("\n");

	// Fragments
	Rprintf("Fragments:\n");
	list<Fragment*>::const_iterator fi;
	for (fi = data.begin(); fi != data.end(); fi++)	{
		Fragment* f = *fi;
		Rprintf("%i\t%i\t%i\n", f->leftc, f->rightc, f->count);
		for (int l = 0; l < f->leftc; l++) Rprintf("%i\n", f->left[l]);
		for (int r = 0; r < f->rightc; r++) Rprintf("%i\n", f->right[r]);
	}

	Rprintf("\n");

}
*/


bool compareF(Fragment* first, Fragment* second)

{ return( first->id == second->id  ); }


bool orderF(Fragment* first, Fragment* second)

{ return( first->id < second->id  ); }


int DataFrame::totCounts()

{

  list<Fragment*> data;

  data = this->data;

  if(this->dataM.size() > 0) 

    {

      data.insert(data.end(), this->dataM.begin(), this->dataM.end());

      data.sort(orderF);

      data.unique(compareF);

    }


  Fragment* f;

  int totC = 0;

  list<Fragment*>::const_iterator fi;                                                                                                                                                                                             

  for (fi = data.begin(); fi != data.end(); fi++) {                                                                                                                                                                               


    f = *fi;                                                                                                                                                                                                      
    
    totC += f->count;

  }

  return(totC);

}
