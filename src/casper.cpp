#include "casper.h"

#include <algorithm>



using namespace std;



int Casper::is_runs = 10000;

int Casper::em_maxruns = 100;

double Casper::em_tol = 0.001;

const double Casper::mh_gammah = 2;

double Casper::priorq = 3;



Casper::Casper(Model* model, DataFrame* frame)

{

	this->model = model;

	this->frame = frame;



	vector<Variant* >::const_iterator vi;

	for (vi = model->items.begin(); vi != model->items.end(); vi++)

	{

	  Variant *v = *vi;

		map<Fragment*, double> table = this->frame->probabilities(v);



		memvprobs[v] = table;

        

		map<Fragment*, double>::const_iterator fi;

		for (fi = table.begin(); fi != table.end(); fi++)

		{

			Fragment* f = fi->first;

			mempprobs[f][v] = fi->second;

		}

	}

}



double* Casper::calculateMode() {

  int n = model->count();

  double* pi = new double[n];

  for (int i = 0; i < n; i++) pi[i] = 1.0 / (double)n;

  this->calculateMode(pi);

  return pi;

}



void Casper::calculateMode(double* pi) {



  double normali = (double)memvprobs.size() * (priorq - 1.0);

  map<Fragment*, map<Variant*, double> >::const_iterator ofi;

  for (ofi = mempprobs.begin(); ofi != mempprobs.end(); ofi++) {

    normali += ofi->first->count;

  }



  double err= 1.0, newpi;

  int r;

  for (r = 0; (r < em_maxruns) & (err>em_tol); r++) {

    map<Fragment*, double> mem = fragdist(pi);

		

    map<Variant*, map<Fragment*, double> >::const_iterator vi;

    for (vi = memvprobs.begin(), err=0; vi != memvprobs.end(); vi++) {

      int i = model->indexOf(vi->first);



      double nsum = 0;

      map<Fragment*, double>::const_iterator fi;

      for (fi = vi->second.begin(); fi != vi->second.end(); fi++) {

	nsum += (double)fi->first->count * fi->second / mem[fi->first];

      }



      newpi = (priorq - 1.0 + nsum * pi[i]) / normali;

      err= max_xy(err, fabs(newpi - pi[i]));

      pi[i]= newpi;

    }

  }

}



void Casper::IPMH(double *pi, double *paccept, double *integralIS, int niter, int burnin) {

  double* mode = calculateMode();

  IPMH(pi, paccept, integralIS, niter, burnin, mode);

  delete [] mode;

}



void Casper::IPMH(double *pi, double *paccept, double *integralIS, int niter, int burnin, double *mode) {

  double **Sinv;

  int n = model->count();

  Sinv= dmatrix(1,n,1,n);

  normapprox(Sinv, mode, n, 1);

  IPMH(pi, paccept, integralIS, niter, burnin, mode, Sinv);

  free_dmatrix(Sinv,1,n,1,n);

}



void Casper::IPMH(double *pi, double *paccept, double *integralIS, int niter, int burnin, double *mode, double **Sinv) {

  //Note: input parameter S is the Hessian at the mode, i.e. the inverse of the covariance matrix

  double det, lold, lnew, l0;

  double *thmode, *thold, *thnew, *piold, *pinew, **Gold, **Gnew, **cholS, **cholSinv;

  int n = model->count();



  //Pre-compute useful quantities

  thmode = new double[n - 1];

  mlogit(thmode, mode, n);



  cholS= dmatrix(1,n-1,1,n-1);

  cholSinv= dmatrix(1,n-1,1,n-1);

  bool posdef;

  choldc(Sinv,n-1,cholSinv,&posdef);

  if (!posdef) {

    int i; double lmin=0, *vals;

    vals= dvector(1,n);

    eigenvals(Sinv,n-1,vals);

    for (i=1; i<n; i++) if (vals[i]<lmin) lmin= vals[i];

    lmin = -lmin + .001;

    for (i=1; i<n; i++) Sinv[i][i] += lmin;

    choldc(Sinv,n-1,cholSinv,&posdef);

    free_dvector(vals,1,n);

  }

  choldc_inv(Sinv,n-1,cholS,&posdef); 

  det= choldc_det(cholSinv,n-1);



  //MCMC

  thold = new double[n - 1];

  thnew = new double[n - 1];

  piold = new double[n];

  pinew = new double[n];

  Gold = dmatrix(0,n,0,n);

  Gnew = dmatrix(0,n,0,n);



  rmvtC(thold-1, n-1, thmode-1, cholS, 3);

  milogit(piold, thold, n);

  l0= lold= priorLikelihoodLn(piold) - dmvtC(thold-1, n-1, thmode-1, cholSinv, det, 3, 1);



  vtGradG(Gold,thold, n);

  lold+= vtGradLogdet(Gold, n);



  (*integralIS)= 0;

  (*paccept)= 0;

  for (int i=0; i<niter; i++) {

    rmvtC(thnew-1, n-1, thmode-1, cholS, 3);

    milogit(pinew, thnew, n);

    lnew= priorLikelihoodLn(pinew) - dmvtC(thnew-1, n-1, thmode-1, cholSinv, det, 3, 1);

    (*integralIS) += exp(lnew-l0);

    vtGradG(Gnew,thnew,n);

    lnew+= vtGradLogdet(Gnew, n);

    double p= exp(lnew - lold);

    double u = runif();

    if (u <= p)	{

      (*paccept) += 1;

      double ltemp= lnew;

      lnew = lold; lold = ltemp;

      double *thtemp;

      thtemp= thnew; thnew= thold; thold= thtemp;

      double *pitemp;

      pitemp= pinew; pinew= piold; piold= pitemp;

      double **Gtemp;

      Gtemp= Gnew; Gnew= Gold; Gold= Gtemp;

    } 



    if (i>=burnin) {

      int idx= i-burnin;

      for (int j=0; j<n; j++) pi[idx+j*(niter-burnin)]= piold[j];

    }

  }

  (*paccept) = (*paccept)/(niter+.0);

  (*integralIS) = l0 + log(*integralIS) - log(niter+.0);



  delete [] thmode;

  delete [] thold;

  delete [] thnew;

  delete [] piold;

  delete [] pinew;

  free_dmatrix(Gold,0,n,0,n);

  free_dmatrix(Gnew,0,n,0,n);

  free_dmatrix(cholS,1,n-1,1,n-1);

  free_dmatrix(cholSinv,1,n-1,1,n-1);

}





double Casper::calculateIntegral(int method) {

  int n = model->count();

  double* mode = calculateMode();

  double ans= calculateIntegral(mode, n, method);

  delete [] mode;

  return ans;

}



double Casper::calculateIntegral(double *mode, int n, int method) {

  double ans;

  if (method==1) {

    ans= LaplaceApprox(mode,n);

  } else if (method==2) {

    double *pi=NULL, paccept;

    IPMH(pi, &paccept, &ans, is_runs, is_runs, mode);  //no samples stored, simply reports average joint / proposal ratio

  }

  return ans;

}



double Casper::LaplaceApprox(double *mode, int n)

{

  if (n == 1) { return priorLikelihoodLn(mode); }

  double *thmode, ***H, **G, **S;



  thmode = new double[n - 1];

  mlogit(thmode, mode, n);



  H= darray3(n,n,n);

  vtHess(H, thmode, n);



  G = dmatrix(0,n,0,n);

  vtGradG(G,thmode, n);



  S= dmatrix(1,n-1,1,n-1);

  normapprox(S, G, H, mode, thmode, n, 1);



  double emlk = priorLikelihoodLn(mode);

  double gdet = vtGradLogdet(G, n);

  double **cholS;

  double sdet;

  cholS= dmatrix(1,n-1,1,n-1);

  bool posdef;
  
  choldc(S,n-1,cholS,&posdef);
  
  if (!posdef) {

    int i; double lmin=0, *vals;

    vals= dvector(1,n);

    eigenvals(S,n-1,vals);

    for (i=1; i<n; i++) if (vals[i]<lmin) lmin= vals[i];

    lmin = -lmin + .001;

    for (i=1; i<n; i++) S[i][i] += lmin;

    choldc(S,n-1,cholS,&posdef);

    free_dvector(vals,1,n);

  }

  sdet= choldc_det(cholS,n-1);

  free_dmatrix(cholS, 1, n-1, 1, n-1);


  double integral = emlk + gdet + (double)(n - 1) / 2.0 * log(2 * M_PI) - 0.5 * sdet;


  delete [] thmode;

  free_darray3(H,n,n,n);

  free_dmatrix(G,0,n,0,n);

  free_dmatrix(S,1,n-1,1,n-1);

	

  return integral;

}





bool Casper::isValid()

{

	list<Fragment*>::const_iterator fi;

	for (fi = frame->data.begin(); fi != frame->data.end(); fi++)

	{

		Fragment* f = *fi;

		if (mempprobs.count(f) == 0 || mempprobs[f].size() == 0)

        {

	        return false;

		}

	}



	return true;

}



double Casper::priorLn(double* pi)

{

	int n = model->count();



	double* alpha = new double[n];

	for (int i = 0; i < n; i++) alpha[i] = priorq;



	int log = 1;

	double ans= ddirichlet(pi, alpha, &n, &log);



	delete [] alpha;



	return ans;

}



double Casper::likelihoodLn(double* pi)

{

	double outer = 0;



	map<Fragment*, map<Variant*, double> >::const_iterator fi;

	for (fi = mempprobs.begin(); fi != mempprobs.end(); fi++)

	{

		double inner = 0;



		map<Variant*, double>::const_iterator vi;

		for (vi = fi->second.begin(); vi != fi->second.end(); vi++)

		{

			int i = model->indexOf(vi->first);

			inner += pi[i] * vi->second;

		}

		outer += fi->first->count * log(inner);

	}

	return outer;

}

double Casper::priorLikelihoodLn(double* pi)

{

	return priorLn(pi) + likelihoodLn(pi);

}



map<Fragment*, double> Casper::fragdist(double* pi)

{

	map<Fragment*, double> mem;



	map<Fragment*, map<Variant*, double> >::const_iterator fi;

	for (fi = mempprobs.begin(); fi != mempprobs.end(); fi++)

	{

		double sum = 0;



		map<Variant*, double>::const_iterator vi;

		for (vi = fi->second.begin(); vi != fi->second.end(); vi++)

		{

			int i = model->indexOf(vi->first);

			sum += pi[i] * vi->second;

		}

		mem[fi->first] = sum;

	}



	return mem;

}







void Casper::asymptoticSE(double *se, double *mode, int n) {

  int Sidx_ini=1;

  double **G, ***H, *thmode, **S, **Sinv;



  thmode = new double[n - 1];

  mlogit(thmode, mode, n);



  H= darray3(n,n,n);

  vtHess(H, thmode, n);



  G = dmatrix(0,n-1,0,n-2);

  vtGradG(G,thmode, n);


  
  S= dmatrix(1,n-1,1,n-1); Sinv= dmatrix(1,n-1,1,n-1);

  normapprox(Sinv, G, H, mode, thmode, n, Sidx_ini);

  bool posdef;

  inv_posdef(Sinv,n-1,S,&posdef);



  double **GtS= dmatrix(0,n-1,1,n-1);

  AB(G, 0, n-1, 0, n-2, S, 1, n-1, 1, n-1, GtS);



  int i,j;

  for (i=0; i<n; i++) {

    for (j=0, se[i]=0; j<n-1; j++) se[i] += GtS[i][j+1] * G[i][j];

    se[i] = sqrt(se[i]);

  }



  delete [] thmode;

  free_darray3(H,n,n,n); free_dmatrix(G,0,n-1,0,n-2);

  free_dmatrix(S,1,n-1,1,n-1); free_dmatrix(Sinv,1,n-1,1,n-1);

  free_dmatrix(GtS,0,n-1,1,n-1);

}





void Casper::normapprox(double **S, double *mode, int n, int Sidx_ini) {

  double **G, ***H, *thmode;



  thmode = new double[n - 1];

  mlogit(thmode, mode, n);



  H= darray3(n,n,n);

  vtHess(H, thmode, n);



  G = dmatrix(0,n,0,n);

  vtGradG(G,thmode, n);



  normapprox(S, G, H, mode, thmode, n, Sidx_ini);



  delete [] thmode;

  free_darray3(H,n,n,n);

  free_dmatrix(G,0,n,0,n);

}





void Casper::normapprox(double **S, double** G, double*** H, double* mode, double* thmode, int n, int Sidx_ini)

{

  map<Fragment*, double> mem = fragdist(mode);



  int rowS, colS;
  double mode_max;

  for (int l = 0; l < n - 1; l++) {

    rowS= l+Sidx_ini;

      for (int m = l; m < n - 1; m++) {

	  colS= m+Sidx_ini;

	  S[rowS][colS] = 0;

	  map<Fragment*, map<Variant*, double> >::const_iterator fi;

	  for (fi = mempprobs.begin(); fi != mempprobs.end(); fi++)

	    {

	      double term1 = 0, term2 = 0, term3 = 0;

	      map<Variant*, double>::const_iterator vi;

	      for (vi = fi->second.begin(); vi != fi->second.end(); vi++)

		{

		  int d = model->indexOf(vi->first);

		  double P = vi->second;

		  term1 += P * H[d][l][m];

		  term2 += P * G[d][l];

		  term3 += P * G[d][m];

		}

	      S[rowS][colS] -= fi->first->count * (term1 * mem[fi->first] - term2 * term3) / pow(mem[fi->first], 2);

	    }

	  for (int d = 0; d < n - 1; d++) {
	    mode_max=max_xy(mode[d], 1e-8);
	    S[rowS][colS] -= (priorq - 1.0) * (H[d][l][m] * mode_max - G[d][l] * G[d][m]) / pow(mode_max, 2);
          }

	  if (l != m) S[colS][rowS] = S[rowS][colS];

	}

    }



}





void Casper::mlogit(double *theta, double* pi, int n) {



  for (int i = 0; i < n - 1; i++) theta[i] = log(pi[i + 1] / pi[0]);



}

void Casper::milogit(double *pi, double* theta, int n)

{

	double sum = 1.0;



	for (int i = 0; i < n - 1; i++) sum += exp(theta[i]);



	pi[0] = 1.0 / sum;

	for (int i = 0; i < n - 1; i++) pi[i + 1] = exp(theta[i]) / sum;



}

void Casper::vtGradG(double **G, double* th, int n)

{

	double sum = 1.0;

	for (int i = 0; i < n - 1; i++)

	{

		sum += exp(th[i]);

	}



	for (int l = 0; l < n - 1; l++)

	{

		G[0][l] = -exp(th[l]) / pow(sum, 2);

	}



	for (int d = 1; d < n; d++)

	{

		for (int l = 0; l < n - 1; l++)

		{

			if (l != d - 1)

			{

				G[d][l] = -exp(th[d - 1] + th[l]) / pow(sum, 2);

			}

			else

			{

				G[d][l] = -exp(th[d - 1] + th[l]) / pow(sum, 2) + exp(th[l]) / sum;

			}

		}

	}

	//return G;

}

double Casper::vtGradLogdet(double** G, int n)

{

  bool posdef;

  double** GS = &G[1];

  double mydet = det(GS, n - 1, &posdef);

  if (mydet < 0) mydet = -mydet;

  double logdet = log(mydet);

  return logdet;

}

void Casper::vtHess(double ***H, double* th, int n)

{

	double sum = 1.0;

	for (int i = 0; i < n - 1; i++)

	{

		sum += exp(th[i]);

	}



	for (int d = 0; d < n; d++)

	{

		for (int l = 0; l < n - 1; l++)

		{

			for (int m = l; m < n - 1; m++)

			{

				if (d == 0)

				{

					if (m == l)

					{

						H[0][l][l] = -exp(th[l]) / pow(sum, 2) + 2.0 * exp(2.0 * th[l]) / pow(sum, 3);

					}

					else

					{

						H[0][l][m] = H[0][m][l] = 2.0 * exp(th[l] + th[m]) / pow(sum, 3);

					}

				}

				else

				{

					if (m == l)

					{

						if (l != d - 1)

						{

							H[d][l][l] = -exp(th[d - 1] + th[l]) / pow(sum, 2) + 2 * exp(th[d - 1] + 2 * th[l]) / pow(sum, 3);

						}

						else

						{

							H[d][l][l] = -2 * exp(2 * th[d - 1]) / pow(sum, 2) + 2 * exp(3 * th[d - 1]) / pow(sum, 3) + exp(th[d - 1]) / sum - exp(2 * th[d - 1]) / pow(sum, 2);

						}

					}

					else

					{

						if (m == d - 1)

						{

							H[d][l][m] = H[d][m][l] = -exp(th[d - 1] + th[l]) / pow(sum, 2) + 2 * exp(th[d - 1] + th[m] + th[l]) / pow(sum, 3);

						}

						else

						{

							if (l == d - 1)

							{

								H[d][l][m] = H[d][m][l] = 2 * exp(th[d - 1] + th[l] + th[m]) / pow(sum, 3) - exp(th[m] + th[l]) / pow(sum, 2);

							}

							else

							{

								H[d][l][m] = H[d][m][l] = 2 * exp(th[d - 1] + th[l] + th[m]) / pow(sum, 3);

							}

						}

					}

				}

			}

		}

	}



}



double Casper::det(double** a, int n, bool *posdef)

{

	double **aout = dmatrix(0, n - 1, 0, n - 1);

	int i,j,k;

	double sum;
	*posdef= true;


	for (i=0;i<n;i++) { for (j=i;j<n;j++) { aout[i][j]= a[i][j]; } }  //copy a into aout

	for (i=0;i<n;i++) {

		for (j=i;j<n;j++) {

			for (sum=aout[i][j],k=i-1;k>=0;k--) sum -= aout[i][k]*aout[j][k];

			if (i == j) {

			  if (sum <= 0.0) *posdef= false;
			  aout[i][i]=sqrt(sum);

			} else aout[j][i]=sum/aout[i][i];

		}

	}

	for (i=0;i<n;i++) { for (j=i+1;j<n;j++) { aout[i][j]= 0; } }  //set upper-diagonal elem to 0





	double det;

	for (det=1,i=0; i<n; i++) { det *= aout[i][i]*aout[i][i]; }

	

	free_dmatrix(aout, 0, n - 1, 0, n - 1);

	return(det);

}

int Casper::totCounts()

{

  int totC = this->frame->totCounts();

  return(totC);

}
