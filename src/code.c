/* History: Sep 08 2017 Initial coding
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void myrpexp(int *pn, double *rate, double *tstar, double *ret);
void myrpexp0(int *pn, double *rate, double *pt, double *ret);
void getef(int *pn, double *tmp, double *pA, double *ret);
void logrnk(int *nn, int *nngroup, int *nstrat,\
	       double *rho, double *time, int *status,\
	       int *group, int *strata, double *obs,\
	       double *exp, double *var, double *risk,\
	       double *kaplan);


static const R_CMethodDef callMethods[] = {
  {"myrpexp", (DL_FUNC)&myrpexp, 4},
  {"myrpexp0", (DL_FUNC)&myrpexp0, 4},
  {"getef", (DL_FUNC)&getef, 4},
  {"logrnk", (DL_FUNC)&logrnk, 13},
  {NULL, NULL, 0}
};

/* Function that R will call */
void myrpexp(int *pn, double *rate, double *tstar, double *ret) {
  /* 
  int *pn;         Number of random numbers to return 
  double *rate;    Vector of length 2 for rate 
  double *tstar;   t.star vector (not including the first 0) 
  double *ret;     Vector of length n 
  */

  double scale1, scale2, *ptstar, re, *pret, ts;
  int i, n;

  scale1  = 1.0/rate[0];
  scale2  = 1.0/rate[1];
  n      = *pn;

  /* For random number generation */
  GetRNGstate();

  /* Note that the c function rexp uses the scale and not rate as the parm */
  for (i=0, pret=ret, ptstar=tstar; i<n; i++, pret++, ptstar++) {
    re = rexp(scale1);
    ts = *ptstar;
    if (re < ts) {
      *pret = re;
    } else {
      *pret = rexp(scale2) + ts;
    }
  }  

  /* For random number generation */
  PutRNGstate();

  return;

} /* END: myrpexp */
 
/* Function for calling rpexp with n > 1, length(rate)=length(t)=2 */
void myrpexp0(int *pn, double *rate, double *pt, double *ret) {
  /*
  int *pn;         Number of random numbers to return 
  double *rate;    Vector of length 2 for rate 
  double *pt;      Time (not the first 0) 
  double *ret;     Vector of length n 
  */

  double scale1, scale2, *pret, time2;
  int i, n;

  scale1  = 1.0/rate[0];
  scale2  = 1.0/rate[1];
  n      = *pn;
  time2  = *pt;

  /* For random number generation */
  GetRNGstate();

  /* Note that the c function rexp uses the scale and not rate as the parm */
  for (i=0, pret=ret; i<n; i++, pret++) *pret = rexp(scale1);
  for (i=0, pret=ret; i<n; i++, pret++) {
    if (*pret >= time2) *pret = rexp(scale2) + time2; 
  }  

  /* For random number generation */
  PutRNGstate();

  return;

} /* END: myrpexp0 */


void getef(int *pn, double *tmp, double *pA, double *ret) {
  /*
  int *pn;      Nupper 
  double *tmp;  Vector of log(runif(Nupper))/a in R code 
  double *pA;   A 
  double *ret;  Vector of length n initialized to 0 
  */

  int k, n;
  double t, *pt, *pret, A;

  n = *pn;
  A = *pA;
  t = 0.0;
  for (k=0, pt=tmp, pret=ret; k<n; k++, pt++, pret++) {
    t = t - *pt;
    if (t>A) break;
    *pret = t;
  }
  
  return;

} /* END: getef */

void logrnk(int *nn, int *nngroup, int *nstrat, double *rho, double *time, int *status, int *group, int *strata, double *obs,\
            double *exp, double *var, double *risk, double *kaplan) {
  /*
  int *nn, *nngroup, *nstrat, *status, *group, *strata;
  double *rho, *time, *obs, *exp, *var, *risk, *kaplan;
  */

    register int i,j,k;
    int kk;
    int n, ngroup, ntot;
    int istart, koff;
    double km, nrisk, wt, tmp;
    double deaths;

    ntot = *nn;
    ngroup = *nngroup;
    istart=0; koff=0;
    for (i=0; i< ngroup*ngroup; i++)  var[i]=0;
    for (i=0; i< *nstrat*ngroup; i++) {
	obs[i]=0;
	exp[i]=0;
	}

    while (istart < ntot) {  /* loop over the strata */
	for (i=0; i<ngroup; i++) risk[i]=0;

	/* last obs of this strata */
	for (i=istart; i<ntot; i++)
	    if (strata[i]==1) break;
	n = i+1;


	/*
	** Compute the k-m, which is only needed if rho!=0
	**   We want it set up as a left-continuous function (unusual)
	*/
	if (*rho !=0){
	    km =1;
	    for (i=istart; i<n; ) {
		kaplan[i] = km;
		nrisk = n-i;
		deaths =status[i];
		for (j=i+1; j<n && time[j]==time[i]; j++) {
		    kaplan[j] = km;
		    deaths += status[j];
		    }
		km = km * (nrisk-deaths)/nrisk;
		i=j;
		}
	    }

	/*
	** Now for the actual test
	*/
	for (i=n-1; i>=istart; i--) {
	    if (*rho ==0) wt=1;
	    else          wt= pow(kaplan[i], *rho);

	    deaths = 0;
	    for (j=i; j>=istart && time[j]==time[i]; j--) {
		k = group[j]-1;
		deaths += status[j];
		risk[k] += 1;
		obs[k + koff] += status[j] *wt;
		}
	    i= j +1;
	    nrisk = n-i;

	    if (deaths>0) {  /* a death time */
		for (k=0; k<ngroup; k++)
		    exp[k+koff] += wt* deaths * risk[k] / nrisk;

		if (nrisk==1) continue;  /*only 1 subject, so no variance */
		kk =0;
		wt = wt*wt;
		for (j=0; j<ngroup; j++) {
		    tmp = wt* deaths* risk[j]* (nrisk-deaths)/(nrisk *(nrisk-1));
		    var[kk+j] += tmp;
		    for (k=0; k<ngroup; k++) {
			var[kk] -= tmp * risk[k] / nrisk;
			kk++ ;
			}
		    }
		}
	    }
	istart = n;
	koff += ngroup;
	}

  return;

} /* END: logrnk */

void R_init_CGEN(DllInfo *dll)
{
    R_registerRoutines(dll, callMethods, NULL, NULL, NULL);
}

