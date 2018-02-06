#ifndef __TV_CORE__
#define __TV_CORE__


void tvi_gpbb_core(double *xk,double *fxk,double *hxk,double *gxk,double *fxkl,double *k,Atype A,double *b,double alpha,double tau,double epsb_rel,int k_max,Dtype D,int ctype,double *d,double *c, int ghxl, int xl,double *hxkl,double *gxkl,double *xlist,int K,double beta,double sigma,int verbose,double* numGrad,double* numBack,double* numFunc);

void tvi_upn_core(double *xkp1,double *fxkp1,double *hxkp1,double *gxkp1,double *fxkp1l,double *k,Atype A,double *b,double alpha,double tau,double bL,double bmu,double epsb_rel,int k_max,Dtype D,int ctype,double *d,double *c, int ghxl, int xl,double *hxkp1l,double *gxkp1l,double *xlist,int qs,int verbose,double *numGrad,double* numBack,double *numFunc,double *NumRest,double *Lklist,double *muklist,listelement *p_rp);


#endif
