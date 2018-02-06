#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "tools.h"

/* Settings which makes the user do a CTRL-C break out of the loop */
#if defined(LIBUT) && (defined(_WIN32) || defined(__WIN32__) )

#define STOPMARK utIsInterruptPending()
#define INITBREAK ;
bool utIsInterruptPending(void);

#else

#include <signal.h>
#define INITBREAK   sigint_cont = 1;  (void) signal(SIGINT , ex_sigint);
#define STOPMARK sigint_cont==0
int sigint_cont = 1;
void ex_sigint(int sig) {
	sigint_cont = 0;
}
#endif

void tvreg_gpbb_core(double *xk,double *fxk,double *hxk,double *gxk,double *fxkl,double *kend,Atype A,double *b,double alpha,double tau,double epsb_rel,int k_max,Dtype D,int ctype,double *d,double *c, int ghxl, int xl,double *hxkl,double *gxkl,double *xlist,int K, double beta, double sigma,int verbose,double* numGrad,double* numBack,double* numFunc){

  double *xkp1,*xkm1,*Nablafxkp1,*Nablafxk,*Nablafxkm1,*dNablaf,*dx,*uijl,*tv,*tv2;

  double fxkp1,hxkp1,gxkp1,fr,nGt,t,den,B,betat,alphat,Nablafc;
  int i, j, k=0,prodDims,stop=0,one=1,counter,start;

  INITBREAK

	prodDims=D.prodDims;
  
	Nablafxk = malloc(prodDims*sizeof(double));
	Nablafxkm1 = malloc(prodDims*sizeof(double));
	Nablafxkp1 = malloc(prodDims*sizeof(double));

	dNablaf = malloc(prodDims*sizeof(double));

	xkp1 = malloc(prodDims*sizeof(double));
	xkm1 = malloc(prodDims*sizeof(double));
	dx = malloc(prodDims*sizeof(double));


	/*temp vectors */
	tv = malloc(prodDims*sizeof(double));
	tv2 = malloc(A.Anr*sizeof(double));
    uijl = malloc(D.dim*sizeof(double));
  

  /* INITIALIZE */
  /* Project solution onto feasible space */
	P(xk,ctype,d,c,prodDims);
	numGrad[0]=0;numBack[0]=0;numFunc[0]=0;


  /* LOOP */
	stop = 0; /*Flag for when to break the for-loop*/

  for(k=0; k<=k_max; k++){

		if(k==0){
			/* Calculate the gradient in xk */									 
			hxk[0] = alpha*DTD(xk,Nablafxk,uijl,tau,D);
		
			for(i=0;i<prodDims;i++){	Nablafxk[i]=alpha*Nablafxk[i]; }

			for(i=0;i<A.Anr;i++){ tv2[i]=0; }
		
			Amul(A,xk,tv2); /*MULTIPLY WITH A*********/
		
			for(i=0;i<A.Anr;i++){	tv2[i] = tv2[i]-b[i]; }
		
			gxk[0] = 0.5*pow( dnrm2_(&A.Anr,tv2,&one),2);
		
			fxk[0] = hxk[0]+gxk[0];
		
			ATmul(A,tv2,Nablafxk) ;/*MULTIPLY WITH AT ********/

			numFunc[0]+=1;numGrad[0]+=1;

		}


		if(k>0){
			for(i=0;i<prodDims;i++){dNablaf[i]=Nablafxk[i]-Nablafxkm1[i];}
		}


		/* append objective function*/
		fxkl[k] = fxk[0];

		if(ghxl){
			hxkl[k] = hxk[0];
			gxkl[k] = gxk[0];
		}
			
		/* store the iterate if requested */
		if(xl){
			for(i=0;i<prodDims;i++){
				xlist[k*prodDims+i] = xk[i];
			}
		}

		/* Backtracking parameters*/
		start = k-K;
		if(start < 0){
			start = 0;
		}
			
		fr = maxf(fxkl,start,k);
		betat = beta;
		counter = 0;

		/* Calculate the initial stepsize via the Barzilai-Borwein strategy*/
		if(k>0){
			den = ddot_(&prodDims,dNablaf,&one,dx,&one);
			if(den < 1e-25){
				den = 1e-25;
				printf("Small denuminator");
			}
			alphat = dnrm2_(&prodDims,dx,&one);
			alphat = alphat*alphat/den;
		}
		else{
			alphat = 1;
		}
			
		/* Take the projected step from xk to xkp1 */
		t = - alphat;
		dcopy_(&prodDims,xk,&one,xkp1,&one);
		daxpy_(&prodDims,&t,Nablafxk,&one,xkp1,&one);

		P(xkp1,ctype,d,c,prodDims);
		
		for(i=0;i<prodDims;i++){ tv[i] = xk[i]-xkp1[i]; }
		
		Nablafc = ddot_(&prodDims,Nablafxk,&one,tv,&one);		
		

		/* Calculate some things for xkp1*/
		hxkp1 = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,D);
		
		for(i=0;i<A.Anr;i++){ tv2[i]=0; }
		
		Amul(A,xkp1,tv2); /*MULTIPLY WITH A*********/
		
		for(i=0;i<A.Anr;i++){	tv2[i] = tv2[i]-b[i]; }
		
		gxkp1 = 0.5*pow( dnrm2_(&A.Anr,tv2,&one),2);
		
		fxkp1 = hxkp1+gxkp1;
		
		numFunc[0]+=1;

		while( (fxkp1 > fr -sigma*Nablafc) && counter < 14){
			numBack[0]+=1;
			counter++;
			betat = betat*betat;

			/* Take the projected step from xk to xkp1 */
			t = - alphat*betat;
			dcopy_(&prodDims,xk,&one,xkp1,&one);
			daxpy_(&prodDims,&t,Nablafxk,&one,xkp1,&one);

			P(xkp1,ctype,d,c,prodDims);
		
			for(i=0;i<prodDims;i++){ tv[i] = xk[i]-xkp1[i]; }
			
			Nablafc = ddot_(&prodDims,Nablafxk,&one,tv,&one);		
		
			/* Calculate some things for xkp1*/
			hxkp1 = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,D);
		
			for(i=0;i<A.Anr;i++){	tv2[i]=0; }
			
			Amul(A,xkp1,tv2); /*MULTIPLY WITH A*********/

			for(i=0;i<A.Anr;i++){ tv2[i] = tv2[i]-b[i]; }

			gxkp1 = 0.5*pow( dnrm2_(&A.Anr,tv2,&one),2);

			numFunc[0]+=1;

			fxkp1 = hxkp1 + gxkp1;		
		}
		
		/* do the remaining computations to obtain the gradient in xkp1*/
		
		for(i=0;i<prodDims;i++){	Nablafxkp1[i]=alpha*Nablafxkp1[i]; }
		ATmul(A,tv2,Nablafxkp1) ;/*MULTIPLY WITH AT ********/
		numGrad[0]+=1;
		

		/* obtain the delta x used to obtain the Barziliai-Borwein strategy */
		for(i=0;i<prodDims;i++){ dx[i]=xkp1[i]-xk[i];}




		/* Check the stopping criteria*/
		nGt = dnrm2_(&prodDims,dx,&one)/(-t);

		
		if(nGt<epsb_rel*prodDims){
			stop=1;
		}

		if(verbose)
			printf("k=%6d  f(x^k+1)=%e  ||G_1/t(x^k+1)||=%e\n",k,fxkp1,nGt);DRAW;

		if(stop || STOPMARK || k==k_max)
			goto cleanup;
				
		
		/* Update the values for the next iteration*/

		dcopy_(&prodDims,xkp1,&one,xk,&one);
		
		/* Save the old gradient*/
		dcopy_(&prodDims,Nablafxk,&one,Nablafxkm1,&one);
		dcopy_(&prodDims,Nablafxkp1,&one,Nablafxk,&one);
		
		hxk[0] = hxkp1;
		gxk[0] = gxkp1;
		fxk[0] = fxkp1;

	}

  cleanup:

  free(Nablafxk);
  free(Nablafxkm1);
  free(Nablafxkp1);
	free(dNablaf);
	free(xkp1);
	free(xkm1);
	free(dx);
	
	free(tv);
	free(tv2);
	free(uijl);

	kend[0]=(double)k;
}


void tvreg_upn_core(double *xkp1,double *fxkp1,double *hxkp1,double *gxkp1,double *fxkp1l,double *kend,Atype A,double *b,double alpha,double tau,double bL,double bmu,double epsb_rel,int k_max,Dtype D,int ctype,double *d,double *c, int ghxl, int xl,double *hxkp1l,double *gxkp1l,double *xlist,int qs,int verbose,double *numGrad,double* numBack,double *numFunc,double *numRest,double *Lklist,double *muklist,listelement *p_rp){

  double *yk,*xk,*Nablafyk,*Nablafxkp1,*uijl,*tv,*tv2;

  double fxk,fyk,q,thetak,thetakp1,betak,nGt,t,s_L=1.3,s_mu=0.7,Lm1,nGtm1,gamma0;
	double cumprod;
  int i, j, k=0,prodDims,stop=0,one=1,kk=-1;

  INITBREAK

	prodDims=D.prodDims;
  
	Nablafyk = malloc(prodDims*sizeof(double));
	Nablafxkp1 = malloc(prodDims*sizeof(double));
	yk = malloc(prodDims*sizeof(double));
	xk = malloc(prodDims*sizeof(double));

	/*temp vectors */
	tv = malloc(prodDims*sizeof(double));
	tv2 = malloc(A.Anr*sizeof(double));
    uijl = malloc(D.dim*sizeof(double));
  
	
    /* INITIALIZE */
	numGrad[0]=0;numBack[0]=0;numFunc[0]=0;numRest[0]=0;



 restart:
	cumprod=1.0;
  /* Project solution onto feasible space */

	P(xkp1,ctype,d,c,prodDims);
	
	if(qs==0)
		thetak=1;
	else
		thetak=sqrt(bmu/bL);

	/* Include a backtracking to initialize everything */
	/* Calculate the gradient in yk=xk+1 */									 
	numGrad[0]+=1;
	dcopy_(&prodDims,xkp1,&one,yk,&one);

	fyk = alpha*DTD(yk,Nablafyk,uijl,tau,D);
	
	for(i=0;i<prodDims;i++){	Nablafyk[i]=alpha*Nablafyk[i]; }
	
	for(i=0;i<A.Anr;i++){ tv2[i]=0; }
	
	Amul(A,yk,tv2); /*MULTIPLY WITH A*********/
	for(i=0;i<A.Anr;i++){	tv2[i] = tv2[i]-b[i]; }
	
	numFunc[0]+=1;
	fyk = fyk + 0.5*pow( dnrm2_(&A.Anr,tv2,&one),2);
	
	for(i=0;i<prodDims;i++){	tv[i]=0; }
	
	ATmul(A,tv2,Nablafyk) ;/*MULTIPLY WITH AT ********/

	/* Take the projected step from yk to xkp1 */
	t = - 1/bL;
	dcopy_(&prodDims,yk,&one,xkp1,&one);
	daxpy_(&prodDims,&t,Nablafyk,&one,xkp1,&one);
	
	P(xkp1,ctype,d,c,prodDims);
	
	/* Backtracking on Lipschitz parameter. */
	for(i=0;i<prodDims;i++){ tv[i] = xkp1[i]-yk[i]; }
	
	hxkp1[0] = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,D);
	
	for(i=0;i<A.Anr;i++){	tv2[i]=0; }
	
	Amul(A,xkp1,tv2); /*MULTIPLY WITH A*********/
	
	for(i=0;i<A.Anr;i++){ tv2[i] = tv2[i]-b[i]; }
	
	gxkp1[0] = 0.5*pow( dnrm2_(&A.Anr,tv2,&one),2);
	
	numFunc[0]+=1;
	fxkp1[0] = hxkp1[0] + gxkp1[0];		
	
	while( fxkp1[0]/(1+1e-14) > fyk + ddot_(&prodDims,Nablafyk,&one,tv,&one) + (bL/2)*pow(dnrm2_(&prodDims,tv,&one),2) ){
		numBack[0]+=1;
		bL = s_L*bL;
				
		/* Take the projected step from yk to xkp1 */
		t = - 1/bL;
		dcopy_(&prodDims,yk,&one,xkp1,&one);
		daxpy_(&prodDims,&t,Nablafyk,&one,xkp1,&one);
		P(xkp1,ctype,d,c,prodDims);
		
		/* Backtracking on Lipschitz parameter. */
		for(i=0;i<prodDims;i++){ tv[i] = xkp1[i]-yk[i]; }
			
		hxkp1[0] = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,D);
		
		for(i=0;i<A.Anr;i++){	tv2[i]=0; }
		
		Amul(A,xkp1,tv2); /*MULTIPLY WITH A*********/
		
		for(i=0;i<A.Anr;i++){ tv2[i] = tv2[i]-b[i]; }
			
		
		gxkp1[0] = 0.5*pow( dnrm2_(&A.Anr,tv2,&one),2);
		
		numFunc[0]+=1;
		fxkp1[0] = hxkp1[0] + gxkp1[0];
	}
	
	/* Calculate initial gradient map */
	if(ctype==1){
		nGt = dnrm2_(&prodDims,Nablafxkp1,&one);
	}
	else{
		t = - 1/bL;
		dcopy_(&prodDims,xkp1,&one,tv,&one);
		daxpy_(&prodDims,&t,Nablafxkp1,&one,tv,&one);
		P(tv,ctype,d,c,prodDims);
		
		for(i=0;i<prodDims;i++){	tv[i] = xkp1[i]-tv[i]; }
		
		nGt = bL*dnrm2_(&prodDims,tv,&one);
	}

	/* save the initial parameter */
	Lm1 = bL;
	nGtm1 = nGt;
		
	/*end initial tracking */
	dcopy_(&prodDims,xkp1,&one,xk,&one);
	dcopy_(&prodDims,xkp1,&one,yk,&one);
	
    /* LOOP */
	stop = 0; /*Flag for when to break the for-loop*/

	/*Calculate fxk */
	fxk = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,D);
		
	for(i=0;i<A.Anr;i++){	tv2[i]=0; }

	Amul(A,xkp1,tv2); /*MULTIPLY WITH A*********/

	for(i=0;i<A.Anr;i++){ tv2[i] = tv2[i]-b[i]; }

	numFunc[0]+=1;
	fxk = fxk + 0.5*pow( dnrm2_(&A.Anr,tv2,&one),2);


    for(k=0; k<=k_max; k++){
		kk+=1;
		Lklist[kk]=bL;
		muklist[kk]=bmu;
		/* Calculate the gradient in yk */
		numGrad[0]+=1;
		fyk = alpha*DTD(yk,Nablafyk,uijl,tau,D);

		for(i=0;i<prodDims;i++){	Nablafyk[i]=alpha*Nablafyk[i]; }

		for(i=0;i<A.Anr;i++){ tv2[i]=0; }

		Amul(A,yk,tv2); /*MULTIPLY WITH A*********/
		for(i=0;i<A.Anr;i++){	tv2[i] = tv2[i]-b[i]; }
		
		numFunc[0]+=1;
		fyk = fyk + 0.5*pow( dnrm2_(&A.Anr,tv2,&one),2);

		for(i=0;i<prodDims;i++){	tv[i]=0; }
		
		ATmul(A,tv2,Nablafyk) ;/*MULTIPLY WITH AT ********/

	
		/* Update estimate of the strong convexity parameter as minimum of current value and the computed value between xk and yk */
		if(k != 0){
			for(i=0;i<prodDims;i++){ tv[i] = xk[i]-yk[i]; }

			bmu = MAX( MIN(2*(fxk*(1+1e-14)-(fyk+ddot_(&prodDims,Nablafyk,&one,tv,&one)))/pow(dnrm2_(&prodDims,tv,&one),2),bmu),0);
		}

		/* Take the projected step from yk to xkp1 */
		t = - 1/bL;
		dcopy_(&prodDims,yk,&one,xkp1,&one);
		daxpy_(&prodDims,&t,Nablafyk,&one,xkp1,&one);

		P(xkp1,ctype,d,c,prodDims);
		
		/* Backtracking on Lipschitz parameter. */
		for(i=0;i<prodDims;i++){ tv[i] = xkp1[i]-yk[i]; }
		
		hxkp1[0] = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,D);
		
		for(i=0;i<A.Anr;i++){	tv2[i]=0; }

		Amul(A,xkp1,tv2); /*MULTIPLY WITH A*********/

		for(i=0;i<A.Anr;i++){ tv2[i] = tv2[i]-b[i]; }

		gxkp1[0] = 0.5*pow( dnrm2_(&A.Anr,tv2,&one),2);
		
		numFunc[0]+=1;
		fxkp1[0] = hxkp1[0] + gxkp1[0];		

		while( fxkp1[0]/(1+1e-14) > fyk + ddot_(&prodDims,Nablafyk,&one,tv,&one) + (bL/2)*pow(dnrm2_(&prodDims,tv,&one),2) ){
			numBack[0]+=1;
			bL = s_L*bL;
				
			/* Take the projected step from yk to xkp1 */
			t = - 1/bL;
			dcopy_(&prodDims,yk,&one,xkp1,&one);
			daxpy_(&prodDims,&t,Nablafyk,&one,xkp1,&one);
			P(xkp1,ctype,d,c,prodDims);
			
			/* Backtracking on Lipschitz parameter. */
			for(i=0;i<prodDims;i++){ tv[i] = xkp1[i]-yk[i]; }
			
			
			hxkp1[0] = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,D);
			
			for(i=0;i<A.Anr;i++){	tv2[i]=0; }
			
			Amul(A,xkp1,tv2); /*MULTIPLY WITH A*********/
				
			for(i=0;i<A.Anr;i++){ tv2[i] = tv2[i]-b[i]; }
			
			
			gxkp1[0] = 0.5*pow( dnrm2_(&A.Anr,tv2,&one),2);

			numFunc[0]+=1;
			fxkp1[0] = hxkp1[0] + gxkp1[0];
		}
		
		fxkp1l[kk] = fxkp1[0];
		
		if(ghxl){
			hxkp1l[kk] = hxkp1[0];
			gxkp1l[kk] = gxkp1[0];
		}
			
		/* store the iterate if requested */
		if(xl){
			for(i=0;i<prodDims;i++){
				xlist[kk*prodDims+i] = xkp1[i];
			}
		}


		if(verbose)
			printf("k=%6d  f(x^k+1)=%e  ||G_L(x^k+1)||=%e  L_k=%.2e  mu_k=%.2e\n",kk,fxkp1[0],nGt,bL,bmu);DRAW;

		/* calculate the gradient in xkp1 */
		numGrad[0]+=1;
		for(i=0;i<prodDims;i++){	Nablafxkp1[i]=alpha*Nablafxkp1[i]; }
		ATmul(A,tv2,Nablafxkp1) ;/*MULTIPLY WITH AT ********/


		/* Check stopping criteria xkp1*/
		if(ctype==1){
			nGt = dnrm2_(&prodDims,Nablafxkp1,&one);
			if(nGt <= epsb_rel*prodDims){
				stop = 1;
   			/*overwrite xkp1 to return with*/
			  t = - 1/bL;
				daxpy_(&prodDims,&t,Nablafxkp1,&one,xkp1,&one); 
			}
		}
		else{
			t = - 1/bL;
			dcopy_(&prodDims,xkp1,&one,tv,&one);
			daxpy_(&prodDims,&t,Nablafxkp1,&one,tv,&one);
			P(tv,ctype,d,c,prodDims);

			for(i=0;i<prodDims;i++){	tv[i] = xkp1[i]-tv[i]; }
			
			nGt = bL*dnrm2_(&prodDims,tv,&one);
			if(nGt <= epsb_rel*prodDims){
				stop = 1;
   			/*overwrite xkp1 to return with*/
				daxpy_(&prodDims,&t,Nablafxkp1,&one,xkp1,&one);
				P(xkp1,ctype,d,c,prodDims);
			}			
		}

		if(stop || STOPMARK || kk==k_max){
			goto cleanup;
		}

		/* Check stopping criteria yk*/
		if(ctype==1){
			if(dnrm2_(&prodDims,Nablafyk,&one) <= epsb_rel*prodDims)
				stop = 1;
		}
		else{
			t = - 1/bL;
			for(i=0;i<prodDims;i++){	tv[i] = yk[i]-xkp1[i]; }
			if(bL*dnrm2_(&prodDims,tv,&one) <= epsb_rel*prodDims)
				stop = 1;
		}

		if(stop || STOPMARK || kk==k_max){
			goto cleanup;
		}

		/* Compute values for accelerated step size*/
		if(qs == 0 || qs == 2)
			q = 0;
		else{
			q = bmu/bL;
		
			if(k!=0)
				cumprod = cumprod*(1-sqrt(q));
			/*tjeck if the convergence rate is fast enough*/
			if(bmu >0){
				if(nGt*nGt> cumprod*(4*bL/bmu-bL/Lm1+4*gamma0*bL/pow(bmu,2))*nGtm1*nGtm1){
				  /*printf("not fast enough %d\n",kk);*/
					/*printf("%f %f %f %f\n",nGt,nGtm1*nGtm1,cumprod,(4*bL/bmu-bL/Lm1+4*gamma0*bL/pow(bmu,2)));*/

					/*list of restart positions */
					p_rp = AddItem(p_rp,kk);

					bmu=bmu*s_mu;
					numRest[0]+=1;
					goto restart;
				}
			}
		}
		/*printf("%f\n",q);DRAW;*/

		thetakp1 = (-(pow(thetak,2)-q)+sqrt(pow( pow(thetak,2)-q,2)+4*pow(thetak,2)))/2.0;
		betak = (thetak*(1-thetak))/(pow(thetak,2)+thetakp1);

		if(k==0){
			gamma0 = thetakp1*(thetakp1*bL-bmu)/(1-thetakp1);
			/*printf("gamma0 %f\n",gamma0);*/
		}

		if(qs==2)
			dcopy_(&prodDims,xkp1,&one,yk,&one);
		else{
			/* accelerated term*/
			/* yk = xkp1 + betak*(xkp1-xk) */
			for(i=0;i<prodDims;i++){ tv[i] = xkp1[i]-xk[i]; }

			dcopy_(&prodDims,xkp1,&one,yk,&one);
			daxpy_(&prodDims,&betak,tv,&one,yk,&one);

			/* Update the values for the next iteration*/
			thetak = thetakp1;
			dcopy_(&prodDims,xkp1,&one,xk,&one);
		}

		fxk = fxkp1[0];

	}

  cleanup:
  free(Nablafyk);
  free(Nablafxkp1);

  free(yk);
  free(xk);

  free(tv);
  free(tv2);
  free(uijl);

  kend[0]=(double)kk;

}
