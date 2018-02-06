#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include "tools.h"

/* Used at trp as the relatice accuracy needed*/
#define NTTOL_REL 1e-4
 
#ifdef BLAS

extern double ddot_(int *n, double *x, int *incx, double *y, int *incy);
extern double dnrm2_(int *n, double *x, int *incx);
extern void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
extern void dcopy_(int *p_n,double *p_x,int *incx,double *p_y,int *incy);

#else

/* Only implements the simplest version of ddot which is used in this code 
   here, with incrx and incry = 1 hardcoded */
double ddot_(int *p_n, double *p_x, int *incx, double *p_y, int *incy)
{
     register int i, n;
     register double *x, *y;
     register double ddot = 0;
     x = p_x;
     y = p_y;
     n = *p_n;
     for (i = 0; i < n; ++i)
          ddot += (*x++)*(*y++);

     return ddot;
}

double dnrm2_(int *p_n, double *p_x, int *incx)
{
	register int i,n;
	register double *x;
	register double nrm2 = 0;
	

	x = p_x;
	n = *p_n;

	for (i = 0; i < n; ++i){
			nrm2 += (*x)*(*x);
			x++; 
	}

	return sqrt(nrm2);
}

/* Only implements the simplest version of daxpy which is used in this code 
   here, with incrx and incry = 1 hardcoded. *y points at the results (overwrite) */
void daxpy_(int *p_n, double *p_alpha, double *p_x, int *incx, double *p_y, int *incy)
{
     register int i, n;
     register double *x, *y,alpha;
     x = p_x;
     y = p_y;
     n = *p_n;
     alpha = *p_alpha;

     for (i = 0; i < n; ++i){
          *y = alpha*(*x) + (*y);
					y++;x++;
    }
}

void dcopy_(int *p_n,double *p_x,int *incx,double *p_y,int *incy){
	register int i, n;
	register double *x, *y;
	x = p_x;
	y = p_y;
	n = *p_n;

	for (i = 0; i < n; ++i){
		*y = *x;
		y++;x++;
 	}
}
#endif



double minf(double *x,int N){
	double c;
	int i;

	c = x[0];

	for(i=1;i<N; i++)
		if(x[i]<c)
			c = x[i];
		
	return c;
}

double maxf(double *x,int n,int N){
	double c;
	int i;

	c = x[n];

	for(i=n+1;i<=N; i++)
		if(x[i]>c)
			c = x[i];
		
	return c;
}

/* Functions used for the inverse problems using Nesterov or BB method 
  Project onto the feasible (convex) set. 
   c==1: Unconstrained. 
   c==2: Lower and upper bounds (elementwise) on x. Inplace. */
void P(double *y,int ctype,double *d,double *c,int mnl){
	if(ctype==2){ /* c <= x <= d (elementwise) */
    register int i=0;
    for(i=0;i<mnl;i++){
    	if(y[i]<c[i])
				y[i]=c[i];
			else if(y[i]>d[i])
				y[i]=d[i];
		}
  }
}

/* Function used to calculate operations involving D and D^T*/
double DTD(double *x,double *Nablafx, double *uijl, double tau, Dtype D){
	
	register int i,u,v,w,i1,i2,i3,i4,mn,mnl,s1,s2,s3,s4;
	double tv_tau_x=0,c1,c2,taud2=tau/2,tau2=1/(tau*2);

	/* Clear the current gradient */
	for (i=0; i<D.prodDims; i++){ Nablafx[i] = 0.0; }

	if(D.dim ==2){
		for (u=0; u<=D.m-1; u++) {
			for (v=0; v<=D.n-1; v++) {
				i1 = (u+1)%D.m + v*D.m;
				i2 = u + ((v+1)%D.n)*D.m;
				i3= u+v*D.m;

				uijl[0] = x[i1]-x[i3];
				uijl[1] = x[i2]-x[i3];
					
				c1 = sqrt(uijl[0]*uijl[0] + uijl[1]*uijl[1]);
					
				if(tau<c1){
					c2 = c1;
					tv_tau_x += c1-taud2;
				}
				else{
					c2 =tau;
					tv_tau_x += c1*c1*tau2;
				}
					
				uijl[0] = uijl[0]/c2;
				uijl[1] = uijl[1]/c2;
					
				Nablafx[i1] += uijl[0];
				Nablafx[i3] -= uijl[0];
				Nablafx[i2] += uijl[1];
				Nablafx[i3] -= uijl[1];
			}
		}

	}

	else if(D.dim ==3){
		mn = D.m*D.n;
    mnl = mn*D.l;

		for (u=0; u<=D.m-1; u++) {
			for (v=0; v<=D.n-1; v++) {
				/*s1= ((u+1)%D.m) + v*D.m;*/
				/*s2 = u + ((v+1)%D.n)*D.m;*/
				s1= ((u+1)%D.m) + v*D.m;
				s2 = u + ((v+1)%D.n)*D.m;
				s3 = u+v*D.m;
				s4 = 0;
				for (w=0; w<=D.l-1; w ++){ 
					i1 = s1 + s4;
					i2 = s2 + s4;
					i4= s3+s4;
					i3= (i4+mn)%mnl;

					s4+=mn;
					
					uijl[0] = x[i1]-x[i4];
					uijl[1] = x[i2]-x[i4];
					uijl[2] = x[i3]-x[i4];
					
					c1 = sqrt(uijl[0]*uijl[0] + uijl[1]*uijl[1] + uijl[2]*uijl[2]);
					
					if(tau>c1){
						c2 =tau;
						tv_tau_x += c1*c1*tau2;
							}
					else{
						c2 = c1;
						tv_tau_x += c1-taud2;
					}
					
					uijl[0] = uijl[0]/c2;
					uijl[1] = uijl[1]/c2;
					uijl[2] = uijl[2]/c2;
					
					Nablafx[i1] += uijl[0];
					Nablafx[i4] -= uijl[0]+uijl[1]+uijl[2];
					Nablafx[i2] += uijl[1];
					Nablafx[i3] += uijl[2];
					
				}
			}
		}
	}
	else{printf("Incorrect dim variable, only dim=2 or dim=3 supported.\n");}

	return tv_tau_x;
}


/* Implements multiplication with A and AT 

 If A is sparse:
 using the definition of sparse structures as in Matlab. The Atype is a struct holding the values Av, and the corresponding rows Ar. The array Ac hold 
 the number of non-zero elements as the difference between two 
 indexes in Ac 

 If A is a PSF:
 implements a two dimentional PSF with reflective boundary conditions

 If A is a function:
 call back    
*/

void Amul(Atype A, double *x, double *y){
	register int j,i,jc,ic,nnzc,ki,kj,ii,jj,m,n,mm,nn;
    int one=1;
    double done=1.0;
	double sum;
    mxArray *lhs;
    mxArray *rhs[3];


	if(A.sparseOrPSFOrMF==0){ /* defines a sparse matrix*/
		j = 0;
		i = 0;
		for(jc=1;jc<=A.Anc;jc++){
			nnzc = A.Ac[jc]-A.Ac[jc-1];
			for(ic=0;ic<nnzc;ic++){
				y[A.Ar[i]] += A.Av[i]*x[j];
				i++;
			}
			j++;
		}
	}
	else if(A.sparseOrPSFOrMF==1){ /* defines a PSF*/
		for(i=0;i<A.xu;i++){
			for(j=0;j<A.xv;j++){
				sum=0;
				m = -1;
				for(mm=A.PSFnr-1;mm>=0;mm--){
					m++;

					/*Reflective boundary conditions*/
					ii = i +(m-A.center[0]+1);
					if(ii<0){ii=-ii;}
					if(ii>=A.xu){ii=A.xu-(ii-A.xu)-1;}
						
					n=-1;
					for(nn=A.PSFnc-1;nn>=0;nn--){
						n++;

						jj = j +(n-A.center[1]+1);

						/*Reflective boundary conditions*/
						if(jj<0){jj=-jj;}
						if(jj>=A.xv){jj=A.xv-(jj-A.xv)-1;}

						sum +=x[ii+jj*A.xu]*A.PSF[mm+nn*A.PSFnr];


						/* zero-boundary condition*/
						/*if(ii>=0 && jj>=0 && ii<A.xu && jj<A.xv)
							sum +=x[ii+jj*A.xu]*A.PSF[mm + nn*A.PSFnr];*/

					}
					
				}
				y[i+j*A.xu] = sum;
			}
		}
	}
    else{ /* matrix free*/
      
      rhs[0] = A.mf;
      rhs[1] = mxCreateDoubleMatrix(A.Anc, 1, mxREAL);
      memcpy(mxGetPr(rhs[1]), x, A.Anc*sizeof(double));
      rhs[2] = mxCreateDoubleScalar(1);

      mexCallMATLAB(1, &lhs, 3, rhs, "feval");      

      if( mxGetM(lhs) != A.Anr && mxGetN(lhs) != 1 )
        printf("mex:Expected an %d x %d but got a %d x %d\n",
                     A.Anr, 1, mxGetM(lhs), mxGetN(lhs));
      else
        daxpy_(&A.Anr, &done, mxGetPr(lhs), &one, y, &one);

      
      /* To make it modular release all mx arrays. */      
      mxDestroyArray(lhs);
      mxDestroyArray(rhs[1]);
    }
      
}

void ATmul(Atype A, double *x,double *y){
	register int j, i, jc, ic, nnzc, ki, kj, ii, jj, m, n, mm, nn;
    int one=1;
    double done=1.0;
	double sum;
    mxArray *lhs;
    mxArray *rhs[3];

	if(A.sparseOrPSFOrMF==0){
		j = 0;
		i = 0;
		for(jc=1;jc<=A.Anc;jc++){
			nnzc = A.Ac[jc]-A.Ac[jc-1];
			for(ic=0;ic<nnzc;ic++){
				y[j] += A.Av[i]*x[A.Ar[i]];
				i++;
			}
			
			j++;
		}
	}
	else if(A.sparseOrPSFOrMF==1){ /* defines a PSF*/
		for(i=0;i<A.xu;i++){
			for(j=0;j<A.xv;j++){
				m = -1;
				for(mm=A.PSFnr-1;mm>=0;mm--){
					m++;

					/*Reflective boundary conditions*/
					ii = i +(m-A.center[0]+1);
					if(ii<0){ii=-ii;}
					if(ii>=A.xu){ii=A.xu-(ii-A.xu)-1;}
						
					n=-1;
					for(nn=A.PSFnc-1;nn>=0;nn--){
						n++;


						jj = j +(n-A.center[1]+1);

						/*Reflective boundary conditions*/
						if(jj<0){jj=-jj;}
						if(jj>=A.xv){jj=A.xv-(jj-A.xv)-1;}

						y[ii+jj*A.xu] +=x[i+j*A.xu]*A.PSF[mm+nn*A.PSFnr];

					}
					
				}
			}
		}
	}
    else{ /* matrix free*/
      
      rhs[0] = A.mf;
      rhs[1] = mxCreateDoubleMatrix(A.Anr, 1, mxREAL);
      memcpy(mxGetPr(rhs[1]), x, A.Anr*sizeof(double));
      rhs[2] = mxCreateDoubleScalar(2);

      mexCallMATLAB(1, &lhs, 3, rhs, "feval");

      if( mxGetM(lhs) != A.Anc && mxGetN(lhs) != 1 )
        printf("mex:Expected an %d x %d but got a %d x %d\n",
              A.Anc, 1, mxGetM(lhs), mxGetN(lhs));
      else
        daxpy_(&A.Anc, &done, mxGetPr(lhs), &one, y, &one);


      /* To make it modular: relese the mx array. */
      mxDestroyArray(lhs);
      mxDestroyArray(rhs[1]);
    }


}

/*Function  for linked lists */


listelement *AddItem(listelement* listpointer, int data){

    listelement * lp = listpointer;

    if (listpointer != NULL) {
			while (listpointer -> link != NULL)
				listpointer = (listelement *) listpointer -> link;

			listpointer -> link = (struct listelement  *) malloc (sizeof (listelement));
			listpointer = (listelement *) listpointer -> link;
			listpointer -> link = NULL;
			listpointer -> dataitem = data;
			return lp;
    }
    else {
			listpointer = (listelement  *) malloc (sizeof (listelement));
			listpointer -> link = NULL;
			listpointer -> dataitem = data;
			return listpointer;
    }
}

listelement *RemoveItem(listelement* listpointer) {

    listelement * tempp;
    tempp = (listelement *) listpointer -> link;
    free (listpointer);
    return tempp;
}


void ClearQueue(listelement* listpointer){

	while (listpointer != NULL) {
		listpointer = RemoveItem(listpointer);
	}
}

void PrintQueue (listelement* listpointer) {

    if (listpointer == NULL)
			printf ("queue is empty!\n");
    else
			while (listpointer != NULL) {
				printf ("%d\t", listpointer -> dataitem);
				listpointer = (listelement *) listpointer -> link;
			}
    printf ("\n");
}


int QueueLength(listelement* listpointer){

	int l=0;

	if (listpointer == NULL)
		return 0;
	else{
		while (listpointer != NULL) {
			l++;
			listpointer = (listelement *) listpointer -> link;
		}
			return l;

	}

   /* if (listpointer == NULL)
			printf ("queue is empty!\n");
    else
			while (listpointer != NULL) {
				printf ("%d\t", listpointer -> dataitem);
				listpointer = (listelement *) listpointer -> link;
			}
    printf ("\n");

		return 0;*/
}

void WriteQueueData(listelement* listpointer, double* rp, int l){

	int i=-1;

	while (listpointer != NULL) {
		if(i != -1){ /*dont write the first element, since its a dummy */
			rp[i] = (double) listpointer->dataitem;
		}
		listpointer = (listelement *) listpointer->link;
		i++;
	}
}
