#ifndef __TOOLS__
#define __TOOLS__

#include <mex.h>

#if (defined(_WIN32) || defined(__WIN32__) )
#define DRAW mexEvalString("drawnow;");
#else
#define DRAW ;
#endif

#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)<(B)?(A):(B))

	
typedef struct{
	int sparseOrPSFOrMF; 
   /* Defines whether its a sparse, PSF, or matrix function,
     sparse=0, PSF=1, MF=2 */
	
	/* Definitions for sparse */
	double *Av; 
	mwIndex *Ar;
	mwIndex *Ac;

	/* Definitions for PSF */
	double *PSF;
	int PSFnr;
	int PSFnc;
	double *center;
	int xu,xv;

    /* Definitions for function handle for matrix free procedure */
    mxArray *mf;

	/* Definitions for all */
	int Anc;
	int Anr;

	
} Atype;

typedef struct{
	int dim;
	int m,n,l;
	int prodDims;
} Dtype;


typedef struct{
    int     dataitem;
    struct listelement *link;
} listelement;


double minf(double *x,int N);

double maxf(double *x,int n,int N);

double ddot_(int *n, double *x, int *incx, double *y, int *incy);

double dnrm2_(int *n, double *x, int *incx);

void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);

void dcopy_(int *p_n,double *p_x,int *incx,double *p_y,int *incy);

void P(double *y,int c,double *l,double *u,int N23);

double DTD(double *x,double *Nablafx, double *uijl, double tau, Dtype D);

void Amul(Atype A, double *x, double *y);

void ATmul(Atype A, double *x, double *y);

listelement *AddItem(listelement *listpointer, int data);

void ClearQueue(listelement *listpointer);

void PrintQueue(listelement * listpointer);

int QueueLength(listelement *listpointer);

void WriteQueueData(listelement *listpointer, double *rp, int l);

listelement *RemoveItem(listelement *listpointer);

#endif
