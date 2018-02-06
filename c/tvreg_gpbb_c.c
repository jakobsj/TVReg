#include <mex.h>
#include "tools.h"
#include "tv_core.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	register double *b,*d,*c,*x,*dims,alpha,tau,gamma_rel,*hxkl,*gxkl,*xlist;
	register double *xk,*fxk,*hxk,*gxk,*fxkl,*k,beta,sigma,*numGrad,*numBack,*numFunc;
    mxArray *dims_lhs;
    mxArray *dims_rhs[3];
    double *dimsp;

	mxArray *M,*S,*Mdims;
	int i,j,k_max,N,dim,prodDims,one=1,ctype,ghxl,xl,K,verbose;  
	Atype A;
	Dtype D;

	if(nrhs != 17){
      mexErrMsgTxt("Should contain 17 input parameters\n");}
    else{		
		if(mxIsStruct(prhs[0])){
			A.sparseOrPSFOrMF = 1;
			M = mxGetField(prhs[0], 0, "center");
			A.center = mxGetPr(M);

			M = mxGetField(prhs[0], 0, "PSF");
			A.PSF = mxGetPr(M);
			A.PSFnr = mxGetM(M);
			A.PSFnc = mxGetN(M);
		}
        else if(mxIsClass(prhs[0] , "function_handle")){
          A.sparseOrPSFOrMF = 2;
          dims_rhs[0] = (mxArray*)prhs[0];
          dims_rhs[1] = mxCreateDoubleScalar(0.0);
          dims_rhs[2] = mxCreateDoubleScalar(0.0);
          
          mexCallMATLAB(1, &dims_lhs, 3, dims_rhs, "feval");
          dimsp = mxGetPr(dims_lhs);
          
          /* save for later */
          A.mf = (mxArray*) prhs[0];
          A.Anr = (int)dimsp[0];
          A.Anc = (int)dimsp[1];
        }
		else{
			A.sparseOrPSFOrMF = 0;
			M = (mxArray*)prhs[0]; /* Pointer to matrix structure*/
			A.Av = mxGetPr(M); /* set the nessecary values in a struct*/
			A.Ar = mxGetIr(M);
			A.Ac = mxGetJc(M);
			A.Anc = mxGetN(M);
			A.Anr = mxGetM(M);
		}


		M = (mxArray*)prhs[1];
		b = mxGetPr(M);

		if(mxIsStruct(prhs[0])){
			A.Anr = mxGetM(M)*mxGetN(M);
		}
		
		S = (mxArray*)prhs[2];
		alpha = (double)(mxGetScalar(S));
		
		S = (mxArray*)prhs[3];
		tau = (double)(mxGetScalar(S));
		
		Mdims = (mxArray*)prhs[4];
		dims = mxGetPr(Mdims);
			
		S = (mxArray*)prhs[5];
		gamma_rel = (double)(mxGetScalar(S));

		S = (mxArray*)prhs[6];
		k_max = (int)(mxGetScalar(S));

		M = (mxArray*)prhs[7];
		x = mxGetPr(M);

		if(mxIsStruct(prhs[0])){
			A.Anc = mxGetM(M)*mxGetN(M);
		}

		S = (mxArray*)prhs[8];
		ctype = (int)(mxGetScalar(S));

		M = (mxArray*)prhs[9];
		d = mxGetPr(M);

		M = (mxArray*)prhs[10];
		c = mxGetPr(M);

		S = (mxArray*)prhs[11];
		ghxl = (int)(mxGetScalar(S));

		S = (mxArray*)prhs[12];
		xl = (int)(mxGetScalar(S));

		S = (mxArray*)prhs[13];
		K = (int)(mxGetScalar(S));

		S = (mxArray*)prhs[14];
		beta = (double)(mxGetScalar(S));

		S = (mxArray*)prhs[15];
		sigma = (double)(mxGetScalar(S));

		S = (mxArray*)prhs[16];
		verbose = (double)(mxGetScalar(S));

		/*obtain the dimensions */

		dim = MAX( mxGetM(Mdims), mxGetN(Mdims) );
		prodDims=1;
		for(i=0;i<dim;i++){
			prodDims = prodDims*(int)dims[i];}

		D.dim=dim;
		D.m=(int)dims[0];D.n=(int)dims[1];
		if(dim==3){
			D.l=(int)dims[2];
		}
		D.prodDims =prodDims;

		if(mxIsStruct(prhs[0])){
			A.xu = (int)dims[0];
			A.xv = (int)dims[1];
		}

		if( mxIsStruct(prhs[0]) && dim==3 ){
			mexErrMsgTxt("Does not support 3 dimensional PSFs\n");}
		else{
			/*Allocate memory and assign output pointer*/
			plhs[0] = mxCreateDoubleMatrix(prodDims, 1, mxREAL); /*mxReal is our data-type*/
			plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
			plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
			plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
			plhs[4] = mxCreateDoubleMatrix(k_max+1, 1, mxREAL);
			plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
 
			if(ghxl){
				plhs[6] = mxCreateDoubleMatrix(k_max+1, 1, mxREAL);
				plhs[7] = mxCreateDoubleMatrix(k_max+1, 1, mxREAL);}
			else{
				plhs[6] = mxCreateDoubleMatrix(1, 1, mxREAL);
				plhs[7] = mxCreateDoubleMatrix(1, 1, mxREAL);
			}
			
			if(xl){
				plhs[8] = mxCreateDoubleMatrix(prodDims*(k_max+1), 1, mxREAL);
			}
			else{
				plhs[8] = mxCreateDoubleMatrix( 1, 1, mxREAL);
			}
			
			plhs[9] = mxCreateDoubleMatrix( 1, 1, mxREAL);
			plhs[10] = mxCreateDoubleMatrix( 1, 1, mxREAL);
			plhs[11] = mxCreateDoubleMatrix( 1, 1, mxREAL);

			/* Get a pointer to the data space in our newly allocated memory */
			xk = mxGetPr(plhs[0]);
			fxk = mxGetPr(plhs[1]);
			hxk = mxGetPr(plhs[2]);
			gxk = mxGetPr(plhs[3]);
			fxkl = mxGetPr(plhs[4]);
			k = mxGetPr(plhs[5]);
			hxkl = mxGetPr(plhs[6]);
			gxkl = mxGetPr(plhs[7]);
			xlist = mxGetPr(plhs[8]);
			numGrad = mxGetPr(plhs[9]);
			numBack = mxGetPr(plhs[10]);
			numFunc = mxGetPr(plhs[11]);

			dcopy_(&prodDims,x,&one,xk,&one);
			
			tvreg_gpbb_core(xk,fxk,hxk,gxk,fxkl,k,A,b,alpha,tau,gamma_rel,k_max,D,ctype,d,c,ghxl,xl,hxkl,gxkl,xlist,K,beta,sigma,verbose,numGrad,numBack,numFunc);
		}
	}
}
