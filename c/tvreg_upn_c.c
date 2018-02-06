#include <mex.h>
#include "tools.h"
#include "tv_core.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	register double *b,*d,*c,*x,*dims,alpha,tau,bL,bmu,epsb_rel,*hxkp1l,*gxkp1l,*xlist;
	register double *xkp1,*fxkp1,*hxkp1,*gxkp1,*fxkp1l,*k,*numGrad,*numBack,*numFunc,*numRest,*Lklist,*muklist,*rklist;
	mxArray *M,*S,*Mdims;
	int i,j,k_max,dim,prodDims,one=1,ctype,ghxl,xl,qs,verbose,rql,temp;  
	Atype A;
	Dtype D;
	listelement *p_rp;
    mxArray *dims_lhs;
    mxArray *dims_rhs[3];
    double *dimsp;

	p_rp = NULL;
	p_rp = AddItem(p_rp, 0); /*Need to initialize to avoid pointer problems when passing to the core function*/

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
		bL = (double)(mxGetScalar(S));
		
		S = (mxArray*)prhs[6];
		bmu = (double)(mxGetScalar(S));
		
		S = (mxArray*)prhs[7];
		epsb_rel = (double)(mxGetScalar(S));

		S = (mxArray*)prhs[8];
		k_max = (int)(mxGetScalar(S));

		M = (mxArray*)prhs[9];
		x = mxGetPr(M);

		if(mxIsStruct(prhs[0])){
			A.Anc = mxGetM(M)*mxGetN(M);
		}

		S = (mxArray*)prhs[10];
		ctype = (int)(mxGetScalar(S));

		M = (mxArray*)prhs[11];
		d = mxGetPr(M);

		M = (mxArray*)prhs[12];
		c = mxGetPr(M);

		S = (mxArray*)prhs[13];
		ghxl = (int)(mxGetScalar(S));

		S = (mxArray*)prhs[14];
		xl = (int)(mxGetScalar(S));

		S = (mxArray*)prhs[15];
		qs = (int)(mxGetScalar(S));

		S = (mxArray*)prhs[16];
		verbose = (int)(mxGetScalar(S));

		/*obtain the dimensions */

		dim = MAX( mxGetM(Mdims), mxGetN(Mdims) );
		prodDims=1;
		for(i=0;i<dim;i++){
			prodDims = prodDims*(int)dims[i];}

		D.dim = dim;
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
			plhs[12] = mxCreateDoubleMatrix( 1, 1, mxREAL);
			plhs[13] = mxCreateDoubleMatrix( k_max+1, 1, mxREAL);
			plhs[14] = mxCreateDoubleMatrix( k_max+1, 1, mxREAL);

			/* Get a pointer to the data space in our newly allocated memory */
			xkp1 = mxGetPr(plhs[0]);
			fxkp1 = mxGetPr(plhs[1]);
			hxkp1 = mxGetPr(plhs[2]);
			gxkp1 = mxGetPr(plhs[3]);
			fxkp1l = mxGetPr(plhs[4]);
			k = mxGetPr(plhs[5]);
			hxkp1l = mxGetPr(plhs[6]);
			gxkp1l = mxGetPr(plhs[7]);
			xlist = mxGetPr(plhs[8]);

			numGrad = mxGetPr(plhs[9]);
			numBack = mxGetPr(plhs[10]);
			numFunc = mxGetPr(plhs[11]);
			numRest = mxGetPr(plhs[12]);
			Lklist = mxGetPr(plhs[13]);
			muklist = mxGetPr(plhs[14]);

			dcopy_(&prodDims,x,&one,xkp1,&one);
			
			tvreg_upn_core(xkp1,fxkp1,hxkp1,gxkp1,fxkp1l,k,A,b,alpha,tau,bL,bmu,epsb_rel,k_max,D,ctype,d,c,ghxl,xl,hxkp1l,gxkp1l,xlist,qs,verbose,numGrad,numBack,numFunc,numRest,Lklist,muklist,p_rp);

			/*write the dynamical allocated restart list to a vector with the correct dimensions*/			
			rql = QueueLength(p_rp);


			plhs[15] = mxCreateDoubleMatrix( rql-1, 1, mxREAL);
			rklist = mxGetPr(plhs[15]);

			WriteQueueData(p_rp,rklist,rql-1); /*write the list, minus the intial*/
			ClearQueue(p_rp);
		}
	}

}
