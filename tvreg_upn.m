function [xkp1 fxkp1 hxkp1 gxkp1 fxkp1l info hxkp1l gxkp1l xlist] = ...
    tvreg_upn(A,b,alpha,tau,dims,constraint,opt)
% 
% Solve the problem 
%
% min f(x) = h(x) + g(x) = alpha* TV(x,tau) + 1/2 ||Ax-b||_2^2
% s.t. x in Q
%
% with implicit definitions of h and g.
%
% The function h(x)=TV(X,tau) is a smooth approximation of the TV.
% In fact TV(x,tau) is the Huber functional
% 
% TV(x,tau) = sum_{i,j,l}^{m,n,l} Phi_tau(D_ijl x)
% with
%
%              { ||y||_2 - tau/2   if  ||y||_2>= tau 
% Phi_tau(y) = {
%              { ||y||_2^2 / (2 tau) else
%    
% and D_ijl a finite difference matrix computing the approximated 
% gradient at coordinate (i,j,l).
%
% Input definitions:
% A:     1. A sparse matrix or 
%        2. A struct with the following fields
%           PSF = a matrix describing the point spread function. 
%                  Should be small
%           center = the center of the PSF.
%
%           In this case A in the problem definition is the matrix 
%           formed using the PSF as spartial invariant with reflective 
%           boundary conditions. The PSF should be small for fast 
%           algorithm, but otherwise there is no requirement on the PSF. 
%           TVReg only supports 2 dimensional PSFs, and dims should 
%           only contain 2 elements.
%
%        3. A function handle that takes two arguments
%           - A(x, 0) should return [m, n], where m and n are the
%             dimensions of A
%           - A(x, 1) should return the result of A*x
%           - A(y, 2) should return the result of A'*y
%
% b:     Observed data
%
% alpha: Regularization parameter
%
% tau:   Smoothing parameter of the TV. Suggestion tau =
%        1e-4*norm(X0(:),'inf') where XO is the true X. 
%         
% dims:  dims=[m,n,l] or dims=[m,n]. A vector describing the
%        dimensionality of x, e.g. dims =  [256,256] in the case 
%        of two dimensional problems.
%
% constraint: A struct with the following fields,
%             type = Constraint type
%             c, d = lower and upper bound on x if type=2, where c and d 
%             has the same dimensions as x.
%
%     type == 1
%          Q = R the reals (no constraints)
% 
%     type == 2
%          Q = {x | d_i =>x_i>= c_i, i=0..prod(dims)}
%
% opt:    (Optional) A struct with one or more of the following fields
%          
%         epsb_rel = stopping criteria (default 1e-4)
%
%           The algorithm stops when the iterate y satisfies
%
%              ||G_t(y)||_2/(m*n*l) \leq epsb_rel 
%    
%           where 
%
%               G_t(x) =  1/t (x - P_Q(x - t \nabla f(x)) ) 
%
%           and returns P_Q(x - t \nabla f(x))
%
%         k_max = the maximum number of iterates (default 10000)
% 
%         x0 = initial estimate 
%             (default A'*b or b if A is a struct defining a PSF x = b)
%
%         qs = setting algorithm type (default qs=1)
%            qs = 0 algorithm UPN_0 where mu_k=0.
%              The algorithm UPN applies Nesterov's optimal first order 
%              method for L-Lipschitz continuous gradient. In the
%              implementation estimates bL is used by backtracking
%              line search. In this case its an optimal method for smooth, 
%              non-strongly convex problems. 
%
%            qs = 1 algorithm UPN (default)
%              The algorithm UPN applies Nesterov's optimal first order 
%              method for mu-strongly convex problem with L-Lipschitz 
%              continuous gradient. In the implementation estimates bmu 
%              and bL are used. 
%
%            qs = 2 algorithm Gradient, which is the classic gradient 
%             method with backtracking line search.      
%
%          bL = Initial setting of L_k. Default value is based on an 
%             an estimate of L.
%
%          bmu = Initial setting of mu_k. Default vaule is
%             bmu = min(bL/50,eigmax), where eigmax is the largest
%             eigenvalue of A.
%
% Output definitions
%         Default return values
%         xkp1 = the last iterate
%         fxkp1 = f(xkp1), objective of the last iterate
%         hxkp1 = h(xkp1), the smooth TV of the last iterate
%         gxkp1 = g(xkp1), the fidelity term of the last iterate
%         fxkp1l = a vector containing f(xkp1) for all iterates
%   
%         Optional
%         info = A struct containing additional information on the
%                the behaviour of the algorithm.
%                numFunc = # the objective function is evaluated
%                numGrad = # the gradient function is evaluated
%                numBack = # of backtrackings
%                numRest = # of restarts to reduce \mu_k
%                Lklist = contains the iterates L_k
%                muklist =  constains the iterates mu_k
%                rklist = contains the restart iteration positions.
%
%         hxkp1l = a vector containing h(xkp1) for all iterates
%         gxkp1l = a vector containing g(xkp1) for all iterates 
%         xlist = a matrix containing the iterates xkp1 in all rows.
%                 Make sure only to operate with small x and 
%                 small k_max.
%
%

% Default values.
epsb_rel = 1e-4;
k_max     = 10000;
qs = 1;
verbose = 0;

% If options is given as input, use these values
if nargin > 6
    if isfield(opt,'epsb_rel')
        epsb_rel = opt.epsb_rel;
    end
    if isfield(opt,'k_max')
        k_max = opt.k_max;
    end
    if isfield(opt,'qs')
        qs = opt.qs; 
    end
    if isfield(opt,'verbose')
        verbose = opt.verbose; 
    end

    
    if isfield(opt,'x0')
        x = opt.x0;
    elseif isstruct(A) 
        x = b;
    elseif isa(A, 'function_handle')
        x = A(b, 2);
    else
        x = A'*b;
    end
else
    if isstruct(A) 
        x = b;
    elseif isa(A, 'function_handle')
        x = A(b, 2);
    else
        x = A'*b;
    end

end

clear X
prodDims = prod(dims);
lenDims  = length(dims);

if(numel(dims)==2)
    nDd = 8;
else
    nDd = 12;
end

% Intitial settings of bmu and bL
if isfield(opt,'bL')
    bL = opt.bL;
else
    svdsopts.tol = 1e-5;
    if(isstruct(A))
        eigmax       = svds(A.PSF,1,'L',svdsopts)^2;
        bL           = (alpha*nDd/tau + eigmax)/1e2;
    elseif(isa(A, 'function_handle'))
        V =          @(x) A(A(x,1),2);
        eigmax       = eigs(V, prodDims, 1);
        bL           = (alpha*nDd/tau + eigmax)/1e2;
    else
        eigmax       = svds(A,1,'L',svdsopts)^2;
        bL           = (alpha*nDd/tau + eigmax)/1e2;
    end
end

if isfield(opt,'bmu')
    bmu = opt.bmu;
else
    bmu = min(bL/50,eigmax);
end

if qs == 0 || qs ==2;
    bmu=0;
end


% Initialize vectors to hold tv and fidelity of the iterates
ghxl = 0;
if nargout > 5
    ghxl=1;
end

% Array to hold iterates in columns
xl=0;
if nargout > 7 & k_max*prodDims<1e7
    xl = 1;
end

if(constraint.type==2)
    [xkp1 fxkp1 hxkp1 gxkp1 fxkp1l k hxkp1l gxkp1l xlist numGrad numBack numFunc numRest Lklist muklist rklist] = tvreg_upn_c(A,b,alpha,tau,dims,bL, ...
                                                  bmu,epsb_rel,k_max,x, ...
                                                  constraint.type, ...
                                                  constraint.d, ...
                                                  constraint.c,ghxl,xl, ...
                                                      qs,verbose);
    
else
    [xkp1 fxkp1 hxkp1 gxkp1 fxkp1l k hxkp1l gxkp1l xlist numGrad numBack numFunc numRest Lklist muklist rklist] = tvreg_upn_c(A,b,alpha,tau,dims,bL, ...
                                                  bmu,epsb_rel,k_max,x, ...
                                                  constraint.type, ...
                                                  0, 0,ghxl,xl, ...
                                                      qs,verbose);
end

    
if(~ghxl)
    clear gxkp1l;
    clear hxkp1l;
end

if(~xl)
    clear xlist;
else
    xlist =reshape(xlist,prodDims,k_max+1);
end
    
    
% Truncate excess zeros in fxkp1l, xlist, gxkp1l and hxkp1l
fxkp1l = fxkp1l(1:k+1);

if nargout > 8 & k_max*prodDims<1e7
        xlist = xlist(:,1:k+1);
end

if nargout >5
    info.numFunc = numFunc;
    info.numGrad = numGrad;
    info.numBack = numBack;
    info.numRest = numRest;
    info.Lklist = Lklist(1:k+1);
    info.muklist = muklist(1:k+1);
    info.rklist = rklist;
end

if nargout > 6
    gxkp1l      = gxkp1l(1:k+1);
    hxkp1l      = hxkp1l(1:k+1);
end



% If the iteration counter reaches the k_max, the algorithm did not
% converge
if k == k_max 
    disp('Did not find a epsb_rel solution in k_max iterations.')
end