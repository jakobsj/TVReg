function [xk fxk hxk gxk fxkl info hxkl gxkl xlist] = ...
    tvreg_gpbb(A,b,alpha,tau,dims,constraint,opt)
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
%        3. A function handle that takes two arguments
%           - A(x, 0) should return [m, n], where m and n are the
%             dimensions of A
%           - A(x, 1) should return the result of A*x
%           - A(y, 2) should return the result of A'*y
%
% b:      Observed data
%
% alpha:  Regularization parameter
%
% tau:    Smoothing parameter of the TV. Suggestion tau =
%         1e-4*norm(X0(:),'inf') where XO is the true X. 
%         
% dims:   dims=[m,n,l] or dims=[m,n]. A vector describing the
%         dimensionality of x, e.g. dims =  [256,256] in the case 
%         of two dimensional problems.
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
%              G_t(x) =  1/t (x - P_Q(x - t \nabla f(x)) ) 
%
%           and returns P_Q(x - t \nabla f(x))
%
%         k_max = the maximum number of iterates (default 10000)
% 
%         x0 = initial estimate 
%             (default A'*b or b if A is a struct defining a PSF x = b)
%
%         K = look back length for sufficient decrease reference value 
%             (default 2)           
%
%         sigma = sufficient decrease parameter (default 0.1)
%
%         beta = reduction parameter of the stepsize in backtracking
%              (default 0.95)
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
%         hxkp1l = a vector containing h(xkp1) for all iterates
%         gxkp1l = a vector containing g(xkp1) for all iterates 
%         xlist = a matrix containing the iterates xkp1 in all rows.
%                 Make sure only to operate with small x and 
%                 small k_max.
%
%
% The algorithm applies steepest descent (the gradient projections 
% algorithm) with Barzilai and Borwein strategy and nonmonotonic line
% search.
%

% Default values.
epsb_rel = 1e-4;
k_max     = 10000;

K     = 2;
beta  = 0.95;
sigma = 0.1;
verbose = 0;

% If options is given as input, use these values
if nargin > 6
    if isfield(opt,'epsb_rel')
        epsb_rel = opt.epsb_rel;
    end
    if isfield(opt,'k_max')
        k_max = opt.k_max;
    end
    if isfield(opt,'K')
        K = opt.K;
    end    
    if isfield(opt,'beta')
        beta = opt.beta;
    end    
    if isfield(opt,'sigma')
        sigma = opt.sigma;
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

    
% Initialize vectors to hold tv and fidelity of the iterates
ghxl = 0;
if nargout > 6
    ghxl=1;
end

% Array to hold iterates in columns
xl=0;
if nargout > 8 & k_max*prodDims<1e7
    xl = 1;
end

    
[xk fxk hxk gxk fxkl k hxkl gxkl xlist numGrad numBack numFunc] = tvreg_gpbb_c(A,b,alpha,tau,dims,epsb_rel,k_max,x, ...
                                                  constraint.type, ...
                                                  constraint.d,constraint.c,ghxl,xl,K,beta,sigma,verbose);
if(~ghxl)
    clear gxkl;
    clear hxkl;
end

if(~xl)
    clear xlist;
else
    xlist =reshape(xlist,prodDims,k_max+1);
end
    
% Truncate excess zeros in fxkl, xlist, gxkl and hxkl

fxkl = fxkl(1:k+1);


if nargout > 8 & k_max*prodDims<1e7
        xlist = xlist(:,1:k+1);
end

if nargout >5
    info.numFunc = numFunc;
    info.numGrad = numGrad;
    info.numBack = numBack;
end

if nargout > 6
    hxkl      = hxkl(1:k+1);
    gxkl      = gxkl(1:k+1);
end

% If the iteration counter reaches the k_max, the algorithm did not
% converge
if k == k_max 
        disp('Did not find a epsb_rel solution in k_max iterations.')
end

