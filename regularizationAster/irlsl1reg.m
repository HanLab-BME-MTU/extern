% Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
% mreg=irlsl1reg(G,d,L,alpha,maxiter,tolx,tolr);
%
% Solves the system Gm=d using sparsity regularization on Lm.
%
% Solves the L1 regularized least squares problem
%
%  min norm(G*m-d,2)^2+alpha*norm(L*m,1)
%
% using iteratively reweighted least squares.
%
% Inputs:
%       G,d,L,alpha            Problem data.
%       maxiter                Maximum number of IRLS iterations.
%       tolx                   Tolerance on successive iterates.
%       tolr                   Tolerance below which we consider an 
%                              element of L*m to be effectively
%                              zero.
% 
% default values are provided for maxiter, tolx, and tolr if they
% are not specified.
%
function mreg=irlsl1reg(G,d,L,alpha,maxiter,tolx,tolr)

% Default for tolr=1.0e-6
if (nargin < 7)
  tolr=1.0e-6;
end

% Default for tolx=1.0e-4;
if (nargin < 6)
  tolx=1.0e-4;
end

% Default for maxiter=100
if (nargin < 5)
  maxiter=100;
end

% unchanging constants in the system that is solved repeatedly
GTG=G'*G;
GTd=G'*d;

% Start with an initial unweighted solution.
m=(2*GTG+alpha*(L'*L))\(2*G'*d);
iter=1;

% iterate until maxiter or we converge and return
while (iter < maxiter)
  iter=iter+1;

  % get get the magnitude of Lm, but don't let any element be less than tolr
  absLm=abs(L*m);
  absLm(absLm<tolr)=tolr;

  % build the diagonal weighting matrix for this iteration
  R=diag(1./absLm);

  mold=m;
  % get the new iterate and check for convergance
  m=(2*GTG+alpha*L'*R*L)\(2*GTd);
  if (norm(m-mold)/(1+norm(mold)) < tolx)
    mreg=m;
    return
  end
end

% Give a warning, if desired, but return best solution.
%warning('irlslreg1 maximum iterations exceeded.');
mreg=m;
