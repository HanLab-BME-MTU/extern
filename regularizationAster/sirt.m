% Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
% x=sirt(A,b,tolx,maxiter) 
%
% Implements the SIRT algorithm to solve a system of equations iteratively.
%
% Input Parameters:
%   A       - Constraint matrix.
%   b       - right hand side.
%   tolx    - difference tolerance for successive iterations (stopping criteria)
%   maxiter - maximum iterations (stopping criteria).
%
% Output Parameters:
%      x - solution.
%
function x=sirt(A,b,tolx,maxiter)

% First, get the size of A.
[m,n]=size(A);

% Alpha is a damping factor.  If alpha<1, then we won't take full steps
% in the SIRT direction.  Using a smaller value of alpha (say alpha=.75)
% can help with convergence on some problems.  
alpha=1.0;

% In the A1 array, we convert all nonzero entries in A to +1.  
A1=(A>0);

% Get transposed copies of the arrays for faster access.
AT=A';
A1T=A1';

% Start with the zero solution.
x=zeros(n,1);

%  Precompute N(i) and L(i) factors.
N=zeros(m,1);
L=zeros(m,1);
NRAYS=zeros(n,1);

for i=1:m
  N(i)=sum(A1T(:,i));
  L(i)=sum(AT(:,i));
end

for i=1:n
  NRAYS(i)=sum(A1(:,i));
end

% Start the iteration count at 0.
iter=0;

% Now, the main loop, don't loop more than maxiter times
while (iter<=maxiter)
  iter=iter+1;

  % Start the next round of updates with the current solution.
  newx=x;

  %  Now, compute the updates for all of the rays and all cells, and put
  %  them in a vector called deltax.  
  deltax=zeros(n,1);
  for i=1:m,
    %  Compute the approximate travel time for ray i.
    q=A1T(:,i)'*newx;

    % We use the following more accurate formula for delta.
    delta=b(i)/L(i)-q/N(i);

    % This formula is less accurate and doesn't work nearly as well.
    %
    %    delta=(b(i)-q)/N(i);
    %
 
    %  Perform updates for those cells touched by ray i.
    deltax=deltax+delta*A1T(:,i);
  end

  %  Now, add the average of the updates to the old model.
  %  Note the "./" here.  This means that the averaging is done with
  %  respect to the number of rays that pass through a particular cell!
  newx=newx+alpha*deltax./NRAYS;
  
  % Check for convergence
  if (norm(newx-x)/(1+norm(x)) < tolx)
    x=newx
    return;
  end

  % No convergence, so setup for the next major iteration.
  x=newx;
end
disp('Max iterations exceeded.');
