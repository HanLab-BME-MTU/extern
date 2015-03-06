% Parameter Estimation and Inverse Problems, 2nd edition, 2011 
% by R. Aster, B. Borchers, C. Thurber
% x=art(A,b,tolx,maxiter) 
%
% Implements the ART algorithm to solve a system of equations iteratively.
%
% Input Parameters:
%   A       - Constraint matrix.
%   b       - right hand side.
%   tolx    - difference tolerance for successive iterations (stopping criteria)
%   maxiter - maximum iterations (stopping criteria).
%
% Output Parameters:
%   x - solution.
function x=art(A,b,tolx,maxiter)

% Alpha is a damping factor.  If alpha<1, then we won't take full steps
% in the ART direction.  Using a smaller value of alpha (say alpha=.75)
% can help with convergence on some problems.  
alpha=1.0;

% First, get the size of A.
[m,n]=size(A);

% In the A1 array, we convert all nonzero entries in A to 1.  
A1=(A>0);

% Get transposed copies of the arrays for faster access.
AP=A';
A1P=A1';

%  Precompute N(i) and L(i) factors.
N=zeros(m,1);
L=zeros(m,1);
for i=1:m
  N(i)=sum(A1(i,:));
  L(i)=sum(A(i,:));
end

% Start with the zero solution.
x=zeros(n,1);

% Start the iteration count at 0.
iter=0;

% Now, the main loop.
while (true)
  % Check to make sure that we haven't exceeded maxiter.
  iter=iter+1;
  if (iter > maxiter)
      disp('Max iterations exceeded.');
      x=newx;
      return;
  end

  % Start the next round of updates with the current solution.
  newx=x;
  
  % Now, update each of the m constraints.
  for i=1:m
    %  Compute the weighted sum for constraint i.
    q=A1P(:,i)'*newx;

    % We use the more accurate formula for delta.
    delta=b(i)/L(i)-q/N(i);

    % This alternative formula is less accurate and doesn't work nearly as well.
    %
    %    delta=(b(i)-q)/N(i);
    %

    % Now do the update.
    newx=newx+alpha*delta*A1P(:,i);
  end

  % Check for convergence
  if (norm(newx-x)/(1+norm(x)) < tolx)
    x=newx;
    return;
  end

  % No convergence, so setup for the next ART iteration.
  x=newx;
end
