% Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
% x=kac(A,b,tolx,maxiter) 
%
% Implements Kaczmarz's algorithm to solve a system of equations iteratively.
%
% Input Parameters:
%   A - Constraint matrix.
%   b - right hand side.
%   tolx - difference tolerence for successive iterations (stopping criteria).
%   maxiter - maximum iterations (stopping criteria).
%
% Output Parameters:
%      x - solution.
function x=kac(A,b,tolx,maxiter)

% First, find the size of the matrix.
[m,n]=size(A);

% Make a copy of A' to speed up some accesses.
AT=A';

% Setup an initial solution of all zeros.
x=zeros(n,1);
iter=0;

%  Precompute the row norms squared.
n2=zeros(m,1);
for i=1:m
  n2(i)=norm(AT(:,i))^2;
end

% The main loop performs iterations of Kaczmarz algorithm until 
% maxiters is exceeded or successive iterates differ by less 
% than tolx.  
while (iter <= maxiter)
  % Update the iteration count.
  iter=iter+1;

  %  Start the update cycle with the current solution.
  newx=x;

  %  Perform a cycle of m updates.
  for i=1:m
    newx=newx-((newx'*AT(:,i)-b(i))/(n2(i)))*AT(:,i);
  end

  %  Check for convergence to fixed solution.
  if (norm(newx-x)/(1+norm(x)) < tolx)
    x=newx;
    return;
  end

  %  Update x for the next major iteration.
  x=newx;
end

% If no convergence
disp('Max iterations exceeded.');
