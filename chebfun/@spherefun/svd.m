function varargout = svd( f )
%SVD    Singular value decomposition of a SPHEREFUN.
%
%   SVD(F) returns the singular values of F. The number of singular values
%   returned is equal to the rank of the SPHEREFUN.
%
%   S = SVD(F) returns the singular values of F. S is a vector of singular
%   values in decreasing order.
%
%   [U, S, V] = SVD(F) returns the SVD of F. V is a quasi-matrix of
%   orthogonal CHEBFUN objects, U is a quasimatrix of CHEBFUN objects that
%   are orthogonal with respect to the sin(th) weight (derived from the
%   measure on the sphere) and S is a diagonal matrix with the singular
%   values on the diagonal.
%
%   The length and rank of a SPHEREFUN are slightly different quantities.
%   LENGTH(F) is the number of pivots used by the constructor, and
%   RANK(F) is the number of significant singular values of F. The relation
%   RANK(F) <= LENGTH(F) should always hold.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    varargout = {[]};
    return
end

% Get CDR decomposition of f:
[C, D, R] = cdr(f);

% Do QR in both variables, one with the weighted inner-product and one with
% the standard L2 inner-product.
[QwC, RwC] = sphereQR(C); 
[QwR, RwR] = qr(R);

% Use the QR factorizations of the columns and rows to make up the SVD of
% the SPHEREFUN object.  Since
%
%        C * D * R = QwC * ( RwC * D * RwR' ) * QwR'
%
% we compute the SVD of ( RwC * D * RwR' ).
[U, S, V] = svd(RwC * D * RwR.');
U = QwC * U;
V = QwR * V;

% Output just like the svd of a matrix.
if ( nargout > 1 )
    varargout = {U, S, V};
else
    varargout = {diag(S)};
end

end


function [Q, R] = sphereQR( A )
% Fast version of the abstractQR code, specifically for the sphere. 

n = max(2*length(A), 30);         % Can probably get away with a smaller n.
[x, w] = legpts(n, [0, pi]);      % Legendre points

% Do a weighted QR, and then unweight the QR: 
WR = spdiags(sqrt(w.' .* sin(x)), 0, n, n);
invWR = spdiags(1./sqrt(w.' .* sin(x)), 0, n, n);

% Discrete QR with inner product <u,v> = sum(sin(th)*conj(u).*v):
[discreteQ, discreteR] = qr(WR*A(x, :), 0);

s = sign(diag(discreteR));    % }
s(~s) = 1;                    %  } Enforce diag(R) >= 0
S = diag(s);                  % }

% Undo the weighting: 
discreteQ = invWR*discreteQ*S;

% Correct the signs in R: 
discreteR = S*discreteR;

% Go back to continuous land: 
Q = chebfun(legvals2chebvals(discreteQ), [0, pi]); 
R = discreteR; 

end