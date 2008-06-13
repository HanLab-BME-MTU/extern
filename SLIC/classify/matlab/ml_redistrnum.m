function n=ml_redistrnum(N,k)
%ML_REDISTRNUM Evenly separate a number into several numbers.
%   N = ML_REDISTRNUM(M,K) separated the number M into a vector with the
%   length K evenly. The sum of N is M.
%   
%   Example:
%       ml_redistrnum(10,4) returns [3 3 2 2].

%   ??-???-???? Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%       - add comments
%   22-MAR-2005 Modified T. Zhao
%       - debug k==1
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

if k==1
    n=N;
    return
end

avgn=floor(N/k);
for i=1:k
    n(i)=avgn;
end

remain=N-sum(n);
k=1;
while remain>0
    n(k)=n(k)+1;
    remain=remain-1;
    k=k+1;
end

% 
% n(k)=N-sum(n(1:k-1));
% 
% [maxn,maxpos]=max(n);
% [minn,minpos]=min(n);
% 
% while(maxn-minn>1)
%     n(maxpos)=n(maxpos)-1;
%     n(minpos)=n(minpos)+1;
%     [maxn,maxpos]=max(n);
%     [minn,minpos]=min(n);
% end