% Parameter Estimation and Inverse Problems, 2nd edition, 2011 
% by R. Aster, B. Borchers, C. Thurber
% generate an m by n array of exponentially distributed random variable
% realizations with expected value mu
function r=exprand(mu,m,n)

% Generate uniform random values, and apply the exponential inverse CDF.
r = -mu * log(rand(m,n));
