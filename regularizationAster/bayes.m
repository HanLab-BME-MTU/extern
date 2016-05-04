% Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
% [covmp,mmap]=bayes(G,mprior,covm,d,covd)
%
%  Given a linear inverse problem Gm=d, a prior mean mprior and covariance
%  matrix covm, data d, and data covariance matrix covd, this function 
%  computes the MAP solution and the corresponding covariance matrix.
%
function [covmp,mmap]=bayes(G,mprior,covm,d,covd)
covmp=inv((G'/covd)*G+inv(covm));
%
% This takes care of any lack of symmetry in covmp.
%
covmp=(covmp+covmp')/2;
covd12=sqrtm(inv(covd));
covm12=sqrtm(inv(covm));
A=[covd12*G; covm12];
rhs=[covd12*d; covm12*mprior];
mmap=A\rhs;
