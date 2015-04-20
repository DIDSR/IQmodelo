% binPropCov.m
% Estimate covariance matrix for a vector of binomial proportions using an
% unbiased estimator, as described in the following reference:
%
% F. Noo, A. Wunderlich, D. Heuscher, K. Schmitt, Z. Yu, “A Nonparametric Approach for Statistical
% Comparison of Results from Alternative Forced Choice Experiments,” Proc. of SPIE Medical Imaging
% Conference, Vol. 8673, 86730F, Feb. 2013.
%
% Note: This function can be used for analyzing MRMC studies with fixed
% readers.  For random-reader analysis, a function is provided in the
% iMRMC_Binary package within the iMRMC project on google code.  
%
% Input: X (n x q matrix of success scores, where n=number of cases and 
% q=number of fixed imaging scenarios or readers)
% Outputs: p (q x 1 vector of percent correct estimates)
% S (q x q covariance matrix estimate)
%
% Adam Wunderlich 5/22/2014

function [p,S] = binPropCov(X)

[n,q] = size(X);
p = mean(X)';
P = X'*X./n;
S = (P-p*p')./(n-1);