% fastDeLong.m
% For ROC analysis, implements a fast rank-based algorithm for the 
% Mann-Whitney AUC estimator and for the DeLong covariance matrix estimator. 
% References:
%
% E. R. DeLong, D. M. DeLong, and D. L. Clarke-Pearson, "Comparing the
% areas under two or more correlated receiver operating characteristic
% curves: A nonparametric approach," Biometrics, vol. 44, no. 3,
% pp. 837-845, Sept. 1988.  
% 
% X. Sun and W. Xu, "Fast implementation of DeLong's algorithm for
% comparing the areas under correlated receiver operating characteristic
% curves", IEEE Signal Processing Letters, vol. 21, no. 11, pp. 1389-1393,
% Nov. 2014.
%
% Note: The meanings of x and y are flipped relative to the above references.
%
% Inputs:
%   X (q x m matrix of class 1 ratings)
%   Y (q x n matrix of class 2 ratings)    
% 
% Outputs: 
%   AUC (q x 1 vector of AUC estimates)
%   S (q x q covariance matrix) 
%
% Adam Wunderlich 9/29/2014

function [AUC,S] = fastDeLong(X,Y)

[q,m] = size(X);
[q,n] = size(Y);
V10 = zeros(m,q);
V01 = zeros(n,q);
xi = (n+1):(m+n); % x indices in z 
yj = 1:n; % y indices in z
for k=1:q
   x = X(k,:);  
   y = Y(k,:);
   z = [y x]; 
   TX = tiedrank(x);
   TY = tiedrank(y);
   TZ = tiedrank(z);  
   %compute the structural components
   V10(:,k) = 1 - (TZ(xi) - TX)/n;
   V01(:,k) = (TZ(yj) - TY)/m;
end
AUC = mean(V01)';
S = cov(V10)/m + cov(V01)/n;



