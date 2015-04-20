% EROCcov.m
% Computes U-statistic estimates of the area under ROC, LROC, or 
% EROC curves, and the corresponding Sen-type covariance matrix
% estimate for q fixed imaging scenarios as described in the following
% references:
%
% E. R. DeLong, D. M. DeLong, and D. L. Clarke-Pearson, "Comparing the
% areas under two or more correlated receiver operating characteristic
% curves: A nonparametric approach," Biometrics, vol. 44, no. 3,
% pp. 837-845, Sept. 1988.  
%
% A. Wunderlich, F. Noo, “A nonparametric procedure for comparing the areas 
% under correlated LROC curves,” IEEE Transactions on Medical Imaging, 
% vol. 31, no. 11, pp. 2050-2061, Nov. 2012.
%
% A. Wunderlich and B. Goossens, “Nonparametric estimation receiver 
% operating characteristic analysis for performance evaluation on combined 
% detection and estimation tasks,” SPIE Journal of Medical Imaging, vol. 1, 
% no. 3, 031002, Oct. 2014. 
% 
% Inputs:
%   X (q x m matrix of class 1 ratings)
%   Y (q x n matrix of class 2 ratings)
%   U (q x n matrix of utilities for class 2 images) 
% Note: For ROC analysis, set U = ones(q,n).  For LROC analysis, each entry 
% of U is one if the lesion was localized correctly, and zero otherwise.  
% For EROC analysis, U is a defined by a utility function for parameter 
% estimation, as described in the third reference above.     
% 
% Outputs: 
%   AUC (q x 1 vector of AUC estimates)
%   S (q x q covariance matrix) 
%
% Adam Wunderlich 9/29/2014


function [AUC,S] = EROCcov(X,Y,U)

[q,m] = size(X);
[q,n] = size(Y);
V10 = zeros(m,q);
V01 = zeros(n,q);
for k=1:q
   x=X(k,:);  
   y=Y(k,:);
   u=U(k,:);
   % compute utility-success matrix for modality k
   US = zeros(m,n); 
   for i=1:m
      % Second term is needed for categorical ratings.
      % For continuous-valued ratings, it has no effect.
      US(i,:)=u.*(y > x(i)) + u.*(.5).*(y==x(i)); 
   end  
   %compute structural components
   V10(:,k) = mean(US,2);
   V01(:,k) = mean(US,1);
end
AUC = mean(V10)';
S = cov(V10)/m + cov(V01)/n;



