% npAEROC_CI.m
% For an LROC/EROC assessment, returns an approximate 1-(alpha1+alpha2) 
% confidence interval for a single AUC or for a difference of AUCs.  
% This function requires the function EROCcov.m, and assumes that 
% variabilty is due to cases only.  
%
% Inputs:  alpha1 (lower significance level), alpha2 (upper significance level) 
% (for conventional two-sided 95% intervals, set alpha1=alpha2=0.025), 
% X (q x m matrix of class 1 ratings)
% Y (q x n matrix of class 2 ratings)
% U (q x n matrix of utilities for class 2 images) 
% where q = 1 or 2, the number of fixed imaging scenarios or readers.
% Note: For LROC analysis, each entry of U is one if the lesion was 
% localized correctly, and zero otherwise.  
%
% Outputs: AUC (AUC point estimates), AUC_CI (confidence interval for 
% a single AUC (q=1) or for a difference (q=2))
% Note: The interval estimate for a difference is for
% performance of the second scenario minus the first.  
% 
% Adam Wunderlich
% 6/4/2014

function [AUC,AUC_CI] = npAEROC_CI(alpha1,alpha2,X,Y,U)

[q,n] = size(Y);
[AUC,S] = EROCcov(X,Y,U);
AUC_CI = zeros(1,2);
if q==1,
    AUC_CI(1) = AUC+norminv(alpha1)*sqrt(S);
    AUC_CI(2) = AUC+norminv(1-alpha2)*sqrt(S);
else % q=2
   AUCdiff = AUC(2)-AUC(1);
   var_diff = S(1,1) + S(2,2) - 2*S(1,2);
   if alpha1==0,
       AUC_CI(1) = -1;
   else
      AUC_CI(1) = AUCdiff + norminv(alpha1)*sqrt(var_diff);
   end
   if alpha2==0,
      AUC_CI(2) = 1;
   else
      AUC_CI(2) = AUCdiff + norminv(1-alpha2)*sqrt(var_diff);
   end
end
