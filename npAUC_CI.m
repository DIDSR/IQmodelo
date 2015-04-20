% npAUC_CI.m
% For an ROC assessment, returns a 1-(alpha1+alpha2) confidence interval 
% for a single AUC or for a difference of AUCs.  This function requires the 
% function fastDeLong.m, and assumes that variabilty is due to cases only.  
%
% The confidence interval for a single AUC is computed using the logit
% transformation method recommended on page 107 in the book: 
% Pepe M.S.,"The statistical evaluation of medical tests for classification 
% and prediction", Oxford Univ. Press, 2003.  For a difference of AUCs, 
% the usual normal approximation is used.  
%
% Inputs:  alpha1 (lower significance level), alpha2 (upper significance level) 
% (for conventional two-sided 95% intervals, set alpha1=alpha2=0.025), 
% X (q x m matrix of class 1 ratings)
% Y (q x n matrix of class 2 ratings)
% where q = 1 or 2, the number of fixed imaging scenarios or readers
%
% Outputs: AUC (AUC point estimates), AUC_CI (confidence interval for 
% a single AUC (q=1) or for a difference (q=2))
% Note: The interval estimate for a difference is for
% performance of the second scenario minus the first.  
% 
% Adam Wunderlich
% 9/26/2014

function [AUC,AUC_CI] = npAUC_CI(alpha1,alpha2,X,Y)

[q,n] = size(Y);
[AUC,S] = fastDeLong(X,Y);
AUC_CI = zeros(1,2);
if q==1,
    % use logit transformation method
    if (alpha1~=0 && alpha2~=0)
       logit_CI(1) = log(AUC/(1-AUC))+norminv(alpha1)*sqrt(S)/(AUC*(1-AUC));
       logit_CI(2) = log(AUC/(1-AUC))+norminv(1-alpha2)*sqrt(S)/(AUC*(1-AUC));
       AUC_CI = exp(logit_CI)./(1+exp(logit_CI)); 
    elseif (alpha1==0 && alpha2~=0)
       logit_CI(2) = log(AUC/(1-AUC))+norminv(1-alpha2)*sqrt(S)/(AUC*(1-AUC));
       AUC_CI(1) = 0;
       AUC_CI(2) = exp(logit_CI(2))./(1+exp(logit_CI(2)));
    elseif (alpha1~=0 && alpha2==0)
       logit_CI(1) = log(AUC/(1-AUC))+norminv(alpha1)*sqrt(S)/(AUC*(1-AUC));
       AUC_CI(1) = exp(logit_CI(1))./(1+exp(logit_CI(1)));
       AUC_CI(2) = 1;
    end
else % q=2
   % use usual normal approximation
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
