% binProp_CI.m
% Returns a 1-(alpha1+alpha2) confidence interval for a single binomial 
% proportion or for a difference of binomial proportions.  This function 
% requires the function binPropCov.m, and assumes that variabilty is due to
% cases only.  
%
% The confidence interval for a single proportion is computed using the logit
% transformation method recommended on page 22 in the book Pepe M.S.,"The statistical
% evaluation of medical tests for classification and prediction", Oxford
% Univ. Press, 2003.  For a difference of porportions, the usual normal
% approximation is used.  
%
% Inputs:  alpha1 (lower significance level), alpha2 (upper significance level) 
% (for conventional two-sided 95% intervals, set alpha1=alpha2=0.025), 
% X (n x q matrix of success scores, where n=number of cases and 
% q = 1 or 2, the number of fixed imaging scenarios or readers.
%
% Outputs: p (point estimates for each proportion), 
% p_CI (confidence interval for single binomial proportion (q=1) or
% for a difference (q=2))
% Note: The interval estimate for a difference is for
% performance of the second scenario minus the first.  
% 
% Adam Wunderlich
% 5/21/2014

function [p,p_CI] = binProp_CI(alpha1,alpha2,X)

[p,S] = binPropCov(X); 
[n,q] = size(X);
p_CI = zeros(1,2);
if q==1,
    % use logit transformation method
    if (alpha1~=0 && alpha2~=0)
       logit_CI(1) = log(p/(1-p))+norminv(alpha1)*sqrt(S)/(p*(1-p));
       logit_CI(2) = log(p/(1-p))+norminv(1-alpha2)*sqrt(S)/(p*(1-p));
       p_CI = exp(logit_CI)./(1+exp(logit_CI));
    elseif (alpha1==0 && alpha2~=0)
       logit_CI(2) = log(p/(1-p))+norminv(1-alpha2)*sqrt(S)/(p*(1-p));
       p_CI(1) = 0;
       p_CI(2) = exp(logit_CI(2))./(1+exp(logit_CI(2)));
    elseif (alpha1~=0 && alpha2==0)
       logit_CI(1) = log(p/(1-p))+norminv(alpha1)*sqrt(S)/(p*(1-p));
       p_CI(1) = exp(logit_CI(1))./(1+exp(logit_CI(1)));
       p_CI(2) = 1;
    end
else % q=2
   % use usual normal approximation
   pdiff = p(2)-p(1);
   var_diff = S(1,1) + S(2,2) - 2*S(1,2);
   if alpha1==0,
       p_CI(1) = -1;
   else
       p_CI(1) = pdiff + norminv(alpha1)*sqrt(var_diff);
   end
   if alpha2==0,
       p_CI(2) = 1;
   else
       p_CI(2) = pdiff + norminv(1-alpha2)*sqrt(var_diff);
   end
end
