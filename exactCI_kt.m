% exactCI_kt.m
% For a linear observer defined by a known (fixed) template, 
% returns exact 1-(alpha1+alpha2) confidence intervals for SNR and AUC as described in
%
% A. Wunderlich and F. Noo, "Confidence Intervals for Performance
% Assessment of Linear Observers," Medical Physics, vol. 38, no. S1, 
% pp. S57-S68, July 2011.  
%
% Inputs: alpha1 (lower significance level), alpha2 (upper significance level) 
% (for conventional two-sided 95% intervals, set alpha1=alpha2=0.025), 
% x (1 x m vector of class 1 ratings), y (1 x n vector of class 2 ratings) 
%
% Output: ret (structure containing point and interval estimates of SNR and AUC)
%
% Adam Wunderlich 5/22/2014

function  [ret] = exactCI_kt(alpha1,alpha2,x,y)

   m = length(x);
   n = length(y);
   xbar = mean(x);
   ybar = mean(y);
   s1_sq = var(x);
   s2_sq = var(y);
   s_sq = ((m-1)*s1_sq + (n-1)*s2_sq)/(m+n-2);
   s = sqrt(s_sq);
   gam = sqrt(2*pi/(m+n-2))/beta((m+n-3)/2,0.5);
   snr = gam*(ybar-xbar)/s;
   auc = normcdf(snr/sqrt(2));

   nu = m+n-2;
   eta = sqrt(m*n/(m+n))/gam;
   t = eta*snr;
   % find confidence interval for noncentral t noncentrality parameter
   % Note: The user may wish to check EXITFLAG and add appropriate error messages, or they may wish to
   % modify the tolerances of the fzero function with the OPTIONS argument.  
   delta0 = 1;  % initial guess
   if (alpha1 ~= 0 && alpha2 ~= 0),
      [delta_L,FVAL,EXITFLAG] = fzero(@(delta) nctcdf(t,nu,delta)-1+alpha1,delta0);
      [delta_U,FVAL,EXITFLAG] = fzero(@(delta) nctcdf(t,nu,delta)-alpha2,delta0);
   elseif (alpha1 ~= 0 && alpha2 == 0),
      [delta_L,FVAL,EXITFLAG] = fzero(@(delta) nctcdf(t,nu,delta)-1+alpha1,delta0);
      delta_U = Inf;
   elseif (alpha1 == 0 && alpha2 ~= 0),
      delta_L = -Inf;
      [delta_U,FVAL,EXITFLAG] = fzero(@(delta) nctcdf(t,nu,delta)-alpha2,delta0);
   else
      disp('Warning: Both alpha1 and alpha2 are zero!')
      delta_U = Inf;
      delta_L = -Inf;
   end
   % find confidence intervals for summary measures
   SNR_CI = zeros(1,2);
   SNR_CI(1) = delta_L*sqrt((m+n)/(m*n));
   SNR_CI(2) = delta_U*sqrt((m+n)/(m*n));
   AUC_CI = normcdf(SNR_CI/sqrt(2));
   ret.m = m;
   ret.n = n;
   ret.alpha1 = alpha1;
   ret.alpha2 = alpha2;
   ret.SNR = snr;
   ret.SNR_CI = SNR_CI;
   ret.AUC = auc;
   ret.AUC_CI = AUC_CI;

