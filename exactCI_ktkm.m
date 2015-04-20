% exactCI_ktkm.m
% For a linear observer defined by a known (fixed) template, and known difference of class means,
% returns exact 1-(alpha1+alpha2) confidence intervals for SNR and AUC as described in
%
% A. Wunderlich and F. Noo, "On Efficient Assessment of Image Quality Metrics Based on Linear Model 
% Observers," IEEE Transactions on Nuclear Science, vol. 59, no. 3, pp. 568-578, June 2012.
%
% Inputs: alpha1 (lower significance level), alpha2 (upper significance level) 
% (for conventional two-sided 95% intervals, set alpha1=alpha2=0.025), 
% delta (known difference of class means), x (1 x m vector of class 1 ratings), 
% y (1 x n vector of class 2 ratings) 
%
% Output: ret (structure containing point and interval estimates of SNR and AUC)
%
% Adam Wunderlich 5/22/2014

function  [ret] = exactCI_ktkm(alpha1,alpha2,delta,x,y)

   omega1=alpha1; % notation from paper
   omega2=alpha2; % notation from paper
   m = length(x);
   n = length(y);
   if m+n < 2,
      error('error: m+n too small!')
   end
   q = m+n-1;
   gam = sqrt(2*pi/q)/beta((q-1)/2,0.5);
   xtilde = (sum(x)+sum(y)-n*delta)/(m+n);
   ytilde = (sum(x)+sum(y)+m*delta)/(m+n);
   z = [(x-xtilde*ones(1,m)),(y-ytilde*ones(1,n))];
   stilde_sq = sum(z.^2)/(m+n-1);
   stilde = sqrt(stilde_sq);
   snr = gam*delta/stilde;
   auc = normcdf(snr/sqrt(2));
   snr_sq = snr^2;

   alpha = (m+n-1)/2; % inv-gamma parameter
   eta = pi/(beta((m+n-2)*.5,.5)^2); 
   % find confidence interval for second parameter of inverted gamma
   % Note: The user may wish to check EXITFLAG and add appropriate error messages, 
   % or they may wish to modify the tolerances of the fzero function with the OPTIONS argument.  
   beta0 = [1e-6 1e6];  % search interval
   if (omega1 ~= 0 && omega2 ~= 0),
      [beta_L,FVAL,EXITFLAG] = fzero(@(beta) gammainc(beta/snr_sq,alpha,'upper')-(1-omega1),beta0);
      [beta_U,FVAL,EXITFLAG] = fzero(@(beta) gammainc(beta/snr_sq,alpha,'upper')-omega2,beta0);  
   elseif (omega1 == 0 && omega2 ~= 0),
      beta_L = 0;
      [beta_U,FVAL,EXITFLAG] = fzero(@(beta) gammainc(beta/snr_sq,alpha,'upper')-omega2,beta0);
   elseif (omega1 ~=0 && omega2 == 0),    
      [beta_L,FVAL,EXITFLAG] = fzero(@(beta) gammainc(beta/snr_sq,alpha,'upper')-(1-omega1),beta0);
      beta_U = Inf;
   else
      disp('Warning: Both omega1 and omega2 are zero!')
      beta_L = 0;
      beta_U = Inf;
   end
   % find confidence intervals for summary measures
   SNR_CI = zeros(1,2);
   SNR_CI(1) = sqrt(beta_L/eta);
   SNR_CI(2) = sqrt(beta_U/eta);
   AUC_CI = normcdf(SNR_CI/sqrt(2));
   ret.m = m;
   ret.n = n;
   ret.alpha1 = alpha1;
   ret.alpha2 = alpha2;
   ret.SNR = snr;
   ret.SNR_CI = SNR_CI;
   ret.AUC = auc;
   ret.AUC_CI = AUC_CI;

