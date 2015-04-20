% exactCI_CHO_km.m
% For a CHO with known difference of class means, returns exact 
% 1-(alpha1+alpha2) confidence intervals for SNR and AUC as described in
%
% A. Wunderlich and F. Noo, "New Theoretical Results On Channelized Hotelling 
% Observer Performance Estimation with Known Difference of Class Means," 
% IEEE Transactions on Nuclear Science, vol. 60, no. 1, pp. 182-193, Feb. 2013.  
%
% Since there is only a marginal difference between the estimators for 
% scenarios 1 and 2 described in the above paper, the more general scenario 2 is assumed.  
% Below, m, n, and p denote the number of class 1 images, class 2 images, and channels, respectively.
%
% Inputs: alpha1 (lower significance level), alpha2 (upper significance level) 
% (for conventional two-sided 95% intervals, set alpha1=alpha2=0.025), 
% delta_mu (p x 1 vector corresponding to known difference of class means in channel space)
% v1 (p x m matrix of class 1 channel output vectors)
% v2 (p x n matrix of class 2 channel output vectors)
%
% Output: ret (structure containing point and interval estimates of SNR and AUC
%
% Adam Wunderlich 5/21/2014

function  [ret] = exactCI_CHO_km(alpha1,alpha2,delta_mu,v1,v2)

[p,m] = size(v1);
[p,n] = size(v2);
omega1=alpha1; % notation from paper
omega2=alpha2; % notation from paper
scenario = 2; % assume scenario 2
l = m+n-p-scenario+1;
if l <= 0,
   error('error: m+n too small! Must have m+n-p-1>0')
end

% estimate covariance matrix
v1tilde = (sum(v1')' + sum(v2')' - n*delta_mu)/(m+n);
v2tilde = (sum(v1')' + sum(v2')' + m*delta_mu)/(m+n);
S1 = zeros(p,p);
S2 = zeros(p,p);
for i=1:m,
   v = v1(:,i);
   S1 = S1 + (v-v1tilde)*(v-v1tilde)';
end
for j=1:n,
   v = v2(:,j);
   S2 = S2 + (v-v2tilde)*(v-v2tilde)';
end
Stilde = (S1+S2)/(m+n-1);

% calculate point estimates

gam = sqrt(2*pi/(l+p))/beta(l*.5,.5);
snr = gam*sqrt(delta_mu'*inv(Stilde)*delta_mu);
snr_sq = snr^2;
auc = normcdf(snr/sqrt(2));

% find confidence intervals
% Note: The user may wish to check EXITFLAG and add appropriate error messages, or they may wish to
% modify the tolerances of the fzero function with the OPTIONS argument. 
eta = pi/(beta(l*.5,.5)^2);
a = (l+1)/2;   % alpha 
b0 = [1e-6 1e6];  % search interval
if (omega1 ~= 0 && omega2 ~= 0),
   [beta_L,FVAL,EXITFLAG] = fzero(@(b) gammainc(b/snr_sq,a,'upper')-(1-omega1),b0);
   [beta_U,FVAL,EXITFLAG] = fzero(@(b) gammainc(b/snr_sq,a,'upper')-omega2,b0);  
elseif (omega1 == 0 && omega2 ~= 0),
   beta_L = 0;
   [beta_U,FVAL,EXITFLAG] = fzero(@(b) gammainc(b/snr_sq,a,'upper')-omega2,b0);
elseif (omega1 ~=0 && omega2 == 0),    
   [beta_L,FVAL,EXITFLAG] = fzero(@(b) gammainc(b/snr_sq,a,'upper')-(1-omega1),b0);
   beta_U = Inf;
else
   disp('Warning: Both omega1 and omega2 are zero!')
   beta_L = 0;
   beta_U = Inf;
end
SNR_CI = zeros(1,2);
SNR_CI(1) = sqrt(beta_L/eta);
SNR_CI(2) = sqrt(beta_U/eta);
AUC_CI = normcdf(SNR_CI/sqrt(2));  
ret.m = m;
ret.n = n;
ret.p = p;
ret.alpha1 = alpha1;
ret.alpha2 = alpha2;
ret.SNR = snr;
ret.SNR_CI = SNR_CI;
ret.AUC = auc;
ret.AUC_CI = AUC_CI;

