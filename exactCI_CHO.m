% exactCI_CHO.m
% For a CHO, returns exact 1-(alpha1+alpha2) confidence intervals for SNR and AUC
% as described in 
%
% A. Wunderlich, F. Noo, B. D. Gallas, and M. E. Heilbrun, “Exact confidence 
% intervals for channelized Hotelling observer performance in image quality studies,” 
% IEEE Trans. Med. Imag., vol. 34, no. 2, pp. 453-464, Feb 2015.
%
% Below, m, n, and p denote the number of class 1 images, class 2 images, 
% and channels, respectively.
%
% Inputs: alpha1 (lower significance level), alpha2 (upper significance level) 
% (for conventional two-sided 95% intervals, set alpha1=alpha2=0.025), 
% v1 (p x m matrix of class 1 channel output vectors)
% v2 (p x n matrix of class 2 channel output vectors)
%
% Output: ret (structure containing point and interval estimates of SNR and AUC)
%
% Adam Wunderlich 2/4/2015

function   [ret] = exactCI_CHO(alpha1,alpha2,v1,v2)

[p,m] = size(v1);
[p,n] = size(v2);
v1bar = mean(v1')';
v2bar = mean(v2')';
S1 = zeros(p,p);
S2 = zeros(p,p);
for i=1:m,
   v = v1(:,i);
   S1 = S1 + (v-v1bar)*(v-v1bar)';
end
for j=1:n,
   v = v2(:,j);
   S2 = S2 + (v-v2bar)*(v-v2bar)';
end
S = (S1+S2)/(m+n-2);
delta_vbar = v2bar-v1bar;
gam = (m+n-p-3)/(m+n-2);
snr_sq = gam*delta_vbar'*inv(S)*delta_vbar;
SNR = sqrt(snr_sq);
AUC = normcdf(SNR/sqrt(2));

% find confidence interval
x = snr_sq*(m+n-p-1)*(m*n)/(p*(m+n-2)*(m+n)*gam);
nu1 = p;
nu2 = m+n-p-1;
d0 = [1e-6 1e5];  % search interval
if ncfcdf(x,nu1,nu2,1e-6) < 1-alpha1,
   delta_L = 0;
else
   [delta_L,FVAL,EXITFLAG] = fzero(@(d) ncfcdf(x,nu1,nu2,d)-(1-alpha1),d0);
end
if ncfcdf(x,nu1,nu2,1e-6) < alpha2,
   delta_U = 0;
else
   [delta_U,FVAL,EXITFLAG] = fzero(@(d) ncfcdf(x,nu1,nu2,d)-alpha2,d0);
end
SNR_CI = zeros(1,2);
SNR_CI(1) = sqrt(delta_L*(m+n)/(m*n));
SNR_CI(2) = sqrt(delta_U*(m+n)/(m*n));
AUC_CI = normcdf(SNR_CI/sqrt(2));
ret.m = m;
ret.n = n;
ret.p = p;
ret.alpha1 = alpha1;
ret.alpha2 = alpha2;
ret.SNR = SNR;
ret.SNR_CI = SNR_CI;
ret.AUC = AUC;
ret.AUC_CI = AUC_CI;