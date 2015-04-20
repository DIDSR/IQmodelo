% diffCI_ktkm.m
% For a linear observer defined by a known (fixed) template, and known 
% difference of class means, returns approximate 1-(alpha1+alpha2) 
% confidence intervals for SNR and AUC differences as described in
%
% A. Wunderlich and F. Noo, "On Efficient Assessment of Image Quality Metrics Based on Linear Model 
% Observers," IEEE Transactions on Nuclear Science, vol. 59, no. 3, pp. 568-578, June 2012.
%
% Inputs: alpha1 (lower significance level), alpha2 (upper significance level) 
% deltaA (known difference of class means for scenario A), 
% deltaB (known difference of class means for scenario B),
% xA (1 x m vector of class 1 ratings for scenario A), 
% yA (1 x n vector of class 2 ratings for scenario A),
% xB (1 x m vector of class 1 ratings for scenario B), 
% yB (1 x n vector of class 2 ratings for scenario B) 
%
% Outputs: ret (structure containing confidence intervals for SNR_B-SNR_A
%               and AUC_B-AUC_A)
%
% Adam Wunderlich 6/4/2014


function  [ret] = diffCI_ktkm(alpha1,alpha2,deltaA,deltaB,xA,yA,xB,yB)

   if (length(xA) ~= length(xB))
       error('error: x vectors are not same length!')
       return
   end
   if (length(yA) ~= length(yB))
       error('error: y vectors are not same length!')
       return
   end
   m = length(xA);
   n = length(yA);
   if m+n < 2,
      error('error: m+n too small!')
      return
   end
   q = m+n-1;
   gam = sqrt(2*pi/q)/beta((q-1)/2,0.5);
   eta = q*(gam^2)/2;
   
   % scenario A
   xtildeA = (sum(xA)+sum(yA)-n*deltaA)/(m+n);
   ytildeA = (sum(xA)+sum(yA)+m*deltaA)/(m+n);
   z = [(xA-xtildeA*ones(1,m)).^2,(yA-ytildeA*ones(1,n)).^2];
   stilde_sq = sum(z)/q;
   sA = sqrt(stilde_sq);
   snrA = gam*deltaA/sA;
   aucA = normcdf(snrA/sqrt(2));
   
   % scenario B
   xtildeB = (sum(xB)+sum(yB)-n*deltaB)/(m+n);
   ytildeB = (sum(xB)+sum(yB)+m*deltaB)/(m+n);
   z = [(xB-xtildeB*ones(1,m)).^2,(yB-ytildeB*ones(1,n)).^2];
   stilde_sq = sum(z)/q;
   sB = sqrt(stilde_sq);
   snrB = gam*deltaB/sB;
   aucB = normcdf(snrB/sqrt(2));

   % estimate covariance of ratings
   z = [(xA-xtildeA*ones(1,m)).*(xB-xtildeB*ones(1,m)),(yA-ytildeA*ones(1,n)).*(yB-ytildeB*ones(1,n))];
   sAB = sum(z)/q;
   r = sAB/(sA*sB);  % estimate of correlation coefficient

   % evaluate the Gauss hypergeometric function 2F1(1/2,1/2,q/2,z) 
   % using a 50-term series approximation
   nterm = 50;  
   z = r^2;
   u=zeros(nterm);
   u(1)=1;  % k=0 term
   for j=1:(nterm-1),
      k=j-1;
      u(j+1)=u(j)*(z/(k+1))*((1/2+k)^2)/(q/2+k); % get (k+1)th term from kth term
   end
   % sum terms from smallest to largest to get best numerical accuracy
   F12=0;
   for j=nterm:-1:1,
      F12=F12+u(j);
   end
   
   % estimate variances and covariance of SNRs using exact expressions
   varSNR_A = (2*eta/(q-2) - 1)*snrA^2;
   varSNR_B = (2*eta/(q-2) - 1)*snrB^2;
   covSNR = snrA*snrB*(F12 - 1);
 
   % estimate variance and covariance of auc estimates
   varAUC_A = 0.5*varSNR_A*normpdf(snrA/sqrt(2))^2;
   varAUC_B = 0.5*varSNR_B*normpdf(snrB/sqrt(2))^2;
   covAUC = 0.5*normpdf(snrA/sqrt(2))*normpdf(snrB/sqrt(2))*covSNR;

   % construct Wald confidence intervals for SNR_B-SNR_A and AUC_B-AUC_A
   SNRdiff = snrB-snrA;
   var_SNRdiff  = varSNR_A + varSNR_B - 2*covSNR;  % variance of estimated difference
   std_SNRdiff = sqrt(var_SNRdiff);
   AUCdiff = aucB-aucA;
   var_AUCdiff  = varAUC_A + varAUC_B - 2*covAUC;  % variance of estimated difference
   std_AUCdiff = sqrt(var_AUCdiff);
   if (alpha1~=0 && alpha2~=0)
      SNRdiff_CI(1) = SNRdiff + norminv(alpha1)*std_SNRdiff;
      SNRdiff_CI(2) = SNRdiff + norminv(1-alpha2)*std_SNRdiff;
      AUCdiff_CI(1) = AUCdiff + norminv(alpha1)*std_AUCdiff;
      AUCdiff_CI(2) = AUCdiff + norminv(1-alpha2)*std_AUCdiff;
   elseif (alpha1==0 && alpha2~=0)
      SNRdiff_CI(1) = -inf;
      SNRdiff_CI(2) = SNRdiff + norminv(1-alpha2)*std_SNRdiff;
      AUCdiff_CI(1) = -1;
      AUCdiff_CI(2) = AUCdiff + norminv(1-alpha2)*std_AUCdiff;
   elseif (alpha1~=0 && alpha2==0), 
      SNRdiff_CI(1) = SNRdiff + norminv(alpha1)*std_SNRdiff;
      SNRdiff_CI(2) = inf;
      AUCdiff_CI(1) = AUCdiff + norminv(alpha1)*std_AUCdiff;
      AUCdiff_CI(2) = 1;
   end   
      
   ret.m=m;
   ret.n=n;
   ret.alpha1=alpha1;
   ret.alpha2=alpha2;
   ret.aucA = aucA;
   ret.aucB = aucB;
   ret.SNRdiff_CI = SNRdiff_CI;
   ret.AUCdiff_CI = AUCdiff_CI;
