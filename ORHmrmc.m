% ORHmrmc.m
% Implements the Obuchowski-Rockette-Hillis method for random-reader MRMC 
% analysis (i.e., the OR method with Hillis' degrees of freedom formula) 
% as described in the following references:
%
% Obuchowski NA, Rockette HE, "Hypothesis testing of diagnostic accuracy
% for mutliple readers and multiple tests: an ANOVA approach with dependent
% observations," Communications in Statistics: Simulation and Computation,
% vol. 24, no. 2, pp. 285-308, 1995.
%
% Hillis SL, "A comparison of denominator degrees of freedom methods for
% multiple observer ROC analysis", Statistics in Medicine, vol. 26, pp.
% 596-619, 2007.
%
% Hillis SL, "A marginal-mean ANOVA approach for analyzing multireader
% multicase radiological imaging data," Statistics in Medicine, vol. 33,
% pp. 330-360, 2014.  
%
% The notation generally follows Hillis' paper, where theta is the figure of merit, 
% t is the number of tests (modalities), and r is the number of readers.
%
% Inputs: 
%    thetaHat (t x r matrix of estimated values for theta), 
%    S (tr x tr fixed-reader covariance matrix for the vector TH(:), where
%        TH is the transpose of thetaHat, i.e., TH = thetaHat')
%       NOTE: This matrix can be estimated with fastDeLong.m, EROCcov.m, or
%              binPropCov.m for ROC, LROC/EROC, and MAFC analysis, respectively.  
%    alpha1,alpha2 (lower and upper significance levels for confidence intervals)
%
% Outputs: ret (structure containing point and interval estimates)
%      fields of ret:
%          ret.mth (tx1 vector of reader-averaged estimates for theta_i)
%          ret.CIdiff (t-1 x 2 matrix of confidence intervals for theta_i-theta_1)
%          ret.CIsingle_all (t x 2 matrix of confidence intervals for theta_i estimated using all data)
%          ret.CIsingle (t x 2 matrix of confidence intervals for theta_i estimated using only corresponding data)
% 
% Adam Wunderlich
% 10/24/2014

function ret=ORHmrmc(thetaHat,S,alpha1,alpha2)

[t,r]=size(thetaHat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate Cov1, Cov2, and Cov3 by averaging appropriate entries of S

% create vectors of linear indices in S corresponding to the three covs
I1 = []; I2 = []; I3 = []; % linear indices for Cov1, Cov2, and Cov3, respectively
I2_i = cell(t,1); % cell array of indices for Cov2 specific to each test
for i=1:t,
   for j=1:r,
       for ip=1:t,
           for jp=1:r,
               IND1 = sub2ind([r,t],j,i); % linear index in thetaHat' for (j,i) entry
               IND2 = sub2ind([r,t],jp,ip); % linear index in thetaHat' for (jp,ip) entry
               IND = sub2ind([t*r,t*r],IND1,IND2); % linear index in S for (IND1,IND2) entry
               if ((i~=ip) && (j==jp)),
                  I1 = [I1,IND];
               elseif ((i==ip) && (j~=jp)),
                  I2 = [I2,IND];
                  I2_i{i} = [I2_i{i},IND];
               elseif ((i~=ip) && (j~=jp)),
                  I3 = [I3,IND];
               end
           end
       end
   end
end
% average entries 
Cov1 = mean(S(I1));
Cov2 = mean(S(I2));
Cov3 = mean(S(I3));

Cov2_i = zeros(t,1);
for i=1:t,
    Cov2_i = mean(S(I2_i{i}));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute point and interval estimates

mth_i=mean(thetaHat,2);  % mean of thetaHat over readers (t x 1 vector)
mth_j=mean(thetaHat,1);  % mean of thetaHat over tests (1 x r vector)
mth = mean(mean(thetaHat));  % grand mean of thetaHat (scalar)
MST=(r/(t-1))*sum((mth_i-mth).^2);
MSR=(t/(r-1))*sum((mth_j-mth).^2);
MSTR= (1/((r-1)*(t-1)))*sum(sum((thetaHat-repmat(mth_i,1,r)-repmat(mth_j,t,1)+mth).^2));

% find CIs for differences theta_i-theta_1 using all data
MSdenOR = MSTR + max(r*(Cov2-Cov3),0);
ddf_H = (t-1)*(r-1)*(MSdenOR^2)/(MSTR^2);
d1 = tinv(alpha1,ddf_H)*sqrt((2/r)*MSdenOR);
d2 = tinv(1-alpha2,ddf_H)*sqrt((2/r)*MSdenOR);
CIdiff = zeros(t-1,2);  
for i=2:t,
   mth_diff = mth_i(i) - mth_i(1);
   CIdiff(i-1,1) = mth_diff + d1;
   CIdiff(i-1,2) = mth_diff + d2;
end

% find CI for each theta_i using all data
MSdenOR_single = (MSR+(t-1)*MSTR+t*r*max(Cov2,0))/t;
ddf_denom = (MSR^2)/(r-1) + ((t-1)*MSTR)^2/((t-1)*(r-1));
ddf_H_single = ((t*MSdenOR_single)^2)/ddf_denom;
d1_single = tinv(alpha1,ddf_H_single)*sqrt(MSdenOR_single/r);
d2_single = tinv(1-alpha2,ddf_H_single)*sqrt(MSdenOR_single/r);
CIsingle_all = zeros(t,2);  
CIsingle_all(:,1) = mth_i + d1_single;
CIsingle_all(:,2) = mth_i + d2_single;

% find CI for each theta_i only using data from i'th test
MSR_i = sum((thetaHat - repmat(mth_i,1,r)).^2,2); % t x 1 vector
MSdenOR_single_i = MSR_i + max(r*Cov2_i,0); % t x 1 vector
ddf_denom_i = (MSR_i.^2)/(r-1); % t x 1 vector
ddf_H_single_i = (MSdenOR_single_i.^2)./ddf_denom_i;  % t x 1 vector
d1_i = tinv(alpha1,ddf_H_single_i).*sqrt(MSdenOR_single_i/r);  % t x 1 vector
d2_i = tinv(1-alpha2,ddf_H_single_i).*sqrt(MSdenOR_single_i/r);  % t x 1 vector
CIsingle = zeros(t,2);  
CIsingle(:,1) = mth_i + d1_i;
CIsingle(:,2) = mth_i + d2_i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put results in output structure
ret.mth = mth_i;
ret.CIdiff = CIdiff;
ret.CIsingle_all = CIsingle_all;
ret.CIsingle = CIsingle;



