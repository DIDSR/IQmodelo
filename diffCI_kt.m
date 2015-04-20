% diffCI_kt.m
% For a linear observer defined by a known (fixed) template, returns 
% conservative 1-(alpha1+alpha2) confidence intervals for a difference of 
% SNR or AUC values.  This method uses exact intervals for each scenario 
% obtained with exactCI_kt.m together with the Bonferroni inequality.   
%
% Inputs: alpha1 (lower significance level), alpha2 (upper significance level) 
% xA (1 x mA vector of class 1 ratings for scenario A) 
% yA (1 x nA vector of class 2 ratings for scenario A)
% xB (1 x mB vector of class 1 ratings for scenario B) 
% yB (1 x nB vector of class 2 ratings for scenario B) 
%
% Outputs: ret (structure containing confidence intervals for SNR_B-SNR_A
%               and AUC_B-AUC_A)
%
% Adam Wunderlich 6/4/2014


function  [ret] = diffCI_kt(alpha1,alpha2,xA,yA,xB,yB)
  
   % find lower bound on SNR and AUC differences
   if alpha1==0,
      SNRdiff_CI(1) = -inf;
      AUCdiff_CI(1) = -1;
   else
      [retA] = exactCI_kt(0,alpha1/2,xA,yA);
      SNR_CI_A  = retA.SNR_CI;
      AUC_CI_A  = retA.AUC_CI;
      [retB] = exactCI_kt(alpha1/2,0,xB,yB);
      SNR_CI_B  = retB.SNR_CI;
      AUC_CI_B  = retB.AUC_CI;
      SNRdiff_CI(1) = SNR_CI_B(1)-SNR_CI_A(2);
      AUCdiff_CI(1) = AUC_CI_B(1)-AUC_CI_A(2);
   end
   % find upper bound on AUC difference
   if alpha2==0,
      SNRdiff_CI(2) = inf;
      AUCdiff_CI(2) = 1;
   else
      [retA] = exactCI_kt(alpha2/2,0,xA,yA);
      SNR_CI_A  = retA.SNR_CI;
      AUC_CI_A  = retA.AUC_CI;
      [retB] = exactCI_kt(0,alpha2/2,xB,yB);
      SNR_CI_B  = retB.SNR_CI;
      AUC_CI_B  = retB.AUC_CI;
      SNRdiff_CI(2) = SNR_CI_B(2)-SNR_CI_A(1);
      AUCdiff_CI(2) = AUC_CI_B(2)-AUC_CI_A(1);
   end
 
   ret.alpha1=alpha1;
   ret.alpha2=alpha2;
   ret.aucA = retA.AUC;
   ret.aucB = retB.AUC;
   ret.SNRdiff_CI = SNRdiff_CI;
   ret.AUCdiff_CI = AUCdiff_CI;
