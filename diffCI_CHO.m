% diffCI_CHO.m
% For CHOs, returns conservative 1-(alpha1+alpha2) confidence intervals 
% for a difference of SNR or AUC values.  This method uses exact intervals 
% for each imaging scenario obtained with exactCI_CHO.m together with the 
% Bonferroni inequality.   
%
% Inputs: alpha1 (lower significance level), alpha2 (upper significance level) 
% vA1 (pA x mA matrix of class 1 channel output vectors for scenario A)
% vA2 (pA x nA matrix of class 2 channel output vectors for scenario A)
% vB1 (pB x mB matrix of class 1 channel output vectors for scenario B)
% vB2 (pB x nB matrix of class 2 channel output vectors for scenario B)
%
% Outputs: ret (structure containing confidence intervals for SNR_B-SNR_A
%               and AUC_B-AUC_A)
%
% Adam Wunderlich 6/4/2014


function  [ret] = diffCI_CHO(alpha1,alpha2,vA1,vA2,vB1,vB2)
  
   % find lower bound on SNR and AUC differences
   if alpha1==0,
      SNRdiff_CI(1) = -inf;
      AUCdiff_CI(1) = -1;
   else     
      [retA] = exactCI_CHO(0,alpha1/2,vA1,vA2);
      SNR_CI_A  = retA.SNR_CI;
      AUC_CI_A  = retA.AUC_CI;
      [retB] = exactCI_CHO(alpha1/2,0,vB1,vB2);
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
      [retA] = exactCI_CHO(alpha2/2,0,vA1,vA2);
      SNR_CI_A  = retA.SNR_CI;
      AUC_CI_A  = retA.AUC_CI;
      [retB] = exactCI_CHO(0,alpha2/2,vB1,vB2);
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
