% demo2a.m
% Demonstrate AUC confidence interval estimators for the case of two 
% imaging scenarios (called A and B), computing a two-sided confidence
% interval for the difference AUC_B - AUC_A.  
%
% Adam Wunderlich
% 5/22/2014

clc
clear all
close all
fclose('all');

alpha = .05;  % significance level
% two-sided intervals
alpha1=alpha/2;  % lower significance level
alpha2=alpha/2;  % upper significance level 

m=150;  % number of lesion-absent images
n=150;  % number of lesion-present images

% image parameters
Nx=96;   % x-dimension
Ny=96;   % y-dimension
params.Nx = Nx;
params.Ny = Nx;
params.zc = [0;0];  % coordinates of signal center
params.sigScale = 2.5; % signal scale parameter (Gaussian stdev)
params.A = 10;  % signal amplitude
params.Rbg1 = 25;  % 2*amplitude of noise component 1
params.Rbg2 = 75; % 2*amplitude of noise component 2

x=-Nx/2:(Nx/2-1);  % x coordinates 
y=-Ny/2:(Ny/2-1);  % y coordinates
dx = .05;  % pixel width (cm) used to define channels


disp('----------------------- demo 2a ----------------------------------');

% get signal-absent images
disp(['computing lesion-absent images (m=',num2str(m),') ...'])
disp('')
gA1 = zeros(Nx*Ny,m);  % modality A
gB1 = zeros(Nx*Ny,m);  % modality B
for i=1:m,
   params.sd = i;   % seed for random number generator
   params.sp = false;
   [gA,gB,s] = create_images2(params);
   gA1(:,i) = gA(:);
   gB1(:,i) = gB(:);
end

% get signal-present images
disp(['computing lesion-present images (n=',num2str(n),') ...'])
disp('')
gA2 = zeros(Nx*Ny,n);  % modality A
gB2 = zeros(Nx*Ny,n);  % modality B
for j=1:n,
   params.sd = m+j;    % seed for random number generator
   params.sp = true;
   [gA,gB,s] = create_images2(params);
   gA2(:,j) = gA(:);
   gB2(:,j) = gB(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(['Displaying two-sided ',num2str((1-alpha)*100,2),'% confidence intervals for AUC difference:'])
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixed template linear observer

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('fixed template linear observer')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% define template (uniform disk of radius sigScale centered at zc)
for k=1:Nx,
   for l=1:Ny,
      wm(k,l) = double(norm([x(k);y(l)] - params.zc) <= params.sigScale);
   end
end
w = wm(:);  % column vector

% compute ratings
xA = w'*gA1; % modality A, class 1
yA = w'*gA2; % modality A, class 2
xB = w'*gB1; % modality B, class 1
yB = w'*gB2; % modality B, class 2

% compute 2AFC outcomes
tA = double(yA>xA);
tB = double(yB>xB);

% 2AFC estimates 
disp('2AFC estimates:')
X = [tA',tB'];
% interval for difference 
[p,pdiff_CI] = binProp_CI(alpha1,alpha2,X)


disp('nonparametric (Mann-Whitney) estimates:')
% nonparametric estimates based on ratings
X = [xA;xB];
Y = [yA;yB];
% interval for difference 
[AUC,AUCdiff_CI] = npAUC_CI(alpha1,alpha2,X,Y)

% parametric results
% find 1-alpha confidence interval for difference using Bonferroni
% correction
disp('parametric estimates (using Bonferroni for CI):')

[ret] = diffCI_kt(alpha1,alpha2,xA,yA,xB,yB);
AUC(1) = ret.aucA;
AUC(2) = ret.aucB;
AUC
AUCdiff_CI = ret.AUCdiff_CI

disp('parametric estimates using known difference of class means:')
delta = w'*s(:); % same for each scenario in this example
deltaA = delta;
deltaB = delta;
[ret] = diffCI_ktkm(alpha1,alpha2,deltaA,deltaB,xA,yA,xB,yB);
AUC(1) = ret.aucA;
AUC(2) = ret.aucB;
AUC
AUCdiff_CI = ret.AUCdiff_CI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHO

disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Channelized Hotelling Observer (CHO) with 3 DOG channels')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% compute channel matrix
[U] = make_channels(2,Nx,Ny,dx); % 3 DOG channels
vA1 = U'*gA1;
vA2 = U'*gA2;
vB1 = U'*gB1;
vB2 = U'*gB2;

disp('parametric estimates (using Bonferroni for CI):')
[ret] = diffCI_CHO(alpha1,alpha2,vA1,vA2,vB1,vB2);
AUC(1) = ret.aucA;
AUC(2) = ret.aucB;
AUC
AUCdiff_CI = ret.AUCdiff_CI

disp('parametric estimates using known difference of class means and Bonferroni:')
delta_mu = U'*s(:);  % same for each modality in this example
delta_muA = delta_mu;
delta_muB = delta_mu;
[ret] = diffCI_CHO_km(alpha1,alpha2,delta_muA,delta_muB,vA1,vA2,vB1,vB2);
AUC(1) = ret.aucA;
AUC(2) = ret.aucB;
AUC
AUCdiff_CI = ret.AUCdiff_CI

disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Channelized Hotelling Observer (CHO) with 40 Gabor channels')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


% compute channel matrix
[U] = make_channels(1,Nx,Ny,dx); % 40 Gabor channels
vA1 = U'*gA1;
vA2 = U'*gA2;
vB1 = U'*gB1;
vB2 = U'*gB2;

disp('parametric estimates (using Bonferroni for CI):')
[ret] = diffCI_CHO(alpha1,alpha2,vA1,vA2,vB1,vB2);
AUC(1) = ret.aucA;
AUC(2) = ret.aucB;
AUC
AUCdiff_CI = ret.AUCdiff_CI

disp('parametric estimates using known difference of class means and Bonferroni:')
delta_mu = U'*s(:);  % same for each modality in this example
delta_muA = delta_mu;
delta_muB = delta_mu;
[ret] = diffCI_CHO_km(alpha1,alpha2,delta_muA,delta_muB,vA1,vA2,vB1,vB2);
AUC(1) = ret.aucA;
AUC(2) = ret.aucB;
AUC
AUCdiff_CI = ret.AUCdiff_CI
