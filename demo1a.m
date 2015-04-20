% demo1a.m
% Demonstrate AUC point estimators and two-sided confidence intervals for a 
% single imaging scenario.
%
% Adam Wunderlich
% 5/27/2014

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
b = (params.Rbg1+params.Rbg2)/2;  % parameter for example image grayscale 

disp('----------------------- demo 1a ----------------------------------');

% get signal-absent images
disp(['computing lesion-absent images (m=',num2str(m),') ...'])
disp('')
g1 = zeros(Nx*Ny,m);
for i=1:m,
   params.sd = i;   % seed for random number generator
   params.sp = false;
   [g,s] = create_images1(params);
   g1(:,i) = g(:);
   if i==1,
      figure
      subplot(2,2,1)
      imagesc(g,[0 b]); colormap(gray); axis image; axis xy; axis off; colorbar
      title('example signal-absent image')
   end
  
end

% get signal-present images
disp(['computing lesion-present images (n=',num2str(n),') ...'])
disp('')
g2 = zeros(Nx*Ny,n);
for j=1:n,
   params.sd = m+j;   % seed for random number generator
   params.sp = true;
   [g,s] = create_images1(params);
   g2(:,j) = g(:);
   if j==1,
      subplot(2,2,2)
      imagesc(g,[0 b]); colormap(gray); axis image; axis xy; axis off; colorbar
      title('example signal-present image')
      subplot(2,2,3)
      imagesc(s,[0 b/2]); colormap(gray); axis image; axis xy; axis off; colorbar
      title('signal')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(['Displaying two-sided ',num2str((1-alpha)*100,2),'% confidence intervals for AUC:'])
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixed template linear observer

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('fixed template linear observer:')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% define template (uniform disk of radius sigScale centered at zc)
for k=1:Nx,
   for l=1:Ny,
      wm(k,l) = double(norm([x(k);y(l)] - params.zc) <= params.sigScale);
   end
end
w = wm(:);  % column vector

% compute ratings
x = w'*g1; % class 1
y = w'*g2; % class 2

% compute 2AFC outcomes
t = double(y>x);

% 2AFC estimates 
disp('2AFC estimates:')
[p,p_CI] = binProp_CI(alpha1,alpha2,t')

disp('nonparametric (Mann-Whitney) estimates:')
% nonparametric estimates based on ratings
[AUC,AUC_CI] = npAUC_CI(alpha1,alpha2,x,y)

% parametric results
disp('parametric estimates:')
[ret] = exactCI_kt(alpha1,alpha2,x,y);
AUC = ret.AUC
AUC_CI  = ret.AUC_CI

delta = w'*s(:);
disp('parametric estimates using known difference of class means:')
[ret] = exactCI_ktkm(alpha1,alpha2,delta,x,y);
AUC = ret.AUC
AUC_CI  = ret.AUC_CI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHO

disp('')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Channelized Hotelling Observer (CHO) with 3 DOG channels')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% compute channel matrix
[U] = make_channels(2,Nx,Ny,dx); % 3 DOG channels
v1 = U'*g1;
v2 = U'*g2;

disp('parametric estimates:')
[ret] = exactCI_CHO(alpha1,alpha2,v1,v2);
AUC = ret.AUC
AUC_CI  = ret.AUC_CI

disp('parametric estimates using known difference of class means:')
delta_mu = U'*s(:);
[ret] = exactCI_CHO_km(alpha1,alpha2,delta_mu,v1,v2);
AUC = ret.AUC
AUC_CI  = ret.AUC_CI

disp('')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Channelized Hotelling Observer (CHO) with 40 Gabor channels')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


% compute channel matrix
[U] = make_channels(1,Nx,Ny,dx); % 40 Gabor channels
v1 = U'*g1;
v2 = U'*g2;

disp('parametric estimates:')
[ret] = exactCI_CHO(alpha1,alpha2,v1,v2);
AUC = ret.AUC
AUC_CI  = ret.AUC_CI

disp('parametric estimates using known difference of class means:')
delta_mu = U'*s(:);
[ret] = exactCI_CHO_km(alpha1,alpha2,delta_mu,v1,v2);
AUC = ret.AUC
AUC_CI  = ret.AUC_CI



   