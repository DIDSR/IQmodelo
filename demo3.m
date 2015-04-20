% demo3.m
% Demonstrate LROC analysis for a single imaging scenario.
%
% Adam Wunderlich
% 7/7/2014

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
Nx=48;  % x-dimension
Ny=48;  % y-dimension
params.Nx = Nx;   
params.Ny = Nx;   
sigScale = 2.5;   % signal scale parameter (Gaussian stdev)
params.sigScale = sigScale;
params.A = 10;   % signal amplitude
params.Rbg1 = 10;  % 2*amplitude of noise component 1
params.Rbg2 = 35;  % 2*amplitude of noise component 1

x=-Nx/2:(Nx/2-1);  % x coordinates 
y=-Ny/2:(Ny/2-1);  % y coordinates 

disp('----------------------- demo 3 ----------------------------------');

% get signal-absent images
disp(['computing lesion-absent images (m=',num2str(m),') ...'])
disp('')
g1 = zeros(Nx,Ny,m);
for i=1:m,
   params.sd = i;   % seed for random number generator
   params.sp = false;
   params.zc = [0;0];
   [g,s] = create_images1(params);
   g1(:,:,i) = g;  
end

% get signal-present images
disp(['computing lesion-present images (n=',num2str(n),') ...'])
disp('')
g2 = zeros(Nx,Ny,n);
zc = zeros(2,n);
for j=1:n,
   params.sd = 1e3+j;   % seed for random number generator for create_images1
   params.sp = true;
   rand('state',m+j);  % set state of random number generator for signal location
   zc(:,j) = (rand(2,1)-.5)*(Nx-sigScale*3);  % lesion center
   params.zc = zc(:,j);
   [g,s] = create_images1(params);
   g2(:,:,j) = g;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixed-template scanning linear observer
disp(['scanning linear observer:'])
disp('')

t1 = zeros(Nx*Ny,m);  % class 1 ratings 
t2 = zeros(Nx*Ny,n);  % class 2 ratings 

% define template (uniform disk of radius sigScale)
for k=1:Nx,
   for l=1:Ny,
      w(k,l) = double(norm([x(k);y(l)]) <= sigScale);
   end
end

b = zeros(Nx,Ny);  % known background
% if background is unknown, can estimate as mean of signal-absent images

% convolve template with each image after background subtraction
for i=1:m,
    t = conv2(w,g1(:,:,i)-b,'same');
    t1(:,i) = t(:);  % convert 2-D image to 1-D column vector
end
for j=1:n,
    t = conv2(w,g2(:,:,j)-b,'same');
    t2(:,j) = t(:);  % convert 2-D image to 1-D column vector
end
[t1max,imax] = max(t1);
[t2max,jmax] = max(t2);

% check localization
u = zeros(1,n); % localization outcomes
for j=1:n,  % loop over lesion-present images
   % form correct localization template
   for k=1:Nx,
      for l=1:Ny,
         L(k,l) = double(norm([x(k);y(l)] - zc(:,j)) <= sigScale);
      end
   end
   LI = find(L>0);  % indices of correct localization template
   ITC = intersect(LI,jmax(j));
   if ~isempty(ITC),
       u(j) = 1;
   end
end

X = t1max;
Y = t2max;
U = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute figures of merit and confidence intervals
disp(' ')
disp(['Displaying point estimates and two-sided ',num2str((1-alpha)*100,2),'% confidence intervals...'])
disp(' ')
disp('Probability of correct localization:')
[PCL,PCL_CI] = binProp_CI(alpha1,alpha2,U')
disp('Area under LROC curve:')
[AL,AL_CI] = npAEROC_CI(alpha1,alpha2,X,Y,U)
