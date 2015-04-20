% create_images2.m
% Generate correlated random images for two imaging scenarios applied to the same
% case, where both scenarios have the same overall task-performance.
%
% Inputs:
% The input is a structure with the following fields:
% sp (boolean variable specifying if signal-present images are desired)
% Nx,Ny (number of pixels in x and y dimensions, respectively)
% Rbg1, Rbg2 (range (twice amplitude) of background noise components)
% zc (2x1 vector of signal center coordinates in the range -Nx/2 to Nx/2)
% sigScale (standard deviation of guassian signal, in pixels)
% A (signal amplitude), sd (seed for random number generator)
%
% Outputs: gA,gB (random Nx x Ny images for modalities A and B)
% signal (Nx x Ny image of signal)
%          
% Adam Wunderlich
% 5/23/2014

function [gA,gB,s] = create_images2(params)

sp = params.sp;
Nx = params.Nx;
Ny = params.Ny;
zc = params.zc;
sigScale = params.sigScale;
A = params.A;
Rbg1 = params.Rbg1;
Rbg2 = params.Rbg2;
sd = params.sd;

% set state of random number generator to get repeatable results
randn('state',sd); 

% generate random background common to both modalities
rb = randn(Nx,Ny);  % uncorrelated Gaussian noise with unit variance
% apply cone-filter to color noise
Rb = fftshift(fft2(rb));
beta=2;
nux=(-Nx/2:Nx/2-1)/Nx;
nuy=(-Ny/2:Ny/2-1)/Ny;
[nuy,nux]=meshgrid(nuy,nux);
conefilter=((nux.^2+nuy.^2).^(-beta/4));
j=find(abs(nux) < 3.0/Nx & abs(nuy) < 3.0/Ny);
conefilter(j)=0;
rb=ifft2(fftshift(Rb.*conefilter));
minrb=min(rb(:));maxrb=max(rb(:));
rb1=((rb-minrb)/(maxrb-minrb)-1/2)*Rbg1; 

% compute signal
x=-Nx/2:(Nx/2-1);
y=-Ny/2:(Ny/2-1);
L = length(x);
Sinv = eye(2,2)/(2*sigScale^2);
for i=1:L,
   for j=1:L,
      z = [x(i);y(j)] - zc;
      s(i,j) = A*exp(-z'*Sinv*z);
   end
end   

for modality=1:2,    
   % generate random "noise" background that is different for each modality
   rb = randn(Nx,Ny);  % uncorrelated Gaussian noise with unit variance
   % apply cone-filter to color noise
   Rb = fftshift(fft2(rb));
   beta=2;
   nux=(-Nx/2:Nx/2-1)/Nx;
   nuy=(-Ny/2:Ny/2-1)/Ny;
   [nuy,nux]=meshgrid(nuy,nux);
   conefilter=((nux.^2+nuy.^2).^(-beta/4));
   j=find(abs(nux) < 3.0/Nx & abs(nuy) < 3.0/Ny);
   conefilter(j)=0;
   rb=ifft2(fftshift(Rb.*conefilter));
   minrb=min(rb(:));maxrb=max(rb(:));
   rb2=((rb-minrb)/(maxrb-minrb)-1/2)*Rbg2; 
   
   if sp,  % generate signal-present image 
      g = rb1+rb2+s;
   else  % signal-absent image
      g = rb1+rb2;
   end
   if modality==1,
       gA = g;
   else
       gB = g;
   end
end

