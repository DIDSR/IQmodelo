% make_channels.m
% Construct (normalized) channel matrix for a 2-D image.
% For details and references see
% A. Wunderlich and F. Noo, "Evaluation of the Impact of Tube Current 
% Modulation on Lesion Detectability using Model Observers," Proc. Intl. 
% Conf. IEEE Eng Med. Biol. Soc., pp. 2705-2708, Aug. 2008.
%
% Inputs: chanType (specifies type of channels), Nx (number of x pixels), 
%         Ny (number of y pixels), dx (pixel size in cm)
%  Denoting the number of channels by p, the channels are specified as 
%  follows:
%  chanType = 1: Gabor (p=40), chanType = 2: DOG (p=3), 
%  chanType = 3: dense DOG (p=10), chanType = 4: square (p=4)
%
% Output: U (Nx*Ny x p channel matrix)
%
% Adam Wunderlich 
% 5/28/2014

function [U] = make_channels(chanType,Nx,Ny,dx)

dy = dx; % y stepsize
x_off = Nx/2;
y_off = Ny/2;
i = 0:(Nx-1);
j = 0:(Ny-1);
x = (i-x_off)*dx;
y = (j-y_off)*dy;
X = [0;0];  % origin for channel grid

switch chanType
  case 1,  % Gabor functions  
     pb = [1/64,1/32,1/16,1/8,1/4];  % edges of passbands (cycles/pixel)
     wf = diff(pb); % bandwidths
     ws = (4/pi*log(2)./wf)*dx;  % spatial channel widths (in cm)
     fc = (pb(1:4)+wf/2)/dx;    % center frequencies (in cycles/cm)
     dtheta = 1/5;
     theta = (0:dtheta:(1-dtheta))*pi;  % orientations
     C = zeros(Nx,Ny);
     U = zeros(Nx*Ny,40);   % channel matrix
     x2 = x - X(1);
     y2 = y - X(2);
     beta = [0;pi/2];
     count = 0;
     for p=1:2, % loop over phases
        for l=1:4, % loop over frequencies
           for k=1:5,  % loop over angles
              count = count +1;
              for i=1:Nx, % x-loop
                 for j=1:Ny, % y-loop
                    C(i,j) = exp(-4*log(2)*(x2(i)^2 + y2(j)^2)/(ws(l)^2))* ...
                       cos(2*pi*fc(l)*(x2(i).*cos(theta(k)) + ...
                       y2(j).*sin(theta(k))) + beta(p)); 
                 end
              end
              C = C./norm(C,'fro');
              U(:,count) = C(1:(Nx*Ny))';       
           end
        end
     end
  case 2,  % sparse DOG    
     Q = 2;
     alpha = 2;
     sigma_0 = 0.015;
     C = zeros(Nx,Ny);
     U = zeros(Nx*Ny,3);
     x2 = (x - X(1))/dx;
     y2 = (y - X(2))/dy;
     for c=1:3,
        for i=1:Nx, % x-loop
           for j=1:Ny, % y-loop
              sigma = (alpha^c)*sigma_0;
              rsq = x2(i)^2 + y2(j)^2;
              C(i,j) =  2*pi*Q^2*sigma^2*exp(-(pi^2)*2*Q^2*sigma^2*rsq) ...
                    - 2*pi*sigma^2*exp(-(pi^2)*2*sigma^2*rsq); 
           end
        end
        C = C./norm(C,'fro');
        U(:,c) = C(1:(Nx*Ny))';
     end
  case 3,  % dense DOG 
     sigma_0 = .005; 
     Q = 1.67;
     alpha = 1.4;
     C = zeros(Nx,Ny);
     U = zeros(Nx*Ny,10);
     x2 = (x - X(1))/dx;
     y2 = (y - X(2))/dy;
     for c=1:10,
        for i=1:Nx, % x-loop
           for j=1:Ny, % y-loop
              sigma = (alpha^c)*sigma_0;
              rsq = x2(i)^2 + y2(j)^2;
              C(i,j) =  2*pi*Q^2*sigma^2*exp(-(pi^2)*2*Q^2*sigma^2*rsq) ...
                    - 2*pi*sigma^2*exp(-(pi^2)*2*sigma^2*rsq); 
           end
        end
        C = C./norm(C,'fro');
        U(:,c) = C(1:(Nx*Ny))';
     end
  case 4,   % square channels
     pb = [1/64,1/32,1/16,1/8,1/4]/dx;  % edges of passbands (cycles/cm)
     C = zeros(Nx,Ny);
     U = zeros(Nx*Ny,4);
     x2 = x - X(1);
     y2 = y - X(2);
     for c=1:4,
        for i=1:Nx, % x-loop
           for j=1:Ny, % y-loop
              r = sqrt(x2(i).^2 + y2(j)^2);
              if r > 0,
                C(i,j) = (pb(c+1)*besselj(1,2*pi*pb(c+1)*r) - ...
                          pb(c)*besselj(1,2*pi*pb(c)*r))/r;
              else
                C(i,j) = pi*(pb(c+1)^2 - pb(c)^2);
              end
           end
        end
        C = C./norm(C,'fro');
        U(:,c) = C(1:(Nx*Ny))';
     end
end  % end switch
