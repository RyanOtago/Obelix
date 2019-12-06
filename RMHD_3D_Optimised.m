function RMHD_3D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                               %%%
%%%   Draft 3D RMHD Solver        %%%
%%%                               %%%
%%%   Uses Elsasser Formulation   %%%
%%%   for an Alfven Wave in a     %%%
%%%   3D Torus                    %%%
%%%                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%% Paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

va = 1;        % Alven velocity
nu = 1e-3;

LX = 2*pi;     % Box-size (x-direction)
LY = 2*pi;     % Box-size (y-direction)
LZ = 2*pi;     % Box-size (z-direction)     %%% !!! Should this scale differently to LX and LY? !!! %%%

NX = 64;      % Resolution in x
NY = 64;      % Resolution in y
NZ = 64;      % Resolution in z
N = NX*NY*NZ;

dt = 1e-3;     % Time Step              %%% !!! Think about CFL conditions !!! %%% (Will eventually implement variable time step)
TF = 1.0;    % Final Time
TSCREEN = 100; % Sreen Update Interval Count (NOTE: plotting is usually slow)
time = [dt:dt:TF];

E_plus = zeros(1,length(time));
E_minus = zeros(1,length(time));

I=sqrt(-1);
dx = LX/NX;
dy = LY/NY;
dz = LZ/NZ;
dV = dx*dy*dz;

grid_int = dV/N;
t=0.0;

%% Initialise wavevector grid %%%%%%%%%%%%%%%%%%

kx = (2*I*pi/LX)*[0:((NX/2)-1)  -(NX/2):-1];    % [0, 1, ..., NX/2-1, -NX/2, -NX/2+1, ..., -1]    % This is a formatting convention
ky = (2*I*pi/LY)*[0:((NY/2)-1)  -(NY/2):-1];
kz = (2*I*pi/LZ)*[0:((NZ/2)-1)  -(NZ/2):-1];
[KX, KY, KZ] = ndgrid(kx, ky, kz);

dealias = abs(KX)<(1/3)*NX & abs(KY)<(1/3)*NY & abs(KZ)<(1/3)*NZ;            % Cutting of frequencies using the 2/3 rule   (2/3)*(N/2)

k2_perp = KX.^2 + KY.^2;      % (Perpendicular) Laplacian in Fourier space
k2_poisson = k2_perp;    
k2_poisson(1,1,:) = 1;          % Fixed Laplacian in F.S. for Poisson's equation   % Because first entry was 0
k2 = k2_perp + KZ.^2;         % Laplacian in F.S.
exp_correct = exp(dt*nu*k2);  % !!! Should this be k2 or k2_perp ??? (Seems like it should be k2 as will have z-dependence)

%% Initial Condition %%%%%%%%%%%%%%%%%%%%%%%%%%%

[i,j,k]=ndgrid((0:NX-1)*dx,(0:NY-1)*dy,(0:NZ-1)*dz);
% Lap_z_plus = k2_perp.*fftn(0.01*sin(2*pi*(i/LX + j/LY + k/LZ)/3));
% Lap_z_minus = k2_perp.*fftn(0.02*cos(2*pi*(i/LX + j/LY + k/LZ)/3));

Lap_z_plus = k2_perp.*fftn(1);
Lap_z_minus = k2_perp.*fftn(0);
tic
%% Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0;
n=1;
for i = [1:length(time)]
    k=k+1;

    z_plus = Lap_z_plus./k2_poisson;        
    z_minus = Lap_z_minus./k2_poisson;

    % Computes Poisson Brackets (RHS of Schekochihin-09 Eq (21))
    zp_x = real(ifft2(KX.*z_plus));
    zm_x = real(ifft2(KX.*z_minus));
    zp_y = real(ifft2(KY.*z_plus));
    zm_y = real(ifft2(KY.*z_minus));
    Lzp_x = real(ifft2(KX.*Lap_z_plus));
    Lzm_x = real(ifft2(KX.*Lap_z_minus));
    Lzp_y = real(ifft2(KY.*Lap_z_plus));
    Lzm_y = real(ifft2(KY.*Lap_z_minus));
    
    PB_zp_Lzm = fft2((zp_x.*Lzm_y) - (zp_y.*Lzm_x));
    PB_zm_Lzp = fft2((zm_x.*Lzp_y) - (zm_y.*Lzp_x));
    PB_zp_zm  = fft2((zp_x.*zm_y)  - (zp_y.*zm_x));
    
    NL_Sup = -(0.5).*(PB_zp_Lzm + PB_zm_Lzp).*dealias;
    NL_Lap = -(0.5).*k2_perp.*PB_zp_zm.*dealias;
    
    NL_plus = NL_Sup - NL_Lap;
    NL_minus = NL_Sup + NL_Lap;
    
    % Compute Linear term
    L_plus = va.*KZ.*Lap_z_plus;
    L_minus = -va.*KZ.*Lap_z_minus;
    
    %%% Compute Solution at the next step %%%
    
    Lap_z_plus_new = dt*(L_plus + NL_plus) + Lap_z_plus;
    Lap_z_minus_new = dt*(L_minus + NL_minus) + Lap_z_minus;
    
    Lap_z_plus_new = Lap_z_plus_new.*exp_correct;
    Lap_z_minus_new = Lap_z_minus_new.*exp_correct;
    
    Lap_z_plus_new(1,1) = 0;
    Lap_z_minus_new(1,1) = 0;
    
    %%% Energy %%%
    
    E_plus_grid = (abs(Lap_z_plus_new).^2)./abs(k2_poisson);
    E_minus_grid = (abs(Lap_z_minus_new).^2)./abs(k2_poisson);
    
    E_plus(n) = sum(E_plus_grid(:))*(grid_int);
    E_minus(n) = sum(E_minus_grid(:))*(grid_int);   
    
    t=t+dt;
    
    %%% Plotting %%%        !!! What variables will I need to plot for RMHD? !!!
%     if (k==TSCREEN)
%         %Go back to real space for plotting
%         zp = real(ifftn(Lap_z_plus_new./k2_poisson));
%         zm = real(ifftn(Lap_z_minus_new./k2_poisson));
%  
%         %subplot(2,1,1)
%         
%         %%% Contour Plot of Zeta_Plus
%         isosurface(permute(zp,[2,1,3]),1); colorbar; shading flat;        %If matrix dimesions don't agree, likely exploded to matrix of NaNs
%         %use imagesc (with transpose matrix) instead
%         title(num2str(t));
%         
%         %%% Contour Plot of Zeta_minus
% %         subplot(2,1,2)
% %         contourslice(permute(zm,[2,1,3]),[],[],10); colorbar; shading flat;
%         drawnow
%        
%        k=0;
%    end
    
    n=n+1;
    Lap_z_plus = Lap_z_plus_new;
    Lap_z_minus = Lap_z_minus_new;
end

%% Energy Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
figure(2)
subplot(1,2,1)
plot(time, E_plus)
title('\zeta^+ "Energy"')
xlabel('Time')
subplot(1,2,2)
plot(time, E_minus)
title('\zeta^- "Energy"')
xlabel('Time')