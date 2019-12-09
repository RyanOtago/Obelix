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

va   = 1;      % Alven velocity
nu   = 1e-3;
beta = 1;         % 

LX = 2*pi;     % Box-size (x-direction)
LY = 2*pi;     % Box-size (y-direction)
LZ = 2*pi;     % Box-size (z-direction)     %%% !!! Should this scale differently to LX and LY? !!! %%%

NX = 32;      % Resolution in x
NY = 32;      % Resolution in y
NZ = 64;      % Resolution in z
N = NX*NY*NZ;

CFL = 0.2;      
dt  = 1e-3;     % InitialTime Step              %%% !!! Think about CFL conditions !!! %%% (Will eventually implement variable time step)
TF  = 2;      % Final Time
TSCREEN = 100;  % Sreen Update Interval Count (NOTE: plotting is usually slow)
time = [dt:dt:TF];

E_zeta_plus  = zeros(1,length(time));
E_zeta_minus = zeros(1,length(time));
E_s_plus  = zeros(1,length(time));
E_s_minus = zeros(1,length(time));

dx = LX/NX;
dy = LY/NY;
dz = LZ/NZ;
x_vec = [0:dx:LX-dx];
y_vec = [0:dy:LY-dy];
z_vec = [0:dz:LZ-dz];
min_mesh = min([dx dy dz]);
dV = dx*dy*dz;

I=sqrt(-1);
bpar = 1/sqrt(1+(1/beta)^2);
grid_int = dV/N;
t=0.0;

%% Initialise wavevector grid %%%%%%%%%%%%%%%%%%

kx = (2*I*pi/LX)*[0:((NX/2)-1)  -(NX/2):-1];    % [0, 1, ..., NX/2-1, -NX/2, -NX/2+1, ..., -1]    % This is a formatting convention
ky = (2*I*pi/LY)*[0:((NY/2)-1)  -(NY/2):-1];
kz = (2*I*pi/LZ)*[0:((NZ/2)-1)  -(NZ/2):-1];
[KX, KY, KZ] = ndgrid(kx, ky, kz);

dealias = abs(KX)<(1/3)*NX & abs(KY)<(1/3)*NY & abs(KZ)<(1/3)*NZ;    % Cutting of frequencies for dealiasing using the 2/3 rule   (2/3)*(N/2)

k2_perp = KX.^2 + KY.^2;      % (Perpendicular) Laplacian in Fourier space
k2_poisson = k2_perp;    
k2_poisson(1,1,:) = 1;        % Fixed Laplacian in F.S. for Poisson's equation   % Because first entry was 0
k2 = k2_perp + KZ.^2;         % Laplacian in F.S.
exp_correct = exp(dt*nu*k2);  % !!! Should this be k2 or k2_perp ??? (Seems like it should be k2 as will have z-dependence)

%% Initial Condition %%%%%%%%%%%%%%%%%%%%%%%%%%%

[i,j,k] = ndgrid((0:NX-1)*dx,(0:NY-1)*dy,(0:NZ-1)*dz);

Lap_zeta_plus  = k2_perp.*fftn(0.001*sin(2*pi*(i/LX + j/LY + k/LZ)/3));
Lap_zeta_minus = k2_perp.*fftn(0.002*cos(2*pi*(i/LX + j/LY + k/LZ)/3));

% Lap_zeta_plus  = k2_perp.*fftn(1);
% Lap_zeta_minus = k2_perp.*fftn(0);

s_plus  = fftn(0.01*cos(2*pi*k/LZ));
s_minus = fftn(0.01*sin(2*pi*j/LZ));

%% Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0;
n=1;
for i = [1:length(time)]
    k=k+1;

    %%% Update zeta p/m for new time-step
    zeta_plus = Lap_zeta_plus./k2_poisson;        
    zeta_minus = Lap_zeta_minus./k2_poisson;
    
%     % Update time-step        (Requires changing time vector!!)
%     u_par = real(ifftn((0.5).*(s_plus + s_minus)));
%     phi = abs((0.5)*(zeta_plus + zeta_minus));
%     u_perp = real(ifftn(sqrt(abs(k2_perp)).*phi));
%     u = sqrt(u_par.^2 + u_perp.^2);
% 
%     dt_grid = (CFL.*min_mesh)./u;
%     dt = abs(min(dt_grid(:)));
%     disp(dt)

    %%% Compute Poisson Brackets
    NL_zeta_Sup = -(0.5).*(Poisson(zeta_plus, Lap_zeta_minus, KX, KY) + Poisson(zeta_minus, Lap_zeta_plus, KX, KY)).*dealias; 
    NL_zeta_Lap = -(0.5).*k2_perp.*Poisson(zeta_plus, zeta_minus, KX, KY).*dealias;
    NL_s_plus   = -(0.5).*((1-bpar).*Poisson(zeta_plus, s_plus, KX, KY) + (1+bpar).*Poisson(zeta_minus, s_plus, KX, KY)).*dealias;
    NL_s_minus  = -(0.5).*((1+bpar).*Poisson(zeta_plus, s_minus, KX, KY) + (1-bpar).*Poisson(zeta_minus, s_minus, KX, KY)).*dealias;
    
    NL_zeta_plus  = NL_zeta_Sup - NL_zeta_Lap;
    NL_zeta_minus = NL_zeta_Sup + NL_zeta_Lap;
    
    %%% Compute Linear terms
    Lin_zeta_plus  = va.*KZ.*Lap_zeta_plus;
    Lin_zeta_minus = -va.*KZ.*Lap_zeta_minus;
    Lin_s_plus     = (va*bpar).*KZ.*s_plus;
    Lin_s_minus    = (va*bpar).*KZ.*s_minus;   
    
    %%% Compute Solution at the next step %%%
    
    Lap_zeta_plus_new  = dt*(Lin_zeta_plus + NL_zeta_plus) + Lap_zeta_plus;
    Lap_zeta_minus_new = dt*(Lin_zeta_minus + NL_zeta_minus) + Lap_zeta_minus;
    s_plus_new         = dt*(Lin_s_plus + NL_s_plus) + s_plus;
    s_minus_new        = dt*(Lin_s_minus + NL_s_minus) + s_minus;
    
    Lap_zeta_plus_new  = Lap_zeta_plus_new.*exp_correct;
    Lap_zeta_minus_new = Lap_zeta_minus_new.*exp_correct;
    s_plus_new  = s_plus_new.*exp_correct;
    s_minus_new = s_minus_new.*exp_correct;

    Lap_zeta_plus_new(1,1)  = 0;
    Lap_zeta_minus_new(1,1) = 0;
    
    %%% Energy %%%
    
    E_zeta_plus_grid  = (abs(Lap_zeta_plus_new).^2)./abs(k2_poisson);
    E_zeta_minus_grid = (abs(Lap_zeta_minus_new).^2)./abs(k2_poisson);
    E_s_plus_grid     = (abs(s_plus_new)).^2;
    E_s_minus_grid    = (abs(s_minus_new)).^2;
    
    E_zeta_plus(n)  = (0.5)*sum(E_zeta_plus_grid(:))*(grid_int);
    E_zeta_minus(n) = (0.5)*sum(E_zeta_minus_grid(:))*(grid_int);
    E_s_plus(n)     = (0.5)*sum(E_s_plus_grid(:))*(grid_int);
    E_s_minus(n)    = (0.5)*sum(E_s_minus_grid(:))*(grid_int);
    
    t=t+dt;
    
    %%% Plotting %%%
%     if (k==TSCREEN)
        %Go back to real space for plotting
%         zetap  = real(ifftn(Lap_zeta_plus_new./k2_poisson));

% %         zppara = squeeze(zetap(NX/2,:,:));
% %         zpperp = squeeze(zetap(:,:,NZ/2));
% 
%         zetam  = real(ifftn(Lap_zeta_minus_new./k2_poisson));
% %         zmpara = zetam(NX/2,:,:);
% %         zmperp = zetam(:,:,NZ/2);
%         
%         sp     = real(ifftn(s_plus_new));
% %         sppara = sp(NX/2,:,:);
% %         spperp = sp(:,:,NZ/2);
%         
%         sm     = real(ifftn(s_minus_new));
% %         smpara = sm(NX/2,:,:);
% %         smperp = sm(:,:,NZ/2);
%         

%         % Contour Plot of Zeta_Plus
%         subplot(2,2,1)
%         VolumePlot(zetap, x_vec, y_vec, z_vec)
% %         contourf(zppara,50,'LineColor','none'); colorbar; shading flat;
%         contourslice(permute(zetap,[2,1,3]),[],[],NZ/2); colorbar; shading flat;   
%         contourslice(permute(zetap,[2,1,3]),NX/2,[],[]);       %If matrix dimesions don't agree, likely exploded to matrix of NaNs
%         title(['\zeta^+' num2str(t)]);
%        view([1 1 1])
%         % Contour Plot of Zeta_minus
%         subplot(2,2,2)
% %         contourf(zpperp,50,'LineColor','none'); colorbar; shading flat;
%         contourslice(permute(zetam,[2,1,3]),[],[],NZ/2); colorbar; shading flat;
%         contourslice(permute(zetam,[2,1,3]),NX/2,[],[]);
%         title('\zeta^-')
%        view([1 1 1])
%         
%         % Contour Plot of z_plus
%         subplot(2,2,3)
%         contourslice(permute(sp,[2,1,3]),[],[],NZ/2); colorbar; shading flat;
%         contourslice(permute(sp,[2,1,3]),NX/2,[],[]);
%         title('z^+')
%         view([1 1 1])
% %         
% %         % Contour Plot of z_minus
%         subplot(2,2,4)
%         contourslice(permute(sm,[2,1,3]),[],[],NZ/2); colorbar; shading flat;
%         contourslice(permute(sm,[2,1,3]),NX/2,[],[]);
%         title('z^-')
%         view([1 1 1])
%         
%         drawnow        
%         k=0;
%     end
    
    n=n+1;
    Lap_zeta_plus  = Lap_zeta_plus_new;
    Lap_zeta_minus = Lap_zeta_minus_new;
    s_plus         = s_plus_new;
    s_minus        = s_minus_new;
end

%% Energy Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
subplot(1,2,1)
plot(time, E_zeta_plus, time, E_zeta_minus)
title('\zeta^{\pm} "Energy"')
legend('\zeta^+', '\zeta^-')
xlabel('Time')
subplot(1,2,2)
plot(time, E_s_plus, time, E_s_minus)
title('z^{\pm} "Energy"')
legend('z^+', 'z^-')
xlabel('Time')