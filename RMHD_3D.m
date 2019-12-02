function RMHD_3D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                               %%%
%%%   Draft 3D RMHD Solver        %%%
%%%                               %%%
%%%   Uses Elsasser Formulation   %%%
%%%   for an Alfven Wave in a     %%%
%%%   2D Torus                    %%%
%%%                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

va = 1;        % Alven velocity         %%% This should be reasonably large?, No everything else should be small, we chose O(v_A)~1
nu = 0.001;

LX = 2*pi;     % Box-size (x-direction)
LY = 2*pi;     % Box-size (y-direction)
LZ = 2*pi;     % Box-size (z-direction)     %Should this scale differently to LX and LY?

NX = 128;      % Resolution in x
NY = 128;      % Resolution in y
NZ = 128;      % Resolution in z

dt = 1e-3;     % Time Step              !!! Think about CFL conditions !!!
TF = 10.0;    % Final Time
TSCREEN = 10; % Sreen Update Interval Count (NOTE: plotting is usually slow)

%%% !!! We need an Initial Condition !!! %%%

I=sqrt(-1);
dx = LX/NX;
dy = LY/NY;
dz = LZ/NZ;
t=0.;

%%%% Initialise wavevector grid %%%%
kx = (2*I*pi/LX)*[0:((NX/2)-1)  -(NX/2):-1];    % [0, 1, ..., NX/2-1, -NX/2, -NX/2+1, ..., -1]    % This is a formatting convention
ky = (2*I*pi/LY)*[0:((NY/2)-1)  -(NY/2):-1];
kz = (2*I*pi/LZ)*[0:((NZ/2)-1)  -(NZ/2):-1];
% kz = 1;     %%% !!! Should be kz >> k_perp_max? In order for
% k_perp/k_para ~ epsilon assumption to hold    Doesn't matter for 2D

[KX, KY, KZ] = ndgrid(kx, ky, kz);

dealias = abs(KX)<(1/3)*NX & abs(KY)<(1/3)*NY & abs(KZ)<(1/3)*NZ;            % Cutting of frequencies using the 2/3 rule   (2/3)*(N/2)

k2_perp = KX.^2 + KY.^2;      % Laplacian in Fourier space
k2_poisson = k2_perp;    
k2_poisson(1,1) = 1;          % Fixed Laplacian in Fourier space for Poisson's equation   % Because first entry is 0
disp(size(k2_perp))
%%% Initial Condition %%% !!!!!!!!!!!!!!!!!!

[i,j,l]=ndgrid((1:NX)*dx,(1:NY)*dy, (1:NZ)*dz);
z_plus=exp(-((i-(LX/2)).^2+(j-(3*LY/8)).^2+(l-(LZ/2)).^2)/(0.4));
z_minus=exp(-((i-(LX/2)).^2+(j-(3*LY/8)).^2)+(l-(LZ/3)).^2/(0.4));
disp(size(z_plus))
k=0;
while t<TF
    k=k+1;
    
    Lap_z_plus = k2_poisson.*z_plus;        %%% Should this be k2_poisson or k2_perp?
    Lap_z_minus = k2_poisson.*z_minus;      %%%                 ""
    
    % Computes Poisson Brackets (RHS of Schekochihin-09 Eq (21))
    % NL -> "Non-Linear"
    NL_Sup = -0.5.*(Poisson(z_plus, Lap_z_minus, KX, KY) + Poisson(z_minus, Lap_z_plus, KX, KY)).*dealias; 
    NL_Lap = k2_perp.*Poisson(z_plus, z_minus, KX, KY).*dealias;
    
    NL_plus = dt.*(NL_Sup - NL_Lap);
    NL_minus = dt.*(NL_Sup + NL_Lap);
    
    % Compute Linear term
    L_plus = dt.*va.*KZ.*k2_perp.*z_plus;
    L_minus = -dt.*va.*KZ.*k2_perp.*z_minus;
    
    %%% Compute Solution at the next step %%%
    
    Lap_z_plus_new = NL_plus + L_plus + Lap_z_plus;
    Lap_z_minus_new = NL_minus + L_minus + Lap_z_minus;
    
    exp_correct = exp(nu*k2_perp);
    z_plus_new = (Lap_z_plus_new./k2_poisson).*exp_correct;
    z_minus_new = (Lap_z_minus_new./k2_poisson).*exp_correct;
    
    disp(t)
    t=t+dt;
    
    %%% Plotting %%%        !!! What variables will I need to plot for RMHD? !!!
    if (k==TSCREEN)
        % Go back to real space for plotting
        zp = real(ifftn(z_plus_new));
        zm = real(ifftn(z_minus_new));
        
        subplot(2,1,1)
        
        %%% Contour Plot of Zeta_Plus
        contourf(zp.',50,'LineColor','none'); colorbar; shading flat;        %If matrix dimesions don't agree, likely exploded to matrix of NaNs
        % use imagesc (with transpose matrix) instead
        title(num2str(t));
        
        %%% Contour Plot of Zeta_minus
        subplot(2,1,2)
        contourf(zm.',50,'LineColor','none'); colorbar; shading flat;
        drawnow
        
        k=0;
    end
    
    z_plus = z_plus_new;
    z_minus = z_minus_new;
end