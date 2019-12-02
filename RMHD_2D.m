function RMHD_2D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                               %%%
%%%   Draft 2D RMHD Solver        %%%
%%%                               %%%
%%%   Uses Elsasser Formulation   %%%
%%%   for an Alfven Wave in a     %%%
%%%   2D Torus                    %%%
%%%                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%%%% Parameters %%%%

va = 5;        % Alven velocity         %%% This should be reasonably large?, No everything else should be small, we chose O(v_A)~1
nu = 0.0001;    % Viscosity

LX = 2*pi;     % Box-size (x-direction)
LY = 2*pi;     % Box-size (y-direction)
NX = 128;      % Resolution in x
NY = 128;      % Resolution in y

dt = 1e-3;     % Time Step              !!! Think about CFL conditions !!!
TF = 10.0;     % Final Time
TSCREEN = 500; % Sreen Update Interval Count (NOTE: plotting is usually slow)

I=sqrt(-1);
dx = LX/NX;
dy = LY/NY;
t=0.;

%%%% Initialise wavevector grid %%%%

kx = (2*I*pi/LX)*[0:((NX/2)-1)  -(NX/2):-1];    % [0, 1, ..., NX/2-1, -NX/2, -NX/2+1, ..., -1]    % This is a formatting convention
ky = (2*I*pi/LY)*[0:((NY/2)-1)  -(NY/2):-1];

[KX, KY] = ndgrid(kx, ky);

dealias = abs(KX) < (1/3)*NX & abs(KY) < (1/3)*NY;            % Cutting of frequencies using the 2/3 rule   (2/3)*(N/2)

k2_perp = KX.^2 + KY.^2;      % Laplacian in Fourier space
k2_poisson = k2_perp;    
k2_poisson(1,1) = 1;          % Fixed Laplacian in Fourier space for Poisson's equation   % Because first entry of k2_perp is 0

%%%% Initial Condition %%%%

        [i,j]=ndgrid((1:NX)*dx,(1:NY)*dy);
        z_plus=fft2(cos(3*pi*j/LX));
        z_minus = fft2(sin(4*pi*i/LX));

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
       
    %%% Compute Solution at the next step %%%
    
    Lap_z_plus_new = NL_plus + Lap_z_plus;
    Lap_z_minus_new = NL_minus + Lap_z_minus;
    
    exp_correct = exp(nu*k2_perp);
    z_plus_new = (Lap_z_plus_new./k2_poisson).*exp_correct;
    z_minus_new = (Lap_z_minus_new./k2_poisson).*exp_correct;
    
    t=t+dt;

    %%% Plotting %%%        !!! What variables will I need to plot for RMHD? !!!
    
    if (k==TSCREEN)
        
        % Go back to real space for plotting
        zp = real(ifft2(z_plus_new));
        zm = real(ifft2(z_minus_new));

        %%% Contour Plot of Zeta_Plus
        subplot(2,1,1)
        contourf(zp',50,'LineColor','none'); colorbar; shading flat;        %If matrix dimesions don't agree, likely exploded to matrix of NaNs
        % use imagesc (with transpose matrix) instead
        title(num2str(t));
        pbaspect([LX LY 1])
        
        %%% Contour Plot of Zeta_minus
        subplot(2,1,2)
        contourf(zm',50,'LineColor','none'); colorbar; shading flat;
        pbaspect([LX LY 1])
        
        drawnow

        
        k=0;
    end
    
    z_plus = z_plus_new;
    z_minus = z_minus_new;
end