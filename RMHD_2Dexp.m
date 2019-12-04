function RMHD_2D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                               %%%
%%%   Draft 2D RMHD Solver        %%%
%%%                               %%%
%%%   Uses Elsasser Formulation   %%%
%%%   for an Alfven Wave in a     %%%
%%%   2D LX x LY Torus            %%%
%%%                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%%% Parameters %%%

nu = 1e-3;     % Viscosity

LX = 2*pi;     % Box-size (x-direction)
LY = 2*pi;     % Box-size (y-direction)

NX = 128;      % Resolution in x
NY = 128;      % Resolution in y

dt = 1e-4;     % Time Step              !!! Think about CFL conditions !!!
TF = 0.1;      % Final Time
TSCREEN = 50; % Sreen Update Interval Count

time = [0:dt:TF];
E_plus = zeros(1,length(time));
E_minus = zeros(1,length(time));

I=sqrt(-1);
dx = LX/NX;
dy = LY/NY;
N = NX*NY;
dA = dx*dy;
grid_int = dA/N;
t=0.0;

%%%% Initialise wavevector grid %%%%

kx = (2*I*pi/LX)*[0:((NX/2)-1)  -(NX/2):-1];    % [0, 1, ..., NX/2-1, -NX/2, -NX/2+1, ..., -1]    % This is a formatting convention
ky = (2*I*pi/LY)*[0:((NY/2)-1)  -(NY/2):-1];
[KX, KY] = ndgrid(kx, ky);
dealias = abs(KX) < (1/3)*NX & abs(KY) < (1/3)*NY;            % Cutting of frequencies using the 2/3 rule   (2/3)*(N/2)

k2_perp = KX.^2 + KY.^2;      % Laplacian in Fourier space
k2_poisson = k2_perp;    
k2_poisson(1,1) = 1;          % Fixed Laplacian in Fourier space for Poisson's equation   % Because first entry of k2_perp is 0

%%%% Initial Condition %%%%
%%%% !!! Ensure z_plus and z_minus are in Fourier space !!! %%%%

[i,j]=ndgrid((1:NX)*dx,(1:NY)*dy);
Lap_z_plus = k2_perp.*fft2(sin((2*pi*j + 2*pi*i)./LY));
Lap_z_minus = k2_perp.*fft2(cos((8*pi*i + 2*pi*j)./LX));         %%% Make sure either is not constant otherwise there will be no time evolution %%%

k=0;
n=1;
while t<TF
    k=k+1;
    
    z_plus = Lap_z_plus./k2_poisson;        
    z_minus = Lap_z_minus./k2_poisson;      
    
    % Computes Poisson Brackets (RHS of Schekochihin-09 Eq (21))
    % NL -> "Non-Linear"
    NL_Sup = -0.5.*((Poisson(z_plus, Lap_z_minus, KX, KY) + Poisson(z_minus, Lap_z_plus, KX, KY)).*dealias); 
    NL_Lap = k2_perp.*(Poisson(z_plus, z_minus, KX, KY).*dealias);
    
    NL_plus = NL_Sup - NL_Lap;
    NL_minus = NL_Sup + NL_Lap;
       
    %%% Compute Solution at the next step %%%
    
    Lap_z_plus_new = dt*NL_plus + Lap_z_plus;
    Lap_z_minus_new = dt*NL_minus + Lap_z_minus;
    
    exp_correct = exp(dt*nu*k2_perp);
    Lap_z_plus_new = Lap_z_plus_new.*exp_correct;
    Lap_z_minus_new = Lap_z_minus_new.*exp_correct;
    
    t=t+dt;
    
    %%% Energy %%%
    
    E_plus_grid = abs(k2_perp).*(abs(Lap_z_plus_new./k2_poisson).^2);
    E_minus_grid = abs(k2_perp).*(abs(Lap_z_minus_new./k2_poisson).^2);
    
    E_plus(n) = sum(sum(E_plus_grid))*(grid_int);
    E_minus(n) = sum(sum(E_minus_grid))*(grid_int);    % Large spikes were a result of unnecesary inverse FT of Energy
    
    %%% Plotting %%%
    
    if (k==TSCREEN)
        
        % Go back to real space for plotting
        zp = real(ifft2(Lap_z_plus_new./k2_poisson));
        zm = real(ifft2(Lap_z_minus_new./k2_poisson));

        %%% Contour Plot of Zeta_Plus
        subplot(1,2,1)
        contourf(zp',50,'LineColor','none'); colorbar; shading flat;        %If matrix dimesions don't agree, likely exploded to matrix of NaNs
        % use imagesc (with transpose matrix) instead
        title(['\zeta^+     ' num2str(t)]);
        pbaspect([LX LY 1])
        
        %%% Contour Plot of Zeta_minus
        subplot(1,2,2)
        contourf(zm',50,'LineColor','none'); colorbar; shading flat;
        title('\zeta^-')
        pbaspect([LX LY 1])
        
        drawnow
        k=0;
    end
    
    Lap_z_plus = Lap_z_plus_new;
    Lap_z_minus = Lap_z_minus_new;
    n=n+1;
end

%%% Energy Plot %%%

figure(2)
plot(time, E_plus, time, E_minus)
legend('\zeta^{+}', '\zeta^{-}')
title('Conserved quantities (in theory)')
xlabel('Time')