function mit18336_spectral_ns2d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Navier-Stokes equations in vorticity/stream function formulation on the torus %
% Version 1.0                                                                   %
% (c) 2008 Jean-Christophe Nave - MIT Department of Mathematics                 %             
%  jcnave (at) mit (dot) edu                                                    %
%                                                                               %
%  Dw/Dt = nu.Laplacian(w)                                                      % 
%  Laplacian(psi) = -w                                                          %
%  u = psi_y                                                                    %
%  v =-psi_x                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

nu = 1.0e-3;   % Viscosity

LX = 2*pi;     % Box-size (x-direction)
LY = 2*pi;     % Box-size (y-direction)

NX = 128;      % Resolution in x
NY = 128;      % Resolution in y

dt = 1e-3;     % Time Step  (Think about CFL conditions)
TF = 10.0;    % Final Time
TSCREEN = 5000; % Sreen Update Interval Time

time = [0:dt:TF];
E = zeros(1,length(time));

I=sqrt(-1);
dx = LX/NX;
dy = LY/NY;
N  = NX*NY;
dA = dx*dy;
t=0.0;

initial_condition = 0;   % 0 = 'Vortices', 1 = 'Simple Vortices' or 2 = 'Random' !!!(Note Random blows up almost immediately)!!!

%%%% Initialise wavevector grid %%%%

kx = (2*I*pi/LX)*[0:((NX/2)-1)  -(NX/2):-1];    % [0, 1, ..., NX/2-1, -NX/2, -NX/2+1, ..., -1]    % This is a formatting convention
ky = (2*I*pi/LY)*[0:((NY/2)-1)  -(NY/2):-1];
[KX, KY] = ndgrid(kx, ky);
dealias = abs(KX) < (1/3)*NX & abs(KY) < (1/3)*NY;            % Cutting of frequencies using the 2/3 rule   (2/3)*(N/2)

k2_perp = KX.^2 + KY.^2;      % Laplacian in Fourier space
k2_poisson = k2_perp;    
k2_poisson(1,1) = 1;     % Fixed Laplacian in Fourier space for Poisson's equation   % Because first entry is 0    %Should this be -1?

%%%%%%% CHECK THIS %%%%%%%%%%
%%% I.C. = Random always blows up?! %%%

%%% Define initial vorticity distribution %%%
switch lower(initial_condition)
    case {0}
        disp('Initial Condition: Vortices')
        [i,j]=ndgrid((0:(NX-LX/NX))*dx,(0:(NY-LY/NY))*dy);
        psi = exp(-((i-(LX/2)).^2+(j-(3*LY/8)).^2)/(0.4))+exp(-((i-(LX/2)).^2+(j-(5*LY/8)).^2)/(0.2))-0.5*exp(-((i-(5*LX/8)).^2+(j-(5*LY/8)).^2)/(0.04));% Change w to psi?
        
        %psi = psi/max(max(psi));    %% Normalise stream function, is this necessary? %%
        
    case {1}
        disp('Initial Condition: Simple Vortices')
        [i,j]=ndgrid((1:NX)*dx,(1:NY)*dy);
        psi=0.1*sin(i).*cos(j).^2;
  
    case {2}
        disp('Initial Condition: Random')
        psi_initial=random('unif',-0.5,0.5,NX,NY);      % Change w to psi?  % Need to make sure the initial condition is real (especially if initialising from k-space)
        wh_initial = k2_perp.*fft2(psi_initial);
         psi = real(ifft2(wh_initial)).*dealias;
    otherwise
        disp('!!! Unknown initial conditions !!!');
        return
end



w_hat = k2_perp.*fft2(psi);    % Finds initial vorticity %%% Why don't we have to dealias here? (Because no multiplication?) %%%

k=0;
n=1;
for i = [1:numel(time)]
    k=k+1;
    
    % Compute the stream function and get the velocity and gradient of vorticity
    psi_hat = -w_hat./k2_poisson;  % Solve Poisson's Equation
    
    u   = real(ifft2( KY.*psi_hat));    % Compute  y derivative of stream function ==> u
    v   = real(ifft2(-KX.*psi_hat));    % Compute -x derivative of stream function ==> v
    w_x = real(ifft2( KX.*w_hat));      % Compute  x derivative of vorticity
    w_y = real(ifft2( KY.*w_hat));      % Compute  y derivative of vorticity
    
    % Important to use 'real' in case numerical errors introduce a complex
    % component through the Fourier Transform
    
    conv     = u.*w_x + v.*w_y;         % evaluate the convective derivative (u,v).grad(w)   
    conv_hat = fft2(conv);              % go back to Fourier space
    conv_hat = dealias.*conv_hat;       % Perform spherical dealiasing 2/3 rule
   
    % Compute Solution at the next step
    
    w_hat_new = -dt*conv_hat + w_hat;
    w_hat_new = w_hat_new.*exp(dt*nu*k2_perp);
    
    % Crank-Nicolson
    %w_hat_new = ((1/dt + 0.5*nu*k2_perp)./(1/dt - 0.5*nu*k2_perp)).*w_hat - (1./(1/dt - 0.5*nu*k2_perp)).*conv_hat;
    
    E_grid = abs((w_hat_new.^2)./abs(k2_poisson));
    E(n) = sum(sum(E_grid))*(dA/N);
    
    t=t+dt;
   
    %%% Plotting %%%
    if (k==TSCREEN)
        % Go back to real space for plotting
        w = real(ifft2(w_hat_new));
        figure(1)
        contourf(w',50,'LineColor','none'); colorbar; shading flat;        %If matrix dimesions don't agree, likely exploded to matrix of NaNs
        % use imagesc (with transpose matrix) instead (What is difference between imagesc and contourf
        title(num2str(t));
        drawnow
        %saveas(gcf,['./gif/' num2str(t,'%0.2f') '.png'])
        k=0;
    end
    n=n+1;    
    w_hat=w_hat_new;
end

figure(2)
plot(time, E)
legend('Energy')
title('Conserved quantities (in theory)')
xlabel('Time')