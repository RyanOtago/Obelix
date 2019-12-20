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

%% Paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

va   = 1;      % Alven velocity
nu   = 1e-3;
beta = 1;      % c_s/v_A

Fullscreen = 0;         % Makes plot figure fullscreen !!! Forces figure to foreground !!!

LX = 2*pi;     % Box-size (x-direction)
LY = 2*pi;     % Box-size (y-direction)
LZ = 2*pi;     % Box-size (z-direction)     %%% !!! Should this scale differently to LX and LY? !!! %%%

NX = 32;      % Resolution in x
NY = 32;      % Resolution in y
NZ = 32;      % Resolution in z
N = NX*NY*NZ;

CFL = 0.01;            %%% !!! Think about CFL conditions !!! %%% (Will eventually implement variable time step)
dt  = 1e-4;     % InitialTime Step
TF  = 1;      % Final Time
TSCREEN = 500;  % Screen Update Interval Count (NOTE: plotting is usually slow)
time = dt:dt:TF;

E_zeta_plus  = zeros(1,length(time));
E_zeta_minus = zeros(1,length(time));
E_s_plus     = zeros(1,length(time));
E_s_minus    = zeros(1,length(time));
u_par_save   = zeros(1,length(time));
u_perx_save  = zeros(1,length(time));
u_pery_save  = zeros(1,length(time));

dx = LX/NX;
dy = LY/NY;
dz = LZ/NZ;
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
k2_poisson(1,1,:) = 1;        % Fixed Laplacian in F.S. for Poisson's equation (Because first entry was 0)
k2 = k2_perp + KZ.^2;         % Laplacian in F.S.

%% Initial Condition %%%%%%%%%%%%%%%%%%%%%%%%%%%

[i,j,k] = ndgrid((1:NX)*dx,(1:NY)*dy,(1:NZ)*dz);
XG = permute(i, [2 1 3]);
YG = permute(j, [2 1 3]);
ZG = permute(k, [2 1 3]);

%%%% Taylor-Green Vortex
% %
% Lap_zeta_plus  = k2_perp.*fftn((1/(pi))*(sin((2*pi/LX)*i).*sin((4*pi/LX)*j).*cos((2*pi/LX)*k) + (1/sqrt(3)).*cos((12*pi/LX)*i).*cos((2*pi/LX)*j).*sin((4*pi/LX)*k)));
% Lap_zeta_minus = k2_perp.*fftn((1/(2*pi))*(sin((2*pi/LX)*i).*sin((2*pi/LX)*j).*cos((8*pi/LX)*k) - (1/sqrt(3)).*cos((2*pi/LX)*i).*cos((4*pi/LX)*j).*sin((2*pi/LX)*k)));
% %
% s_plus  = fftn((-2/sqrt(12*pi))*(1/bpar).*sin((2*pi/LX)*i).*sin((2*pi/LX)*j).*cos((2*pi/LX)*k));
% s_minus = fftn((2/sqrt(12*pi))*(1/bpar).*sin((2*pi/LX)*i).*sin((2*pi/LX)*j).*cos((2*pi/LX)*k));

Lap_zeta_plus  = k2_perp.*fftn(sin(4*pi*(i/LX - 2*j/LY - k/LZ)));
Lap_zeta_minus = k2_perp.*fftn(sin(2*pi*(i/LX + j/LY + k/LZ)));

% Lap_zeta_plus  = k2_perp.*fftn(1);
% Lap_zeta_minus = k2_perp.*fftn(0);

s_plus  = fftn(0.01*cos(-2*pi*(k/LZ+i/LX)));
s_minus = fftn(0.01*sin( 2*pi*(j/LY+k/LZ)));

%% Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
k=0;
n=1;
for i = (1:length(time))
    k=k+1;
    
    %% Update time-step        (!!Requires changing time vector!!)
    u_par = real(ifftn((0.5).*(s_plus + s_minus)));
    
    phi = abs((0.5)*(Lap_zeta_plus + Lap_zeta_minus)./k2_poisson);
    phi_x = real(ifftn(KX.*phi));    % y component of u_perp
    phi_y = real(ifftn(KY.*phi));    %-x component of u_perp
    
    psi = abs((0.5)*(Lap_zeta_plus - Lap_zeta_minus)./k2_poisson);
    psi_x = real(ifftn(KX.*psi));    % y component of u_perp
    psi_y = real(ifftn(KY.*psi));    %-x component of u_perp
    
    u_perp = sqrt((phi_x).^2 + (phi_y).^2);
    u = sqrt(u_par.^2 + u_perp.^2);
    
    va_perp = sqrt(psi_x.^2 + psi_y.^2);
    va_par = real(ifftn(bpar*(0.5).*(s_plus - s_minus)));
    va_time = sqrt(va_par.^2 + va_perp.^2);
    
    dt_grid_phi = (CFL*min_mesh)./u;
    dt_grid_psi = (CFL*min_mesh)./va_time;
    
    dt = min( [min(dt_grid_phi(:)) min(dt_grid_psi(:))])/2;
    exp_correct = exp(dt*nu*k2);  % (Romain thinks k2) %%% !!! Can include linear term !!!
    
    %%% Update zeta p/m for new time-step
    zeta_plus = Lap_zeta_plus./k2_poisson;
    zeta_minus = Lap_zeta_minus./k2_poisson;
    
    %%% Compute Poisson Brackets
    
    % Calculates derivatives for PB
    zp_x  = real(ifftn(KX.*zeta_plus));
    zm_x  = real(ifftn(KX.*zeta_minus));
    zp_y  = real(ifftn(KY.*zeta_plus));
    zm_y  = real(ifftn(KY.*zeta_minus));
    
    sp_x  = real(ifftn(KX.*s_plus));
    sm_x  = real(ifftn(KX.*s_minus));
    sp_y  = real(ifftn(KY.*s_plus));
    sm_y  = real(ifftn(KY.*s_minus));
    
    Lzp_x = real(ifftn(KX.*Lap_zeta_plus));
    Lzm_x = real(ifftn(KX.*Lap_zeta_minus));
    Lzp_y = real(ifftn(KY.*Lap_zeta_plus));
    Lzm_y = real(ifftn(KY.*Lap_zeta_minus));
    
    % Calculates PB
    % Alfven
    PB_zp_Lzm = fftn((zp_x.*Lzm_y) - (zp_y.*Lzm_x));
    PB_zm_Lzp = fftn((zm_x.*Lzp_y) - (zm_y.*Lzp_x));
    PB_zp_zm  = fftn((zp_x.*zm_y)  - (zp_y.*zm_x));
    % Compressive
    PB_zp_sp  = fftn((zp_x.*sp_y)  - (zp_y.*sp_x));
    PB_zp_sm  = fftn((zp_x.*sm_y)  - (zp_y.*sm_x));
    PB_zm_sp  = fftn((zm_x.*sp_y)  - (zm_y.*sp_x));
    PB_zm_sm  = fftn((zm_x.*sm_y)  - (zm_y.*sm_x));
    
    NL_zeta_Sup = -(0.5).*(PB_zp_Lzm + PB_zm_Lzp).*dealias;
    NL_zeta_Lap = -(0.5).*k2_perp.*PB_zp_zm.*dealias;
    NL_zeta_plus  = NL_zeta_Sup - NL_zeta_Lap;
    NL_zeta_minus = NL_zeta_Sup + NL_zeta_Lap;
    
    NL_s_plus   = -(0.5).*((1-bpar).*PB_zp_sp + (1+bpar).*PB_zm_sp).*dealias;
    NL_s_minus  = -(0.5).*((1+bpar).*PB_zp_sm + (1-bpar).*PB_zm_sm).*dealias;
    
    %%% Compute Linear terms
    
    Lin_zeta_plus  =   va.*KZ.*Lap_zeta_plus;
    Lin_zeta_minus =  -va.*KZ.*Lap_zeta_minus;
    Lin_s_plus     =  (va*bpar).*KZ.*s_plus;
    Lin_s_minus    = -(va*bpar).*KZ.*s_minus;
    
    %%% Compute Solution at the next step %%%
    
    Lap_zeta_plus_new  = dt*(Lin_zeta_plus + NL_zeta_plus) + Lap_zeta_plus;
    Lap_zeta_minus_new = dt*(Lin_zeta_minus + NL_zeta_minus) + Lap_zeta_minus;
    s_plus_new         = dt*(Lin_s_plus + NL_s_plus) + s_plus;
    s_minus_new        = dt*(Lin_s_minus + NL_s_minus) + s_minus;
    
    Lap_zeta_plus_new  = Lap_zeta_plus_new.*exp_correct;
    Lap_zeta_minus_new = Lap_zeta_minus_new.*exp_correct;
    s_plus_new  = s_plus_new.*exp_correct;
    s_minus_new = s_minus_new.*exp_correct;
    
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
    
    %% Plotting %%%
    if (k==TSCREEN)
        
        %Go back to real space for plotting
        zetap  = double(permute(real(ifftn(Lap_zeta_plus_new./k2_poisson)),[2,1,3]));
        zetam  = double(permute(real(ifftn(Lap_zeta_minus_new./k2_poisson)),[2,1,3]));
        sp     = double(permute(real(ifftn(s_plus_new)),[2,1,3]));
        sm     = double(permute(real(ifftn(s_minus_new)),[2,1,3]));
        
        figure(1)
        if Fullscreen == 1
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])        % Makes figure fullscreen
        end
        subplot(2,2,1)
        hold on
        hx = slice(XG, YG, ZG, zetap, LX, [], []);
        set(hx,'FaceColor','interp','EdgeColor','none')
        hy = slice(XG, YG, ZG, zetap, [], dy, []);
        set(hy,'FaceColor','interp','EdgeColor','none')
        hz = slice(XG, YG, ZG, zetap, [], [], LZ);
        set(hz,'FaceColor','interp','EdgeColor','none')
        hold off
        
        daspect([1,1,1])
        axis tight
        box on
        view(42,16)
        camproj perspective
        set(gcf,'Renderer','zbuffer')
        title([num2str(t,'%0.3f') '  \zeta^+'])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        colorbar
        
        subplot(2,2,2)
        hold on
        hx = slice(XG, YG, ZG, zetam, LX, [], []);
        set(hx,'FaceColor','interp','EdgeColor','none')
        hy = slice(XG, YG, ZG, zetam, [], dy, []);
        set(hy,'FaceColor','interp','EdgeColor','none')
        hz = slice(XG, YG, ZG, zetam, [], [], LZ);
        set(hz,'FaceColor','interp','EdgeColor','none')
        hold off
        
        daspect([1,1,1])
        axis tight
        box on
        view(42,16)
        camproj perspective
        set(gcf,'Renderer','zbuffer')
        title('\zeta^-')
        xlabel('x')
        ylabel('y')
        zlabel('z')
        colorbar
        
        subplot(2,2,3)
        hold on
        hx = slice(XG, YG, ZG, sp, LX, [], []);
        set(hx,'FaceColor','interp','EdgeColor','none')
        hy = slice(XG, YG, ZG, sp, [], dy, []);
        set(hy,'FaceColor','interp','EdgeColor','none')
        hz = slice(XG, YG, ZG, sp, [], [], LZ);
        set(hz,'FaceColor','interp','EdgeColor','none')
        hold off
        
        daspect([1,1,1])
        axis tight
        box on
        view(42,16)
        camproj perspective
        set(gcf,'Renderer','zbuffer')
        title('z^+')
        xlabel('x')
        ylabel('y')
        zlabel('z')
        colorbar
        
        subplot(2,2,4)
        hold on
        hx = slice(XG, YG, ZG, sm, LX, [], []);
        set(hx,'FaceColor','interp','EdgeColor','none')
        hy = slice(XG, YG, ZG, sm, [], dy, []);
        set(hy,'FaceColor','interp','EdgeColor','none')
        hz = slice(XG, YG, ZG, sm, [], [], LZ);
        set(hz,'FaceColor','interp','EdgeColor','none')
        hold off
        
        daspect([1,1,1])
        axis tight
        box on
        view(42,16)
        camproj perspective
        set(gcf,'Renderer','zbuffer')
        title('z^-')
        xlabel('x')
        ylabel('y')
        zlabel('z')
        colorbar
        %         saveas(gcf, ['./gif/' num2str(t) '.jpg'])
        drawnow
        k=0;
    end
    
    % Update variables for next timestep
    n = n+1
    Lap_zeta_plus  = Lap_zeta_plus_new;
    Lap_zeta_minus = Lap_zeta_minus_new;
    s_plus         = s_plus_new;
    s_minus        = s_minus_new;
end
disp(toc)
disp(t)
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

end


    