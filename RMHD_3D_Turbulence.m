function RMHD_3D_Turbulence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Test for RMHD     %%%
%%%   Turbulence in code   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SlowModes        = 0;         % Calculate evolution of compressive modes in run
TF               = 4;         % Final Time
NormalisedEnergy = 0;         % Scales initial condition so u_perp ~ B_perp ~ 1

SaveOutput       = 1;         % Writes energies, u and B components for each time step to a .txt file
TOutput          = 10;        % Number of iterations before output
OutputDirectory  = './Run';   % Directory .txt file above is saved to

VariableTimeStep = 1;         % Enable variable time step, else dt must be defined below
% VARIABLE time step
Cutoff           = 100000;    % Maximum number of iterations for variable time step
CFL              = 0.10;      % Courant Number
% FIXED time step
dt               = 1e-5;      % Time Step (For fixed time step runs)

% PLOTTING
TScreen          = 0;         % Screen Update Interval Count (NOTE: plotting is usually slow) (Set to 0 for no plotting)
Fullscreen       = 1;         % Makes plot figure fullscreen (Recommended if saving plots) !!! Forces figure to foreground through run !!!
SavePlot         = 0;         % Saves figure as a .jpg file everytime a new plot is created
PlotDirectory    = './gif/';  % Directory the plot is saved to
EnergyPlot       = 1;         % Plots energy when run has completed

%% Paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
va   = 1;      % Alfven velocity
nu   = 1e-3;   % Viscosity          !!! Check which type it is? (Just for label really) !!!
beta = 1;      % c_s/v_A

LX = 2*pi;     % Box-size (x-direction)
LY = 2*pi;     % Box-size (y-direction)
LZ = 2*pi;     % Box-size (z-direction)

NX = 128;       % Resolution in x
NY = 128;       % Resolution in y
NZ = 128;       % Resolution in z
N  = NX*NY*NZ;

if VariableTimeStep == 1
    time         = zeros(1, Cutoff+1);    % +1 accounts for t=0 (incase we reach the cutoff)
    dt_save      = zeros(1, Cutoff);
    E_z_plus     = zeros(1, Cutoff);
    E_z_minus    = zeros(1, Cutoff);
    E_s_plus     = zeros(1, Cutoff);
    E_s_minus    = zeros(1, Cutoff);
    
else
    time = dt:dt:TF;
    
    E_z_plus     = zeros(1, length(time));
    E_z_minus    = zeros(1, length(time));
    E_s_plus     = zeros(1, length(time));
    E_s_minus    = zeros(1, length(time));
end

dx = LX/NX;
dy = LY/NY;
dz = LZ/NZ;
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

dealias = (LX/(2*pi))*abs(KX) < (1/3)*NX & (LY/(2*pi))*abs(KY) < (1/3)*NY & (LZ/(2*pi))*abs(KZ) < (1/3)*NZ;    % Cutting of frequencies for dealiasing using the 2/3 rule   (2/3)*(N/2)

k2_perp = KX.^2 + KY.^2;      % (Perpendicular) Laplacian in Fourier space
k2_poisson = k2_perp;
k2_poisson(1,1,:) = 1;        % Fixed Laplacian in F.S. for Poisson's equation (Because first entry was 0)
k2 = k2_perp + KZ.^2;         % Laplacian in F.S.
kperpmax = max([abs(kx) abs(ky)]);
kzmax    = max(abs(kz));

if VariableTimeStep == 0
    exp_correct = exp(dt*nu*k2);  % (Romain thinks k2) %%% !!! Can include linear term !!!
end

%% Initial Condition %%%%%%%%%%%%%%%%%%%%%%%%%%%

[i,j,k] = ndgrid((1:NX)*dx,(1:NY)*dy,(1:NZ)*dz);
% Grids for plotting
XG = permute(i, [2 1 3]);
YG = permute(j, [2 1 3]);
ZG = permute(k, [2 1 3]);

% Lap_z_plus  = k2_perp.*fftn(0.1*cos(2*pi*(4*i/LX - 2*j/LY - k/LZ)));
% Lap_z_minus = k2_perp.*fftn(0.1*sin(2*pi*(i/LX + j/LY + k/LZ)));

% Lap_zeta_plus  = k2_perp.*fftn(1);
% Lap_zeta_minus = k2_perp.*fftn(0);

if SlowModes == 1
    s_plus  = fftn(0.01*cos(-2*pi*(k/LZ+4*i/LX)));
    s_minus = fftn(0.01*sin( 2*pi*(3*j/LY+k/LZ)));
end

% Another way to create initial condition
k2filter = sqrt(abs(k2)) < 8*pi/LX;
Lap_z_plus = k2_perp.*k2filter.*fftn(randn(NX,NY,NZ));
Lap_z_minus = k2_perp.*k2filter.*fftn(randn(NX,NY,NZ));

% s_plus = k2filter.*fftn(0.1*randn(NX,NY,NZ));
% s_minus = k2filter.*fftn(0.1*randn(NX,NY,NZ));

% if NormalisedEnergy == 1
%     z_plus  = Lap_z_plus./k2_poisson;
%     z_minus = Lap_z_minus./k2_poisson;
%     
%     z_px = ifftn(z_plus.*KX);
%     z_py = ifftn(z_plus.*KY);
%     z_mx = ifftn(z_minus.*KX);
%     z_my = ifftn(z_minus.*KY);
%     
%     E_u_grid = (1/4)*(abs(z_py + z_my).^2 + abs(z_px + z_mx).^2);
%     E_b_grid = (1/(4*va^2))*(abs(-z_py + z_my).^2 + abs(z_px - z_mx).^2);
% 
%     E_u = (0.5)*sum(E_u_grid(:))*(grid_int);
%     E_b = (0.5)*sum(E_b_grid(:))*(grid_int);
% 
% %     disp([E_u E_b])
% end


if SaveOutput == 1
    m=0;
    RunFolder = datestr(datetime('now'), 31);
    RunFolder = strrep(RunFolder, ':', '-');
    mkdir(OutputDirectory, RunFolder)
    
    % Save initial conditions and parameter values
    input.Lzp = Lap_z_plus;
    input.Lzm = Lap_z_minus;
    if SlowModes == 1
        input.sp = s_plus;
        input.sm = s_minus;
    end
    input.KX = KX;
    input.KY = KY;
    input.KZ = KZ;    
    
    Parameters = struct('va', va, 'nu', nu, 'beta', beta, 'LX', LX, 'LY', LY, 'LZ', LZ, 'NX', NX, 'NY', NY, 'NZ', NZ, 'dtFixed', dt, 'CFL', CFL, 'TF', TF, 'TOutput', TOutput);
    input.Parameters = Parameters;
    save([OutputDirectory '/' RunFolder '/' num2str(m)], 'input')
end

%% Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0;    % Plotting interval
l=0;    % Output interval
n=1;

while t<TF && n<Cutoff
    k=k+1;
    l=l+1;
    
    if VariableTimeStep == 1
        %% Update time-step        (!!Requires changing time vector!!)
        
        phi = abs((0.5)*(Lap_z_plus + Lap_z_minus)./k2_poisson);
        phi_x = real(ifftn(KX.*phi));    % y component of u_perp
        phi_y = real(ifftn(KY.*phi));    %-x component of u_perp
        
        u_perp = sqrt((phi_x).^2 + (phi_y).^2);
        u_time = abs(u_perp);
        
        psi = abs((0.5)*(Lap_z_plus - Lap_z_minus)./k2_poisson);
        psi_x = real(ifftn(KX.*psi));    % y component of b_perp
        psi_y = real(ifftn(KY.*psi));    %-x component of b_perp
        
        va_perp = sqrt(psi_x.^2 + psi_y.^2);
        va_time = abs(va_perp);
        
        gamma_NL = kperpmax.*(u_time + va_time);
        gamma_L  = kzmax*va;
        dt_grid = CFL./(gamma_NL + gamma_L);
        
        dt = min(dt_grid(:));
        
        dt_save(n) = dt;
        time(n+1)  = time(n) + dt;
        exp_correct = exp(dt*nu*k2);  % Romain thinks k2 %%%                        !!! Can include linear term !!!
    end
    
    %%% Update zeta p/m for new time-step
    z_plus = Lap_z_plus./k2_poisson;
    z_minus = Lap_z_minus./k2_poisson;
    
    %%% Compute Poisson Brackets
    
    % Calculates derivatives for PB
    zp_x  = real(ifftn(KX.*z_plus));
    zm_x  = real(ifftn(KX.*z_minus));
    zp_y  = real(ifftn(KY.*z_plus));
    zm_y  = real(ifftn(KY.*z_minus));
    
    if SlowModes == 1
        sp_x  = real(ifftn(KX.*s_plus));
        sm_x  = real(ifftn(KX.*s_minus));
        sp_y  = real(ifftn(KY.*s_plus));
        sm_y  = real(ifftn(KY.*s_minus));
    end
    
    Lzp_x = real(ifftn(KX.*Lap_z_plus));
    Lzm_x = real(ifftn(KX.*Lap_z_minus));
    Lzp_y = real(ifftn(KY.*Lap_z_plus));
    Lzm_y = real(ifftn(KY.*Lap_z_minus));
    
    % Calculates PB
    % Alfven
    PB_zp_Lzm = fftn((zp_x.*Lzm_y) - (zp_y.*Lzm_x));
    PB_zm_Lzp = fftn((zm_x.*Lzp_y) - (zm_y.*Lzp_x));
    PB_zp_zm  = fftn((zp_x.*zm_y)  - (zp_y.*zm_x));
    
    % Compressive
    if SlowModes == 1
        PB_zp_sp  = fftn((zp_x.*sp_y)  - (zp_y.*sp_x));
        PB_zp_sm  = fftn((zp_x.*sm_y)  - (zp_y.*sm_x));
        PB_zm_sp  = fftn((zm_x.*sp_y)  - (zm_y.*sp_x));
        PB_zm_sm  = fftn((zm_x.*sm_y)  - (zm_y.*sm_x));
    end
    
    NL_z_Sup = -(0.5).*(PB_zp_Lzm + PB_zm_Lzp).*dealias;
    NL_z_Lap = -(0.5).*k2_perp.*PB_zp_zm.*dealias;
    NL_z_plus  = NL_z_Sup - NL_z_Lap;
    NL_z_minus = NL_z_Sup + NL_z_Lap;

    if SlowModes == 1
        NL_s_plus   = -(0.5).*((1-bpar).*PB_zp_sp + (1+bpar).*PB_zm_sp).*dealias;
        NL_s_minus  = -(0.5).*((1+bpar).*PB_zp_sm + (1-bpar).*PB_zm_sm).*dealias;
    end
    %%% Compute Linear terms
    
    Lin_z_plus  =   va.*KZ.*Lap_z_plus;
    Lin_z_minus =  -va.*KZ.*Lap_z_minus;

    if SlowModes == 1
        Lin_s_plus     =  (va*bpar).*KZ.*s_plus;
        Lin_s_minus    = -(va*bpar).*KZ.*s_minus;
    end
    %%% Compute Solution at the next step %%%
    
    Lap_z_plus_new  = dt*(Lin_z_plus + NL_z_plus) + Lap_z_plus;
    Lap_z_minus_new = dt*(Lin_z_minus + NL_z_minus) + Lap_z_minus;

    Lap_z_plus_new  = Lap_z_plus_new.*exp_correct;
    Lap_z_minus_new = Lap_z_minus_new.*exp_correct;

    if SlowModes == 1
        s_plus_new         = dt*(Lin_s_plus + NL_s_plus) + s_plus;
        s_minus_new        = dt*(Lin_s_minus + NL_s_minus) + s_minus;
        
        s_plus_new  = s_plus_new.*exp_correct;
        s_minus_new = s_minus_new.*exp_correct;
    end
    %%% Energy %%%
    
    E_z_plus_grid  = (abs(Lap_z_plus_new).^2)./abs(k2_poisson);
    E_z_minus_grid = (abs(Lap_z_minus_new).^2)./abs(k2_poisson);
    
    E_z_plus(n)  = (0.5)*sum(E_z_plus_grid(:))*(grid_int);
    E_z_minus(n) = (0.5)*sum(E_z_minus_grid(:))*(grid_int);
    
    if SlowModes == 1
        E_s_plus_grid     = (abs(s_plus_new)).^2;
        E_s_minus_grid    = (abs(s_minus_new)).^2;
        
        E_s_plus(n)     = (0.5)*sum(E_s_plus_grid(:))*(grid_int);
        E_s_minus(n)    = (0.5)*sum(E_s_minus_grid(:))*(grid_int);
    end
    t=t+dt;
    
    %% Plotting %%%         %%% Can add better slow mode switch (get rid of 2 subplots)
    if k == TScreen
        
        %Go back to real space for plotting
        zp  = double(permute(real(ifftn(Lap_z_plus_new./k2_poisson)),[2,1,3]));
        zm  = double(permute(real(ifftn(Lap_z_minus_new./k2_poisson)),[2,1,3]));
        if SlowModes == 1
            sp     = double(permute(real(ifftn(s_plus_new)),[2,1,3]));
            sm     = double(permute(real(ifftn(s_minus_new)),[2,1,3]));
        end
        figure(1)
        if Fullscreen == 1
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])        % Makes figure fullscreen
        end
        if SlowModes == 1
            subplot(2,2,1)
        else
            subplot(1,2,1)
        end
        hold on
        hx = slice(XG, YG, ZG, zp, LX, [], []);
        set(hx,'FaceColor','interp','EdgeColor','none')
        hy = slice(XG, YG, ZG, zp, [], dy, []);
        set(hy,'FaceColor','interp','EdgeColor','none')
        hz = slice(XG, YG, ZG, zp, [], [], LZ);
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
        
        if SlowModes == 1
            subplot(2,2,2)
        else
            subplot(1,2,2)
        end
        hold on
        hx = slice(XG, YG, ZG, zm, LX, [], []);
        set(hx,'FaceColor','interp','EdgeColor','none')
        hy = slice(XG, YG, ZG, zm, [], dy, []);
        set(hy,'FaceColor','interp','EdgeColor','none')
        hz = slice(XG, YG, ZG, zm, [], [], LZ);
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
        
        if SlowModes == 1
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
        end
        drawnow
        
        if SavePlot == 1
            saveas(gcf, [PlotDirectory num2str(t) '.jpg'])
        end
        k=0;
    end
    
%% Save Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if l == TOutput && SaveOutput == 1   % To save a new variable simply add line, "output.variable = variable"
        m=m+1;
     
        output.time = t;
        % Zeta and Z
        output.Lzp = Lap_z_plus_new;
        output.Lzm = Lap_z_minus_new;
        if SlowModes == 1
            output.sp = s_plus_new;
            output.sm = s_minus_new;
        end
        % Save Energies
        output.Ezp = E_z_plus;
        output.Ezm = E_z_minus;
        if SlowModes == 1
            output.Esp = E_s_plus;
            output.Esm = E_s_minus;
        end
        
        save([OutputDirectory '/' RunFolder '/' num2str(m)], 'output')
        l=0;
    end

%% Update variables for next timestep %%%%%%%%%%%%%%%%%%%%%
    n = n+1
    Lap_z_plus  = Lap_z_plus_new;
    Lap_z_minus = Lap_z_minus_new;
    if SlowModes == 1
        s_plus         = s_plus_new;
        s_minus        = s_minus_new;
    end
end

%% Energy Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if EnergyPlot == 1
    if VariableTimeStep == 1
        % Cut off trailing zeros from time and energy vectors
        time = time(2:find(time,1,'last'));
        E_z_plus  = E_z_plus(1:length(time));
        E_z_minus = E_z_minus(1:length(time));
        E_s_plus  = E_s_plus(1:length(time));
        E_s_minus = E_s_minus(1:length(time));
    end
    
    figure(2)
    if SlowModes == 1
        subplot(1,2,1)
        plot(time, E_z_plus, time, E_z_minus)
        title('\zeta^{\pm} "Energy"')
        legend('\zeta^+', '\zeta^-', 'Location', 'Best')
        xlabel('Time')
        axis([0 TF 0 1.1*max([E_z_plus E_z_minus])])
        
        subplot(1,2,2)
        plot(time, E_s_plus, time, E_s_minus)
        title('z^{\pm} "Energy"')
        legend('z^+', 'z^-', 'Location', 'Best')
        xlabel('Time')
        axis([0 TF 0 1.1*max([E_s_plus E_s_minus])])
    else
        plot(time, E_z_plus, time, E_z_minus)
        title('\zeta^{\pm} "Energy"')
        legend('\zeta^+', '\zeta^-', 'Location', 'Best')
        xlabel('Time')
        axis([0 TF 0 1.1*max([E_z_plus E_z_minus])])
    end
end

end