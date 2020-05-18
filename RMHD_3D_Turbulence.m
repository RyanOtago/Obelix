function RMHD_3D_Turbulence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Test for RMHD     %%%
%%%   Turbulence in code   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SlowModes        = 1;         % Calculate evolution of compressive modes in run
TF               = 20;        % Final Time
NormalisedEnergy = 1;         % Scales initial condition so u_perp ~ B_perp ~ 1
HyperViscosity   = 1;         % Use nu*(k^6) instead of nu*(k^2) for dissipation

SaveOutput       = 0;         % Writes energies, u and B components for each time step to a .mat file
TOutput          = 200;       % Number of iterations before output
OutputDirectory  = './Turbulence';   % Directory .mat file above is saved to

% Time step
VariableTimeStep = 1;         % Enable variable time step, else dt must be defined below
% Variable
Cutoff           = 1000000;   % Maximum number of iterations for variable time step
dtCutoff         = 1e-4;      % If dt gets smaller than dtCutoff, run will stop
CFL              = 0.13;      % Courant Number
% Fixed
dt               = 1e-4;      % Time Step (For fixed time step runs)

% PLOTTING
TScreen          = 0;         % Screen Update Interval Count (NOTE: plotting is usually slow) (Set to 0 for no plotting)
Fullscreen       = 0;         % Makes plot figure fullscreen (Recommended if saving plots) !!! Forces figure to foreground through run !!!
SavePlot         = 0;         % Saves figure as a .jpg file everytime a new plot is created
PlotDirectory    = './gif/';  % Directory the plot is saved to
EnergyPlot       = 0;         % Plots energy when run has completed

%% Paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
va   = 1;      % Alfven velocity
nu2   = 0.036;          % Viscosity coefficient for N = 128
nu6 = (1/180)*2e-10;    % Hyperviscosity coefficient for N = 128
beta = 1;      % c_s/v_A
sigma_A = 240;   % Alfven-wave Forcing Strength, sigma=0 turns off forcing
sigma_S = 240;   % Slow-mode Forcing Strength
% Filter for forcing, k2filter, is defined in initial condition
init_energy = 0;

LX = 1;     % Box-size (x-direction)
LY = 1;     % Box-size (y-direction)
LZ = 1;     % Box-size (z-direction)

NX = 128;       % Resolution in x
NY = 128;       % Resolution in y
NZ = 128;       % Resolution in z
N  = NX*NY*NZ;

if VariableTimeStep == 1
    time         = zeros(1, Cutoff+1);    % +1 accounts for t=0 (incase we reach the cutoff)
    dt_save      = zeros(1, Cutoff);
    E_z_plus     = zeros(1, Cutoff);
    E_z_minus    = zeros(1, Cutoff);
    if SlowModes == 1
        E_s_plus     = zeros(1, Cutoff);
        E_s_minus    = zeros(1, Cutoff);
    end
%     E_zp_diss    = zeros(1, Cutoff);
%     E_zm_diss    = zeros(1, Cutoff);
    
else
    time = dt:dt:TF;
    
    E_z_plus     = zeros(1, length(time));
    E_z_minus    = zeros(1, length(time));
    if SlowModes == 1
        E_s_plus     = zeros(1, length(time));
        E_s_minus    = zeros(1, length(time));
    end
%     E_zp_diss    = zeros(1, length(time));
%     E_zm_diss    = zeros(1, length(time));
end

dx = LX/NX;
dy = LY/NY;
dz = LZ/NZ;
dV = dx*dy*dz;

I=sqrt(-1);
bpar = 1/sqrt(1+(1/beta)^2);
grid_int = dV/N;
t=0.0;

% Scale nu according to mesh resolution
if HyperViscosity ==1
    nu = nu6*(128/NX)^(16/3);
else
    nu = nu2*(128/NX)^(4/3);
end

%% Initialise wavevector grid %%%%%%%%%%%%%%%%%%

kx = (2*I*pi/LX)*[0:((NX/2)-1) -(NX/2):-1];    % [0, 1, ..., NX/2-1, -NX/2, -NX/2+1, ..., -1]    % This is a formatting convention
ky = (2*I*pi/LY)*[0:((NY/2)-1) -(NY/2):-1];
kz = (2*I*pi/LZ)*[0:((NZ/2)-1) -(NZ/2):-1];
[KX, KY, KZ] = ndgrid(kx, ky, kz);

dealias = (LX/(2*pi))*abs(KX) < (1/3)*NX & (LY/(2*pi))*abs(KY) < (1/3)*NY & (LZ/(2*pi))*abs(KZ) < (1/3)*NZ;    % Cutting of frequencies for dealiasing using the 2/3 rule   (2/3)*(N/2)

k2_perp = KX.^2 + KY.^2;      % (Perpendicular) Laplacian in Fourier space
k2_poisson = k2_perp;
k2_poisson(1,1,:) = 1;        % Fixed Laplacian in F.S. for Poisson's equation (Because first entry was 0)
k2 = k2_perp + KZ.^2;         % Laplacian in F.S.
k6 = k2.^3;
kperpmax = max([abs(kx) abs(ky)]);
kzmax    = max(abs(kz));


if VariableTimeStep == 0
    if HyperViscosity == 1
        exp_correct = exp(dt*nu*k6);
    else
        exp_correct = exp(dt*nu*k2);
    end
end

%% Initial Condition %%%%%%%%%%%%%%%%%%%%%%%%%%%

[i,j,k] = ndgrid((1:NX)*dx,(1:NY)*dy,(1:NZ)*dz);
% Grids for plotting
XG = permute(i, [2 1 3]);
YG = permute(j, [2 1 3]);
ZG = permute(k, [2 1 3]);

if SlowModes == 0
    s_plus_new  = 1;
    s_minus_new = 1;
end
% Another way to create initial condition
k2filter    = sqrt(abs(k2)) < 5*pi/LX & pi/LX < sqrt(abs(k2)) & abs(KX) > 0  & abs(KY) > 0  & abs(KZ) > 0;
Lap_z_plus  = k2_perp.*k2filter.*fftn(randn(NX,NY,NZ));
Lap_z_minus = k2_perp.*k2filter.*fftn(randn(NX,NY,NZ));

if SlowModes == 1
    s_plus = k2filter.*fftn(0.1*randn(NX,NY,NZ));
    s_minus = k2filter.*fftn(0.1*randn(NX,NY,NZ));
end

if NormalisedEnergy == 1
    kz_plus  = Lap_z_plus./sqrt(k2_poisson);
    kz_minus = Lap_z_minus./sqrt(k2_poisson);
    
    E_zplus_grid  = abs(kz_plus).^2;
    E_zminus_grid = abs(kz_minus).^2;
    E_splus_grid  = abs(s_plus).^2;
    E_sminus_grid = abs(s_minus).^2;
    
    E_zp = (0.5)*sum(E_zplus_grid(:))*(grid_int);
    E_zm = (0.5)*sum(E_zminus_grid(:))*(grid_int);
    E_sp = (0.5)*sum(E_splus_grid(:))*grid_int;
    E_sm = (0.5)*sum(E_sminus_grid(:))*grid_int;
    
    kz_plus  = init_energy*(1/sqrt(E_zp))*kz_plus;
    kz_minus = init_energy*(1/sqrt(E_zm))*kz_minus;
    s_plus   = init_energy*(1/sqrt(E_sp))*s_plus;
    s_minus  = init_energy*(1/sqrt(E_sm))*s_minus;

    Lap_z_plus  = sqrt(k2_perp).*kz_plus;
    Lap_z_minus = sqrt(k2_perp).*kz_minus;
end


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
    
    Parameters = struct('va', va, 'nu', nu, 'beta', beta, 'LX', LX, 'LY', LY, 'LZ', LZ, 'NX', NX, 'NY', NY, 'NZ', NZ, 'dtFixed', dt, 'CFL', CFL, 'TF', TF, 'TOutput', TOutput, 'VariableTimeStep', VariableTimeStep, 'HyperViscosity', HyperViscosity, 'SlowModes', SlowModes);
    input.Parameters = Parameters;
    save([OutputDirectory '/' RunFolder '/' num2str(m)], 'input')
end

clear kz_plus kz_minus E_u_grid E_b_grid input kx ky kz

%% Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0;    % Plotting interval
l=0;    % Output interval
n=1;

while t<TF && n<Cutoff
    k=k+1;
    l=l+1;
    
    if VariableTimeStep == 1
        %% Update time-step
        
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
        
        if dt<dtCutoff
           disp(['Time step got too small at t = ' num2str(t)]) 
           disp('Stopping run...')
           return
        end
        
        dt_save(n) = dt;
        time(n+1)  = time(n) + dt;
        
        if HyperViscosity == 1
            exp_correct = exp(dt*nu*k6);
        else
            exp_correct = exp(dt*nu*k2);
        end
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
    
    %%% Forcing %%%
    
    if sigma_A>0
        force_Ap = k2filter.*fftn(randn(NX,NY,NZ));
        force_Am = k2filter.*fftn(randn(NX,NY,NZ));
        
        force_Ap = sigma_A*(force_Ap./sqrt((sum(abs(force_Ap(:)).^2))*(grid_int)));
        force_Am = sigma_A*(force_Am./sqrt((sum(abs(force_Am(:)).^2))*(grid_int)));
    else 
        force_Ap=0;force_Am=0;
    end
                
    if sigma_S>0 && SlowModes == 1
            force_Sp = sigma_S * k2filter.*fftn(randn(NX,NY,NZ));
            force_Sm = sigma_S * k2filter.*fftn(randn(NX,NY,NZ));   
    else
        force_Sp=0;force_Sm=0;
    end
    
    %%% Compute Solution at the next step %%%
    
    Lap_z_plus_new  = dt*(Lin_z_plus + NL_z_plus) + Lap_z_plus + force_Ap*sqrt(dt); 
    Lap_z_minus_new = dt*(Lin_z_minus + NL_z_minus) + Lap_z_minus + force_Am*sqrt(dt);
    
    Lap_z_plus_new  = Lap_z_plus_new.*exp_correct;
    Lap_z_minus_new = Lap_z_minus_new.*exp_correct;

    if SlowModes == 1
        s_plus_new         = dt*(Lin_s_plus + NL_s_plus) + s_plus + force_Sp*sqrt(dt);
        s_minus_new        = dt*(Lin_s_minus + NL_s_minus) + s_minus + force_Sm*sqrt(dt);
        
        s_plus_new  = s_plus_new.*exp_correct;
        s_minus_new = s_minus_new.*exp_correct;
    end
    
    %%% Conserved Quantities %%%
    
    E_z_plus_grid  = (abs(Lap_z_plus_new).^2)./abs(k2_poisson);
    E_z_minus_grid = (abs(Lap_z_minus_new).^2)./abs(k2_poisson);
    
%     E_zp_diss_grid = abs(k2.*Lap_z_plus_new).^2;
%     E_zm_diss_grid = abs(k2.*Lap_z_minus_new).^2;
    
    E_z_plus(n)  = (0.5)*sum(E_z_plus_grid(:))*(grid_int);
    E_z_minus(n) = (0.5)*sum(E_z_minus_grid(:))*(grid_int);
    
%     E_zp_diss(n) = nu*sum(E_zp_diss_grid(:))*grid_int;
%     E_zm_diss(n) = nu*sum(E_zm_diss_grid(:))*grid_int;
    
    if SlowModes == 1
        E_s_plus_grid     = (abs(s_plus_new)).^2;
        E_s_minus_grid    = (abs(s_minus_new)).^2;
        
        E_s_plus(n)     = (0.5)*sum(E_s_plus_grid(:))*(grid_int);
        E_s_minus(n)    = (0.5)*sum(E_s_minus_grid(:))*(grid_int);
    end
             
    t=t+dt;
    
%% Plotting %%%
    if k == TScreen
        PlotGrid(Lap_z_plus_new, Lap_z_minus_new, s_plus_new, s_minus_new, k2_poisson, Fullscreen, SlowModes, SavePlot, PlotDirectory, XG, YG, ZG, LX, LZ, dy, t, KX)
        figure(2)

        if VariableTimeStep == 1
            % Cut off trailing zeros from time and energy vectors
            timep = time(2:find(time,1,'last'));
            E_z_plusp  = E_z_plus(1:length(timep));
            E_z_minusp = E_z_minus(1:length(timep));
        end
        
        plot(timep, E_z_plusp, timep, E_z_minusp)
        title('|\nabla_{\perp}\zeta^{\pm}|^2  "Energy"')
        legend('\zeta^+', '\zeta^-', 'Location', 'Best')
        xlabel('Time')
        drawnow
        
        k=0;
    end
    
%% Save Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if l == TOutput && SaveOutput == 1   % To save a new variable simply add line, "output.variable = variable"
        m=m+1;  % Output counter
        
        output.time = t;
        if VariableTimeStep == 1
            output.timevec = time;
        end
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
    n = n+1;
    Lap_z_plus  = Lap_z_plus_new;
    Lap_z_minus = Lap_z_minus_new;
    if SlowModes == 1
        s_plus  = s_plus_new;
        s_minus = s_minus_new;
    end
end

%% Energy Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if EnergyPlot == 1
    if VariableTimeStep == 1
        % Cut off trailing zeros from time and energy vectors
        time = time(2:find(time,1,'last'));
        E_z_plus  = E_z_plus(1:length(time));
        E_z_minus = E_z_minus(1:length(time));
        if SlowModes == 1
            E_s_plus  = E_s_plus(1:length(time));
            E_s_minus = E_s_minus(1:length(time));
        end
%         E_zp_diss = E_zp_diss(1:length(time));
%         E_zm_diss = E_zm_diss(1:length(time));
    end
    
    figure(2)
    if SlowModes == 1
        subplot(1,2,1)
        plot(time, E_z_plus, time, E_z_minus)
        title('|\nabla_{\perp}\zeta^{\pm}|^2  "Energy"')
        legend('\zeta^+', '\zeta^-', 'Location', 'Best')
        xlabel('Time')
        axis([0 TF 0 1.1*max([E_z_plus E_z_minus])])
        
        subplot(1,2,2)
        plot(time, E_s_plus, time, E_s_minus)
        title('|\nabla_{\perp}z^{\pm}|^2  "Energy"')
        legend('z^+', 'z^-', 'Location', 'Best')
        xlabel('Time')
        axis([0 TF 0 1.1*max([E_s_plus E_s_minus])])
    else
        plot(time, E_z_plus, time, E_z_minus)
        title('|\nabla_{\perp}\zeta^{\pm}|^2  "Energy"')
        legend('\zeta^+', '\zeta^-', 'Location', 'Best')
        xlabel('Time')
%         axis([0 TF 0 1.1*max([E_z_plus E_z_minus])])
    end
end

end

function PlotGrid(Lap_z_plus_new, Lap_z_minus_new, s_plus_new, s_minus_new, k2_poisson, Fullscreen, SlowModes, SavePlot, PlotDirectory, XG, YG, ZG, LX, LZ, dy, t, KX)

        %Go back to real space for plotting
        zp  = double(permute(real(ifftn(KX.*Lap_z_plus_new./k2_poisson)),[2,1,3]));
        zm  = double(permute(real(ifftn(KX.*Lap_z_minus_new./k2_poisson)),[2,1,3]));
        
        zp = (0.5)*(zp + zm);
%         zp = 0.5*(zp + zm);
        if SlowModes == 1
            sp     = double(permute(real(ifftn(s_plus_new)),[2,1,3]));
            sm     = double(permute(real(ifftn(s_minus_new)),[2,1,3]));
        end
        figure(1)
        clf
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
            saveas(gcf, [PlotDirectory num2str(t, '%10f') '.jpg'])
        end
        
end
