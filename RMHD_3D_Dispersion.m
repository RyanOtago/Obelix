function [omegaout] = RMHD_3D_Dispersion(kn, LX, NX, bpar, va, nu, TFrun, no_perp, no_alfven, strengthpx, strengthmx, strengthpy, strengthmy)

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

Cutoff = 1000000;
CFL    = 0.13;      % Courant Number

LY = LX;     % Box-size (y-direction)
LZ = LX;     % Box-size (z-direction)     %%% !!! Should this scale differently to LX and LY? !!! %%%

NY = NX./2;      % Resolution in y
NZ = NX./2;      % Resolution in z

for stnn = 1:length(strengthpx)
    for bparnn = 1:length(bpar)
        for meshnn = 1:length(NX)
            for knnn = 1:length(kn)
                
                N = NX(meshnn)*NY(meshnn)*NZ(meshnn);
                TF = TFrun(bparnn, knnn, stnn);
                TScreen = 100;
                SlowModes = 1;
                Fullscreen = 1;
                
                time         = zeros(1, Cutoff+1);    % +1 accounts for t=0 (incase we reach the cutoff)
                dt_save      = zeros(1, Cutoff);
                E_z_plus     = zeros(1, Cutoff);
                E_z_minus    = zeros(1, Cutoff);
                E_s_plus     = zeros(1, Cutoff);
                E_s_minus    = zeros(1, Cutoff);
                u_par_save   = zeros(1, Cutoff);
                u_perx_save  = zeros(1, Cutoff);
                u_pery_save  = zeros(1, Cutoff);
                
                dx = LX(knnn)/NX(meshnn);
                dy = LY(knnn)/NY(meshnn);
                dz = LZ(knnn)/NZ(meshnn);
                dV = dx*dy*dz;
                
                I=sqrt(-1);
                grid_int = dV/N;
                t=0.0;
                
                %% Initialise wavevector grid %%%%%%%%%%%%%%%%%%
                
                kx = (2*I*pi/LX(knnn))*[0:((NX(meshnn)/2)-1)  -(NX(meshnn)/2):-1];    % [0, 1, ..., NX/2-1, -NX/2, -NX/2+1, ..., -1]    % This is a formatting convention
                ky = (2*I*pi/LY(knnn))*[0:((NY(meshnn)/2)-1)  -(NY(meshnn)/2):-1];
                kz = (2*I*pi/LZ(knnn))*[0:((NZ(meshnn)/2)-1)  -(NZ(meshnn)/2):-1];
                [KX, KY, KZ] = ndgrid(kx, ky, kz);
                idmat = zeros(size(KX));
                idmat(1,1,1) = 1;
                dealias = abs(KX)<(1/3)*NX(meshnn) & abs(KY)<(1/3)*NY(meshnn) & abs(KZ)<(1/3)*NZ(meshnn);    % Cutting of frequencies for dealiasing using the 2/3 rule   (2/3)*(N/2)
                
                k2_perp = KX.^2 + KY.^2;      % (Perpendicular) Laplacian in Fourier space
                k2_poisson = k2_perp;
                k2_poisson(1,1,:) = 1;        % Fixed Laplacian in F.S. for Poisson's equation (Because first entry was 0)
                k2 = k2_perp + KZ.^2;         % Laplacian in F.S.
                kperpmax = max([abs(kx) abs(ky)]);
                kzmax    = max(abs(kz));
                dt = min(CFL/(kzmax*va), TF/1e4);
%                 disp(TF/dt)
%                 dt = TF/1e4
                exp_correct = exp(dt*nu*k2);  % Romain thinks k2 %%%                        !!! Can include linear term !!!
                
                %% Initial Condition %%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                [i,j,k] = ndgrid((1:NX(meshnn))*dx,(1:NY(meshnn))*dy,(1:NZ(meshnn))*dz);
                XG = permute(i, [2 1 3]);
                YG = permute(j, [2 1 3]);
                ZG = permute(k, [2 1 3]);
                
                stpx = (idmat*strengthpx(stnn));
                stpy = (idmat*strengthpy(stnn));
                stmx = (idmat*strengthmx(stnn));
                stmy = (idmat*strengthmy(stnn));
                
%                 if no_alfven == 0
%                     Lap_z_plus  = k2_perp.*fftn(sin(kn(knnn)*i));
%                     Lap_z_minus = k2_perp.*fftn(cos(kn(knnn)*(i+j)));
%                 else
%                     Lap_z_plus  = k2_perp.*fftn(0);
%                     Lap_z_minus = k2_perp.*fftn(0);
%                 end
                
                
                
                if no_perp == 1
                    s_plus  = fftn(0.01*sin(kn(knnn)*k));
                    s_minus = fftn(0);
                else
                    s_plus  = fftn(0.01*sin((kn(knnn))*(i + j + k)));
                    s_minus = fftn(0);
                end
                
                
                
                %% Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                k=1;
                n=1;
                m=0;
                
                
%                 zeta_plus = Lap_z_plus./k2_poisson;
%                 zeta_minus = Lap_z_minus./k2_poisson;
                
                %%% Compute Poisson Brackets
                
                % Calculates derivatives for PB
                if no_perp == 0 && no_alfven == 0
                    zp_x  = real(ifftn(stpx));
                    zm_x  = real(ifftn(stmx));
                    zp_y  = real(ifftn(stpy));
                    zm_y  = real(ifftn(stmy));
                else
                    zp_x  = real(ifftn(KX.*zeta_plus));
                    zm_x  = real(ifftn(KX.*zeta_minus));
                    zp_y  = real(ifftn(KY.*zeta_plus));
                    zm_y  = real(ifftn(KY.*zeta_minus));
                end
                
%                 Lzp_x = real(ifftn(KX.*Lap_z_plus));
%                 Lzm_x = real(ifftn(KX.*Lap_z_minus));
%                 Lzp_y = real(ifftn(KY.*Lap_z_plus));
%                 Lzm_y = real(ifftn(KY.*Lap_z_minus));
                Lap_z_plus  = KX.*stpx + KY.*stpy;
                Lap_z_minus = KX.*stmx + KY.*stmy;


                Lzp_x = real(ifftn(KX.*Lap_z_plus));
                Lzm_x = real(ifftn(KX.*Lap_z_minus));
                Lzp_y = real(ifftn(KY.*Lap_z_plus));
                Lzm_y = real(ifftn(KY.*Lap_z_minus));
                
                
                while t<TF && n<Cutoff
                    %                     tic
                    %%% Update zeta p/m for new time-step
                    
                    
                    sp_x  = real(ifftn(KX.*s_plus));
                    sm_x  = real(ifftn(KX.*s_minus));
                    sp_y  = real(ifftn(KY.*s_plus));
                    sm_y  = real(ifftn(KY.*s_minus));
                    

                    
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
                    
                    NL_s_plus   = -(0.5).*((1-bpar(bparnn)).*PB_zp_sp + (1+bpar(bparnn)).*PB_zm_sp).*dealias;
                    NL_s_minus  = -(0.5).*((1+bpar(bparnn)).*PB_zp_sm + (1-bpar(bparnn)).*PB_zm_sm).*dealias;
                    
                    %%% Compute Linear terms
                    
                    Lin_zeta_plus  =  va.*KZ.*Lap_z_plus;
                    Lin_zeta_minus = -va.*KZ.*Lap_z_minus;
                    Lin_s_plus     = (va*bpar(bparnn)).*KZ.*s_plus;
                    Lin_s_minus    = -(va*bpar(bparnn)).*KZ.*s_minus;
                    
                    %%% Compute Solution at the next step %%%
                    if t == dt
                        Lap_zeta_plus_new  = dt*(Lin_zeta_plus + NL_zeta_plus);% + Lap_z_plus;
                        Lap_zeta_minus_new = dt*(Lin_zeta_minus + NL_zeta_minus);% + Lap_z_minus;
                    else
                        Lap_zeta_plus_new  = 0*dt*(Lin_zeta_plus + NL_zeta_plus) + Lap_z_plus;
                        Lap_zeta_minus_new = 0*dt*(Lin_zeta_minus + NL_zeta_minus) + Lap_z_minus;
                        
                    end
                    s_plus_new         = dt*(Lin_s_plus + NL_s_plus) + s_plus;
                    s_minus_new        = dt*(Lin_s_minus + NL_s_minus) + s_minus;
                    
                    %                     Lap_zeta_plus_new  = Lap_zeta_plus_new.*exp_correct;
                    %
                    %                     %%% Forgot to disable this, so field was changing
                    %                     over time. Could have (probably) screwed with smallscale
                    %                     wavenumbers
                    %
                    %                     Lap_zeta_minus_new = Lap_zeta_minus_new.*exp_correct;
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
                    if (k == TScreen) && m == 0
                        
                        %Go back to real space for plotting
                        zetap = Lap_zeta_plus_new./k2_poisson;
                        zetam = Lap_zeta_minus_new./k2_poisson;
                        zpx = KX.*zetap;
                        zpy = KY.*zetap;
                        zmx = KX.*zetam;
                        zmy = KY.*zetam;
                        zp  = double(permute(real(ifftn(0.5*(- zpy - zmy))),[2,1,3]));
                        zm  = double(permute(real(ifftn(0.5*(  zpx + zmx))),[2,1,3]));
                        disp([mean(zp(:)) mean(zm(:))])
                        disp([(strengthpy(stnn)-strengthmy(stnn)) strengthpx(stnn)+strengthmx(stnn)])
%                         disp(zpx)
% %                                                 zp  = double(permute(real(ifftn(Lap_zeta_plus_new./k2_poisson)),[2,1,3]));
% %                                                 zm  = double(permute(real(ifftn(Lap_zeta_minus_new./k2_poisson)),[2,1,3]));
%                         if SlowModes == 1
%                             sp     = double(permute(real(ifftn(s_plus_new)),[2,1,3]));
%                             sm     = double(permute(real(ifftn(s_minus_new)),[2,1,3]));
%                         end
%                         figure(1)
%                         if Fullscreen == 1
%                             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])        % Makes figure fullscreen
%                         end
%                         if SlowModes == 1
%                             subplot(2,2,1)
%                         else
%                             subplot(1,2,1)
%                         end
%                         hold on
%                         hx = slice(XG, YG, ZG, zp, LX, [], []);
%                         set(hx,'FaceColor','interp','EdgeColor','none')
%                         hy = slice(XG, YG, ZG, zp, [], dy, []);
%                         set(hy,'FaceColor','interp','EdgeColor','none')
%                         hz = slice(XG, YG, ZG, zp, [], [], LZ);
%                         set(hz,'FaceColor','interp','EdgeColor','none')
%                         hold off
%                         daspect([1,1,1])
%                         axis tight
%                         box on
%                         view(42,16)
%                         camproj perspective
%                         set(gcf,'Renderer','zbuffer')
%                         title([num2str(t,'%0.3f') '  \zeta^+'])
%                         xlabel('x')
%                         ylabel('y')
%                         zlabel('z')
%                         colorbar
%                         
%                         if SlowModes == 1
%                             subplot(2,2,2)
%                         else
%                             subplot(1,2,2)
%                         end
%                         hold on
%                         hx = slice(XG, YG, ZG, zm, LX, [], []);
%                         set(hx,'FaceColor','interp','EdgeColor','none')
%                         hy = slice(XG, YG, ZG, zm, [], dy, []);
%                         set(hy,'FaceColor','interp','EdgeColor','none')
%                         hz = slice(XG, YG, ZG, zm, [], [], LZ);
%                         set(hz,'FaceColor','interp','EdgeColor','none')
%                         hold off
%                         daspect([1,1,1])
%                         axis tight
%                         box on
%                         view(42,16)
%                         camproj perspective
%                         set(gcf,'Renderer','zbuffer')
%                         title('\zeta^-')
%                         xlabel('x')
%                         ylabel('y')
%                         zlabel('z')
%                         colorbar
%                         
%                         if SlowModes == 1
%                             subplot(2,2,3)
%                             hold on
%                             hx = slice(XG, YG, ZG, sp, LX, [], []);
%                             set(hx,'FaceColor','interp','EdgeColor','none')
%                             hy = slice(XG, YG, ZG, sp, [], dy, []);
%                             set(hy,'FaceColor','interp','EdgeColor','none')
%                             hz = slice(XG, YG, ZG, sp, [], [], LZ);
%                             set(hz,'FaceColor','interp','EdgeColor','none')
%                             hold off
%                             
%                             daspect([1,1,1])
%                             axis tight
%                             box on
%                             view(42,16)
%                             camproj perspective
%                             set(gcf,'Renderer','zbuffer')
%                             title('z^+')
%                             xlabel('x')
%                             ylabel('y')
%                             zlabel('z')
%                             colorbar
%                             
%                             subplot(2,2,4)
%                             hold on
%                             hx = slice(XG, YG, ZG, sm, LX, [], []);
%                             set(hx,'FaceColor','interp','EdgeColor','none')
%                             hy = slice(XG, YG, ZG, sm, [], dy, []);
%                             set(hy,'FaceColor','interp','EdgeColor','none')
%                             hz = slice(XG, YG, ZG, sm, [], [], LZ);
%                             set(hz,'FaceColor','interp','EdgeColor','none')
%                             hold off
%                             
%                             daspect([1,1,1])
%                             axis tight
%                             box on
%                             view(42,16)
%                             camproj perspective
%                             set(gcf,'Renderer','zbuffer')
%                             title('z^-')
%                             xlabel('x')
%                             ylabel('y')
%                             zlabel('z')
%                             colorbar
%                         end
%                         saveas(gcf, ['./gif/' num2str(t) '.jpg'])
                        drawnow
                        k=0;
                        m=1;
                    end
                    
                    u_par = real(ifftn((0.5).*(s_plus_new + s_minus_new)));
                    u_par_save(n)  =  u_par(NX(meshnn)/2, NY(meshnn)/2, NZ(meshnn)/2);
                    
                    phi = abs((0.5)*((Lap_zeta_plus_new + Lap_zeta_minus_new)./k2_poisson));
                    phi_x = real(ifftn(KX.*phi)).*dealias;    % y component of u_perp
                    phi_y = real(ifftn(KY.*phi)).*dealias;    %-x component of u_perp
                    u_perx_save(n) = -phi_y(NX(meshnn)/2, NY(meshnn)/2, NZ(meshnn)/2);
                    u_pery_save(n) =  phi_x(NX(meshnn)/2, NY(meshnn)/2, NZ(meshnn)/2);
%                     disp(toc)
                    k=k+1;
                    n=n+1;
                    Lap_z_plus  = Lap_zeta_plus_new;
                    Lap_z_minus = Lap_zeta_minus_new;
                    s_plus         = s_plus_new;
                    s_minus        = s_minus_new;
                end
                %% Calculate omega %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 disp(n)
                [velxmax, tvx_max] = findpeaks(abs(u_perx_save),'MinPeakProminence',0.2e-12);
                
                [velymax, tvy_max] = findpeaks(abs(u_pery_save),'MinPeakProminence', 0.2e-12);
                
                [velzmax, tvz_max] = findpeaks(abs(u_par_save),'MinPeakProminence', 0.2e-12);
                
                %%% Calculate Omega from vel_x
                
                T_nvx = dt*(abs(diff(tvx_max)));
                Tvx = 2*mean(T_nvx);
                
                %%% Calculate Omega from vel_y
                
                T_nvy = dt*(abs(diff(tvy_max)));
                Tvy = 2*mean(T_nvy);
                
                %%% Calculate Omega from vel_z
                
                T_nvz = dt*(abs(diff(tvz_max)));
                Tvz = 2*mean(T_nvz);
                
                %%% Find average Omega from 3 velocity components
                if no_perp == 1
                    Tvect = [Tvz];    %Manually remove nonlinear components from this vector
                    T = nanmean(Tvect);
                    omega = (2*pi)/T;
                else
                    Tvect = [Tvx, Tvy, Tvz];    %Manually remove nonlinear components from this vector
                    T = nanmean(Tvect);
                    omega = (2*pi)/T;
                end
                
                if isnan(omega) == 0
                    w(bparnn,knnn,stnn) = omega;
                else
                    w(bparnn,knnn,stnn) = 0;
                end
                %                 disp(TF - t)
                disp([' - Just Finished Wavenumber: ' num2str(kn(knnn))])
                
            end
        end
        disp(['Just Finished bpar = ' num2str(bpar(bparnn))])
    end
    disp(['Just finished strength #' num2str(stnn)])
end
omegaout = w;
end
