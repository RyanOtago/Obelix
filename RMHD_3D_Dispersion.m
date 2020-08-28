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
CFL    = 0.1;      % Courant Number

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
                TScreen = 0;
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
                dealias = (LX(knnn)/(2*pi))*abs(KX)<(1/3)*NX(meshnn) & (LX(knnn)/(2*pi))*abs(KY)<(1/3)*NY(meshnn) & (LX(knnn)/(2*pi))*abs(KZ)<(1/3)*NZ(meshnn);    % Cutting of frequencies for dealiasing using the 2/3 rule   (2/3)*(N/2)
                
                k2_perp = KX.^2 + KY.^2;      % (Perpendicular) Laplacian in Fourier space
                k2_poisson = k2_perp;
                k2_poisson(1,1,:) = 1;        % Fixed Laplacian in F.S. for Poisson's equation (Because first entry was 0)
                k2 = k2_perp + KZ.^2;         % Laplacian in F.S.
                kperpmax = max([abs(kx) abs(ky)]);
                kzmax    = max(abs(kz));
%                 dt = min([CFL/kzmax TF/1e4]);
%                 disp([TF dt])
%                 dt = TF/1e4
%                 exp_correct = exp(dt*nu*k2);  % Romain thinks k2 %%%                        !!! Can include linear term !!!
                
                %% Initial Condition %%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                [i,j,k] = ndgrid((1:NX(meshnn))*dx,(1:NY(meshnn))*dy,(1:NZ(meshnn))*dz);
                XG = permute(i, [2 1 3]);
                YG = permute(j, [2 1 3]);
                ZG = permute(k, [2 1 3]);
        
                if no_perp == 1
                    s_plus  = fftn(sin(kn(knnn)*k));
                    s_minus = fftn(0);
                else
                    s_plus  = fftn(sin((kn(knnn))*(i + j + k)));
                    s_minus = fftn(0);
                end
                
                %% Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                k=1;
                n=1;
                m=0;
                
                %%% Compute Poisson Brackets
                
                while t<TF && n<Cutoff
                    
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
                    
                    sp_x  = real(ifftn(KX.*s_plus));
                    sm_x  = real(ifftn(KX.*s_minus));
                    sp_y  = real(ifftn(KY.*s_plus));
                    sm_y  = real(ifftn(KY.*s_minus));
                    
                    % Calculates PB
                    % Compressive


                    PB_zp_sp  = fftn((strengthpx(stnn).*sp_y)  - (strengthpy(stnn).*sp_x));

                    PB_zp_sm  = fftn((strengthpx(stnn).*sm_y)  - (strengthpy(stnn).*sm_x));

                    PB_zm_sp  = fftn((strengthmx(stnn).*sp_y)  - (strengthmy(stnn).*sp_x));

                    PB_zm_sm  = fftn((strengthmx(stnn).*sm_y)  - (strengthmy(stnn).*sm_x));
                                      
                    NL_s_plus   = -(0.5).*((1-bpar(bparnn)).*PB_zp_sp + (1+bpar(bparnn)).*PB_zm_sp).*dealias;

                    NL_s_minus  = -(0.5).*((1+bpar(bparnn)).*PB_zp_sm + (1-bpar(bparnn)).*PB_zm_sm).*dealias;
%                     disp(sum(abs(NL_s_plus(:)).^2))
                    %%% Compute Linear terms
                    
                    Lin_s_plus     = (va*bpar(bparnn)).*KZ.*s_plus;

%                     disp(sum(abs(Lin_s_plus(:)).^2))
                    Lin_s_minus    = -(va*bpar(bparnn)).*KZ.*s_minus;
                    
                    %%% Compute Solution at the next step %%%

                    s_plus_new         = dt*(Lin_s_plus + NL_s_plus) + s_plus;

                    s_minus_new        = dt*(Lin_s_minus + NL_s_minus) + s_minus;
                    
                    s_plus_new  = s_plus_new.*exp_correct;
                                                                                                   
                    s_minus_new = s_minus_new.*exp_correct;
                    
                    %%% Energy %%%
                    E_s_plus_grid     = (abs(s_plus_new)).^2;
                    E_s_minus_grid    = (abs(s_minus_new)).^2;

                    E_s_plus(n)     = (0.5)*sum(E_s_plus_grid(:))*(grid_int);
                    E_s_minus(n)    = (0.5)*sum(E_s_minus_grid(:))*(grid_int);
                    
                    t=t+dt;
                    
                    %% Plotting %%%
                    if (k == TScreen) && m == 0
                        
                        %Go back to real space for plotting
                        

                        if SlowModes == 1
                            sp     = double(permute(real(ifftn(s_plus_new)),[2,1,3]));
                            sm     = double(permute(real(ifftn(s_minus_new)),[2,1,3]));
                        end
                        figure(1)
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
                        if SlowModes == 1
%                             subplot(2,2,3)
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
                        end
%                         saveas(gcf, ['./gif/' num2str(t) '.jpg'])
                        drawnow
                        k=0;
                        
                    end
                    
                    u_par = real(ifftn((0.5).*(s_plus_new + s_minus_new)));
                    u_par_save(n)  =  u_par(NX(meshnn)/2, NY(meshnn)/2, NZ(meshnn)/2);
                    
%                     phi = abs((0.5)*((Lap_zeta_plus_new + Lap_zeta_minus_new)./k2_poisson));
%                     phi_x = real(ifftn(KX.*phi)).*dealias;    % y component of u_perp
%                     phi_y = real(ifftn(KY.*phi)).*dealias;    %-x component of u_perp
%                     u_perx_save(n) = -phi_y(NX(meshnn)/2, NY(meshnn)/2, NZ(meshnn)/2);
%                     u_pery_save(n) =  phi_x(NX(meshnn)/2, NY(meshnn)/2, NZ(meshnn)/2);
%                     disp(toc)
                    k=k+1;
                    n=n+1;
                    
                    s_plus         = s_plus_new;
                    s_minus        = s_minus_new;
                end
                %% Calculate omega %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 disp(n)
%                 [velxmax, tvx_max] = findpeaks(abs(u_perx_save),'MinPeakProminence',0.2e-12);
%                 
%                 [velymax, tvy_max] = findpeaks(abs(u_pery_save),'MinPeakProminence', 0.2e-12);
                
                [velzmax, tvz_max] = findpeaks(abs(u_par_save),'MinPeakProminence', 0.2e-12);
                
                %%% Calculate Omega from vel_x
                
%                 T_nvx = dt*(abs(diff(tvx_max)));
%                 Tvx = 2*mean(T_nvx);
%                 
%                 %%% Calculate Omega from vel_y
%                 
%                 T_nvy = dt*(abs(diff(tvy_max)));
%                 Tvy = 2*mean(T_nvy);
                
                %%% Calculate Omega from vel_z
                
                T_nvz = dt*(abs(diff(tvz_max)));
                Tvz = 2*mean(T_nvz);
                
                %%% Find average Omega from 3 velocity components
                if no_perp == 1
                    Tvect = [Tvz];
                    T = nanmean(Tvect);
                    omega = (2*pi)/T;
                else
                    Tvect = Tvz;%[Tvx, Tvy, Tvz];
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
