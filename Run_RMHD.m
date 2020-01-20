function Run_RMHD
clear all

% Automatically cycles through wavenumbers to test the dispersion
% relationship of linear waves propagating through a 3D torus

no_alfven = 0;      % Only test w = v_a*bpar*k_parallel (Makes u_perp and B_perp 0)
no_perp   = 0;      % Only test z-axis propagation
va   = 1;
nu = 1e-3;
kn = (2*pi)*[1, 10];%, 50, 100];%, 500, 1000];
NX = 32;        % Mesh resolution
beta = [1];%, 5, 10];%, 10];%, 10];%[0.01, 0.1, 1, 10, 100];

if no_alfven == 0 && no_perp == 0
    stpx = [10];%, 0];%, 10];%,  10];          % Note max of 7 (unless add more colours to 'colour' vector)
    stmx = [0];%, 0];%, 10];%,  10];
    stpy = [0];%, -2];%,  0];%,  10];
    stmy = [0];%, 2];%,  0];%,   0];
else
    stpx = 1;
    stmx = 1;
    stpy = 1;
    stmy = 1;
end

bpar = 1./sqrt(1+(1./beta).^2);


if no_alfven == 1
    disp('RMHD Dispersion Test !! No Alfven !!')
    TF = 6*pi./(va.*(bpar')*kn);          % 3 times the expected period for each wavenumber to allow (hopefully) atleast 2 complete periods for analysis
    
    if no_perp == 1
        disp('!! z-axis propagation only !!')
        LX = (2*pi)./(kn);
        [omega] = RMHD_3D_Dispersion(kn, LX, NX, bpar, va, nu, TF, no_perp, no_alfven, stpx, stmx, stpy, stmy);
    else
        disp('Diagonal Propagation')
        LX = (2*pi)./kn;
        [omega] = RMHD_3D_Dispersion(kn, LX, NX, bpar, va, nu, TF, no_perp, no_alfven, stpx, stmx, stpy, stmy);
    end
else
    disp('RMHD Dispersion Test with Alfven Waves')
    
    if no_perp == 1
        disp('!! z-axis propagation only !!')
        TF = 6*pi./(va.*(bpar')*kn);
        LX = (2*pi)./(kn);
        [omega] = RMHD_3D_Dispersion(kn, LX, NX, bpar, va, nu, TF, no_perp, no_alfven, stpx, stmx, stpy, stmy);
    else
        disp('Diagonal Propagation')
        for ii = 1:length(stpx)
            for jj = 1:length(kn)
                for kk = 1:length(bpar)
                    TF(kk, jj, ii) = 6*pi./((kn(jj)/(2*sqrt(3)))*(stpx(ii) + stmx(ii) - stpy(ii) - stmy(ii)) + va*bpar(kk)*kn(jj) + (bpar(kk)*kn(jj))/(2*sqrt(3))*(stpx(ii) + stmy(ii) - stmx(ii) - stpy(ii)));
                end
            end
        end
        LX = (2*pi)./kn;
        [omega] = RMHD_3D_Dispersion(kn, LX, NX, bpar, va, nu, TF, no_perp, no_alfven, stpx, stmx, stpy, stmy);
    end
end

%% Analytic Dispersion

if no_alfven == 1
    kAn = linspace(kn(1),kn(end));
    omegaAn = va.*(bpar')*kAn;
    plot(kAn, omegaAn, kn, squeeze(omega), 'o')
    legend(num2str(beta'))
else
    if no_perp == 1
        kAn = linspace(kn(1),kn(end));
        omegaAn = va.*(bpar')*kAn;
        plot(kAn, omegaAn, kn, squeeze(omega), 'o')
        legend(num2str(beta'))
    else
        kAn = linspace(kn(1),kn(end));
        for ii = 1:length(stpx)
            for jj = 1:length(kAn)
                for kk = 1:length(bpar)
                    omegaAn(kk, jj, ii) = ((kAn(jj)/(2*sqrt(3)))*(stpx(ii) + stmx(ii) - stpy(ii) - stmy(ii)) + va*bpar(kk)*kAn(jj) + (bpar(kk)*kAn(jj))/(2*sqrt(3))*(stpx(ii) + stmy(ii) - stmx(ii) - stpy(ii)));
                end
            end
        end
    end
end

if no_alfven == 0 && no_perp == 0
    colourAn = {'b-', 'r-', 'm-', 'y-', 'c-', 'g-', 'k-'};
    colour   = {'bo', 'ro', 'mo', 'yo', 'co', 'go', 'ko'};
    close all
    figure(1)
    hold on
    for kk = 1:length(bpar)
        subplot(1, length(bpar), kk)
        plot(kAn, squeeze(omegaAn(kk, :, :)), kn, squeeze(omega(kk, :, :)),'o')
        legend(num2str([1:length(stpx)]'), 'Location', 'Best')
        title(num2str(bpar(kk)))
    end
end