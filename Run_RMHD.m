function Run_RMHD

% Automatically cycles through wavenumbers to test the dispersion
% relationship of linear waves propagating through a 3D torus

no_alfven = 0;      % Only test w = v_a*bpar*k_parallel (Makes u_perp and B_perp 0)
no_perp   = 0;      % Only test z-axis propagation
va   = 1;
nu = 1e-3;
kn = (2*pi)*[0.1, 1, 10, 20, 50, 100];
NX = 32;        % Mesh resolution
beta = [1];
stpx = [0.1, 0.1, 1, 0,  10];
stmx = [0.1, 0.1, 3, 0, -10];
stpy = [0.1, 0  , 0, 10, 10];   
stmy = [0.1, 0  , 0, 10, 0];

bpar = 1./sqrt(1+(1./beta).^2);
        

if no_alfven == 1
    disp('RMHD Dispersion Test !! No Alfven !!')
    TF = 6*pi./(va.*(bpar')*kn);          % 3 times the expected period for each wavenumber to allow (hopefully) atleast 2 complete periods for analysis
    
    if no_perp == 1
        disp('!! z-axis propagation only !!')
        LX = (2*pi)./(kn);
        [omega] = RMHD_3D_Dispersion(kn, LX, NX, bpar, va, nu, TF, no_perp, no_alfven);
    else
        disp('Diagonal Propagation')
        LX = (2*pi)./kn;
        [omega] = RMHD_3D_Dispersion(kn, LX, NX, bpar, va, nu, TF, no_perp, no_alfven);
    end
else
    disp('RMHD Dispersion Test with Alfven Waves')
    
    if no_perp == 1
        disp('!! z-axis propagation only !!')
        TF = 6*pi./(va.*(bpar')*kn);
        LX = (2*pi)./(kn);
        [omega] = RMHD_3D_Dispersion(kn, LX, NX, bpar, va, nu, TF, no_perp, no_alfven);
    else
        disp('Diagonal Propagation')
        TF = 6*pi./((kn/(2*sqrt(3))) + va.*(bpar')*kn);
        LX = (2*pi)./kn;
        [omega] = RMHD_3D_Dispersion(kn, LX, NX, bpar, va, nu, TF, no_perp, no_alfven, stpx, stmx, stpy, stmy);
    end
end


%% Analytic Dispersion

for meshnn = 1:length(NX)
    if no_alfven == 1
        kAn = linspace(kn(1),kn(end),100);
        omegaAn = va.*(bpar')*kAn;
        plot(kAn, omegaAn, kn, squeeze(omega), 'o')
        legend(num2str(beta'))
    else
        if no_perp == 1
            kAn = linspace(kn(1),kn(end),100);
            omegaAn = va.*(bpar')*kAn;
            plot(kAn, omegaAn, kn, squeeze(omega), 'o')
            legend(num2str(beta'))
        else
            kAn = linspace(kn(1),kn(end),100);
            omegaAn = va.*(bpar')*kAn;
            plot(kAn, omegaAn, kn, squeeze(omega), 'o')
            legend(num2str(beta'))
            return
        end
    end
end
