function spectrumMATLAB()

%%% Takes output from RMHD_3D

%%% Data Directory %%%
Directory = './Turbulence/';
Folder    = '2020-05-28 15-07-08/';

filename = @(n) [Directory Folder sprintf('%u',n) '.mat'];

%%% Read initialdata from 0.mat %%%
Dinit = dir([Directory Folder '*.mat']);
Nfiles = length(Dinit)-1;       % '-1' accounts for 0.mat

Init = load(filename(0));
input = Init.input;

KX = input.KX; KY = input.KY; KZ = input.KZ;
[NX, NY, NZ] = size(KX);
Kperp = sqrt(abs(KX).^2 + abs(KY).^2); % |K_perp|
Kspec = Kperp;
k2_poisson = KX.^2 + KY.^2;
k2_poisson(1,1,:) = 1;
% Bins for k
kgrid = (0:(2*pi/(input.Parameters.LY)):max(abs(KY(:)))).'+1e-4; %

% To hold the spectrum
S.Nk = length(kgrid)-1;
S.kgrid = (kgrid(1:end-1) +  kgrid(2:end))/2;   % Finds mean wavenumber over bin

% Count the number of modes in each bin to normalize later -- this gives a
% smoother result, since we want the average energy in each bin.
oneG = ones(size(KX));
S.nbin = spect1D(oneG,oneG,Kspec,kgrid)*numel(oneG)^2;
S.nnorm = S.nbin./S.kgrid.^2; % k^2 accounts for the fact that in 3D, number of modes in shell increases with k^2
S.nnorm = S.nnorm/mean(S.nnorm); % Normalization by number of modes
S.EK = 0;
fields = {'Lzp','Lzm','EK', 'KX', 'KY'};

%%% Initialise Plot %%%
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.4, 0.3, 0.6]);
loglog(S.kgrid, S.kgrid.^(-5/3),'k:')

%%% Calculate Spectra %%%
for nn = 1:50:Nfiles
    clear('S.EK')
    for var = fields;S.(var{1}) = 0;end
    
    try
        D = load(filename(nn));
        disp(['    - ' num2str(nn) ' of ' num2str(Nfiles)])
    catch
        warning(['Didnt find the file ' filename(nn) '.mat'])
        break
    end
    
    for var = {'KX', 'KY'}
        zp = D.output.Lzp./k2_poisson;
        zm = D.output.Lzm./k2_poisson;
        uperp = (0.5)*input.(var{1}).*(zp + zm);
        S.(var{1}) = S.(var{1}) + spect1D(uperp,uperp,Kspec,kgrid);
        S.EK = S.EK + S.(var{1})/2; % Total spectrum is the sum of each component
        S.EK = S.EK.*S.nnorm;
    end
    
    hold on
    loglog(S.kgrid, S.EK)
    ylabel('$E_K$','interpreter','latex')
    xlabel('$k$','interpreter','latex')
    drawnow
    
end
plot([2*pi*NX/3 2*pi*NX/3], [1e-10, 1e2], 'k:')
end

function out = spect1D(v1,v2,K,kgrid)
% Function to find the spectrum <v1 v2>, 
% K is the kgrid associated with v1 and v2
% kgrid is the grid for spectral shell binning

nk = length(kgrid)-1;
out = zeros(nk,1);
NT2 = numel(K)^2;
for kk = 1:nk
    out(kk) = sum( real(v1(K<kgrid(kk+1) & K>kgrid(kk)).*conj(v2(K<kgrid(kk+1) & K>kgrid(kk)))) )/NT2;
end
end