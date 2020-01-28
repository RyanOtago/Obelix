function spectrumMATLAB()

Directory = './Run/';
Folder    = '2020-01-28 13-00-00/';

filename = @(n) [Directory Folder sprintf('%u',n) '.mat'];

Dinit = dir([Directory Folder '*.mat']);
Nfiles = length(Dinit);

Init = load(filename(0));
input = Init.input;

KX = input.KX;
KY = input.KY;
KZ = input.KZ;

Kmag = sqrt(abs(KX).^2 + abs(KY).^2 + abs(KZ).^2); % |K|
Kperp = sqrt(abs(KY).^2 + abs(KZ).^2); % |K_perp|
Kpois = Kperp;
Kpois(1,1,:) = 1;
Kprl = abs(KX);
Kspec = Kperp; % Choose Kmag or Kperp
% Bins for k
kgrid = [0:(2*pi/(input.Parameters.LY)):max(abs(KY(:)))].'+1e-4; %

% To hold the spectrum
S.Nk = length(kgrid)-1;
S.kgrid = (kgrid(1:end-1) +  kgrid(2:end))/2;   % Finds mean wavenumber over bin

% Count the number of modes in each bin to normalize later -- this gives a
% smoother result, since we want the average energy in each bin.
oneG = ones(size(KX));
S.nbin = spect1D(oneG,oneG,Kspec,kgrid)*numel(oneG)^2;
S.nnorm = S.nbin./S.kgrid.^2; % k^2 accounts for the fact that in 3D, number of modes in shell increases with k^2
S.nnorm = S.nnorm/mean(S.nnorm); % Normalization by number of modes

m3 = @(a) mean(mean(mean(a)));

% Average over all the snapshots
ns = 0;

fields = {'Lzp','Lzm','EK'};%,'vel3','Bcc1','Bcc2','Bcc3','EK','EM','B','rho'};
for var = fields;S.(var{1}) = 0;end

for nnn = 1: Nfiles-1
%     disp(['Doing ' Folder ' nnn = ' num2str(nnn)])      %<<<<< Change wording
    try 
        D = load(filename(nnn));
    catch 
        warning(['Didnt find the file ' filename(nnn)])
        break
    end
   
    for var = {'Lzp', 'Lzm'}
        %         ft = fftn(D.output.(var{1}));
        if strcmp(var{1}, 'Lzp')
            ft = -(KY.*(D.output.Lzp)./(Kpois.^2)) - (KY.*(D.output.Lzm)./(Kpois.^2));
        else
            ft = (KY.*(D.output.Lzp)./(Kpois.^2)) + (KY.*(D.output.Lzm)./(Kpois.^2));
        end
        S.(var{1}) = S.(var{1}) + spect1D(ft,ft,Kspec,kgrid);
        S.EK = S.EK + S.(var{1}); % Total spectrum is the sum of each component
    end
%     if MHD
%         for var = {'Bcc1','Bcc2','Bcc3'}
%             ft = fftn(D.(var{1}));
%             S.(var{1}) = S.(var{1}) + spect1D(ft,ft,Kspec,kgrid);
%             S.EM = S.EM + S.(var{1});
%         end
%         Bmag = sqrt(D.Bcc1.^2+D.Bcc2.^2+D.Bcc3.^2);
%         ft = fftn(Bmag);
%         S.B = S.B + spect1D(ft,ft,Kspec,kgrid); 
%     end
%     ft = fftn( D.rho - m3(D.rho) );
%     S.rho = S.rho + spect1D(ft,ft,Kspec,kgrid);
%     
    ns = ns+1;
end
for var = fields;S.(var{1}) = S.(var{1})/ns;end
% S.nums = nums;
save(['spectrum.mat'],'S');


%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.4, 0.3, 0.6]);
% if MHD;subplot(211);end
loglog(S.kgrid, S.EK, S.kgrid, S.kgrid.^(-5/3),'k:' )
ylabel('$E_K$','interpreter','latex')
xlabel('$k$','interpreter','latex')
% if MHD
%     hold on 
%     loglog(S.kgrid, S.EM)
%     legend({'$E_K$','$k^{-5/3}$','$E_M$'},'interpreter','latex')
%     subplot(212)
%     loglog(S.kgrid, S.B, S.kgrid, S.rho, S.kgrid,S.kgrid.^(-5/3),'k:' )
%     legend({'$E_{|B|}$','$E_{n}$','$E_{\Delta p}$','$k^{-5/3}$'},'interpreter','latex')
%     xlabel('$k$','interpreter','latex')
%     title('$\beta=10$','interpreter','latex')
% end

% oldpath = path;
% path(oldpath,'~/Research/export-fig')
% set(gcf,'color','w')
% export_fig(['saved-states/spectrum-' folder '.pdf']) 
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