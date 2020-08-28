%% Data Directory %%%
Directory = './Turbulence/';
Folder    = '';%'2020-07-24 11-46-45/';
Number = 98;     % Number of file we're looking at
rng('shuffle')
filename = @(n) [Directory Folder sprintf('%u',n) '.mat'];

%% Read initialdata from 0.mat %%%
Dinit = dir([Directory Folder '*.mat']);
Nfiles = length(Dinit)-1;       % '-1' accounts for 0.mat

Init = load(filename(0));
input = Init.input;

gamma = 5/3;
va = input.Parameters.va; beta = input.Parameters.beta;
KX = input.KX; KY = input.KY; KZ = input.KZ;
NX = input.Parameters.NX; NY = input.Parameters.NY; NZ = input.Parameters.NZ;
LX = input.Parameters.LX; LY = input.Parameters.LY; LZ = input.Parameters.LZ;

dx = LX/NX; dy = LY/NY; dz = LZ/NZ;

k2_perp = KX.^2 + KY.^2;      % (Perpendicular) Laplacian in Fourier space
k2_poisson = k2_perp; k2_poisson(1,1,:) = 1;

[i,j,k] = ndgrid((1:NX)*dx,(1:NY)*dy,(1:NZ)*dz);
[ip,jp,kp] = ndgrid((0:(NX+1))*dx,(0:(NY+1))*dy,(0:(NZ+1))*dz);
XG = permute(i, [2 1 3]); YG = permute(j, [2 1 3]); ZG = permute(k, [2 1 3]);

%% Reading Data from Simulation
Init1 = load(filename(Number));
output = Init1.output;

Lap_z_plus  = output.Lzp;
Lap_z_minus = output.Lzm;
t = output.time;

try
    s_plus  = output.sp; %randn(size(Lap_z_plus));
    s_minus = output.sm; %randn(size(Lap_z_plus));
catch
    SlowModes = 0;
    s_plus  = 0;
    s_minus = 0;
end

%% Extracting Variables from  Data

sp = real(ifftn(s_plus));
sm = real(ifftn(s_minus));

delB =  (1/(2*va))*(1+(1/beta)^2)^(-0.5)*(sp - sm);  % deltaB_parallel/B_0
delU =  (1/(2*va))*(sp + sm);                        % deltaU_parallel/v_A
delR = -(1/beta^2)*delB;                                      % Density
delP =  gamma*delR;                                           % Pressure

%% Linear Interpolation
lstart = [0.4 0.5 0.5];
ldir = [-0.1 sqrt(0.2) (pi^2)]; 
L = [LX, LY, LZ];
% [delU_line, length_line] = Interpolate(delU, dx, LX, lstart, ip, jp, kp, L, ldir);
% [delB_line, length_line] = Interpolate(delB, dx, LX, lstart, ip, jp, kp, L, ldir);
% [delR_line, length_line] = Interpolate(delR, dx, LX, lstart, ip, jp, kp, L, ldir);
% [delP_line, length_line] = Interpolate(delP, dx, LX, lstart, ip, jp, kp, L, ldir);

% Coefficient Magnitudes
% XiMag_line  = rms(delR_line)./rms(delB_line);
% ChiMag = rms(delU_line)./rms(delB_line);
% ChiMagWT = cwt(delU_line)./cwt(delB_line);
% PsiMag_line = rms(delP_line)./rms(delB_line);
% plot(abs(ChiMag_line))
% Coefficient Phases

% n = length(delB_line); 
% time = [0:(n-1)]*dx;        % 'time' array for wavelet spectrum (note actually a length for us)
% pad = 1;                    % Pads time-series with zeroes for transform
% dj = 0.25;                  % Fraction of sub-octaves per octave
% s0 = 2*dx;                  % Starting scale 
% j1 = 7/dj;                  % Do 7 powers-of-two with dj sub-octaves each
% mother = 'Morlet';          % Wavelet form

%% Wavelet Transform

% Plot of Wavelet Coherence
% wcoherence(delU_line,delB_line, seconds(dx), 'PhaseDisplayThreshold', 0.5)
% wcoherence(delU_line,delB_line, 'PhaseDisplayThreshold', 0.2);
% axis([])
% [wcohUB, wcsUB, period, coi] = wcoherence(delU_line,delB_line, seconds(dx), 'PhaseDisplayThreshold', 0.2);

% Trim off scales larger than LX/2 to avoid measuring coherence of 'same
% point'
% period = seconds(period);
% smallscale = find(period<(LX/2));
% wcsUBsmall = wcsUB(smallscale, :);
% wcohUBsmall = wcohUB(smallscale, :);
% periodsmall = period(smallscale);
% 
% UBphase = angle(wcsUBsmall);
% bins = linspace(-pi, pi, 51);
% % histUB = zeros(51, length(smallscale));
% % figure
% hold on
for i = smallscale'
    disp(i)
    UBphasei = UBphase(i,:);
    wcohUBi  = wcohUBsmall(i,:);
    histUB = histwv(UBphasei, wcohUBi, -pi, pi, 51);
    histUB = histUB./max(histUB);
    barh(bins,histUB) 
    drawnow
    pause(0.1)
end
% figure
% subplot(1,2,1)
% semilogx(period, UBphaseW, 'o')
% title(['\beta = ' num2str(beta)])
% subplot(1,2,2)
% 
% [histw histv] = histwv(meanscaleUB, WEIGHT, -pi, pi, 31);
% barh(bins, histw)

% histogram(meanscaleUB)% figure(2)

% axis([0 2 -pi pi])
% wcoherence(delR_line, delB_line)
% figure(3)
% wcoherence(delP_line, delB_line)
%% Interpolation Function

function [out, lvec] = Interpolate(delA, dx, LX, lstart, i, j, k, L, ldir)

ldir = ldir/norm(ldir); % Direction along which to look
dl = dx;
total_length = 10000*LX^2 ;
lvec = [(-total_length/2):dl:(total_length/2 -dl)];
clear plist
plist = zeros(1,length(lvec));
for kkk=1:3
    plist(kkk,:) = lstart(kkk) + lvec*ldir(kkk);
    % Make sure all points are within the grid
    plist(kkk,:) = mod(plist(kkk,:),L(kkk));
end
% Interpolate to points
itype = 'linear';
out = interpn(i, j, k,padBoundaries(delA), plist(1,:), plist(2,:), plist(3,:), itype, NaN).';
end

%% Weighted Histogram

function histw = histwv(v, w, min, max, bins)
%Inputs:
% v - values
% w - weights
% min - minimum value
% max - max value
% bins - number of bins (inclusive)

%Outputs:
%histw - wieghted histogram
%histv (optional) - histogram of values

delta = (max-min)/(bins-1);
subs = round((v-min)/delta)+1;
histw = accumarray(subs',w',[bins,1]);
end

