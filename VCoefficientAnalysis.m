%% Data Directory %%%
Directory = '';%'./Turbulence/';
Folder    = '';%'2020-02-14 15-38-23/';
Number = 366;     % Number of file we're looking at

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
    s_plus  = output.sp;
    s_minus = output.sm;
catch
    SlowModes = 0;
    s_plus  = 0;
    s_minus = 0;
end

%% Extracting Variables from  Data

sp = ifftn(s_plus);
sm = ifftn(s_minus);

delB =  (1/(2*va))*(1+(1/beta)^2)^(-0.5)*(sp - sm);  % deltaB_parallel/B_0
delU =  (1/(2*va))*(sp + sm);                        % deltaU_parallel/v_A
delR = -(1/beta^2)*delB;                                      % Density
delP =  gamma*delR;                                           % Pressure

%% Linear Interpolation
lstart = [0.8 0.5 0.5];
ldir = [-0.1 sqrt(0.5) -pi/4]; 
L = [LX, LY, LZ];

[delB_line, length_line] = Interpolate(delB, dx, LX, lstart, ip, jp, kp, L, ldir);
[delU_line, length_line] = Interpolate(delU, dx, LX, lstart, ip, jp, kp, L, ldir);
[delR_line, length_line] = Interpolate(delR, dx, LX, lstart, ip, jp, kp, L, ldir);
[delP_line, length_line] = Interpolate(delP, dx, LX, lstart, ip, jp, kp, L, ldir);

% Coefficient Magnitudes
XiMag_line  = abs(delR_line)./abs(delB_line);
ChiMag_line = abs(delU_line)./abs(delB_line);
PsiMag_line = abs(delP_line)./abs(delB_line);

% Coefficient Phases

% n = length(delB_line);
% time = [0:(n-1)]*dx;        % 'time' array for wavelet spectrum
% pad = 1;                    % Pads time-series with zeroes for transform
% dj = 0.25;                  % Fraction of sub-octaves per octave
% s0 = 2*dt;                  % Starting scale
% j1 = 7/dj;                  % Do 7 powers-of-two with dj sub-octaves each
% mother = 'Morlet';          % Wavelet form
% 
% % Wavelet Transform
% [delB_WT, periodB, scaleB, coiB] = wavelet(delB_line, dx, pad, dj, s0, j1, mother);
% [delU_WT, periodU, scaleU, coiU] = wavelet(delU_line, dx, pad, dj, s0, j1, mother);
% [delR_WT, periodR, scaleR, coiR] = wavelet(delR_line, dx, pad, dj, s0, j1, mother);
% [delP_WT, periodP, scaleP, coiP] = wavelet(delP_line, dx, pad, dj, s0, j1, mother);
% 
% C_UB = conj(delU_WT).*delB_WT;
% C_RB = conj(delR_WT).*delB_WT;
% C_PB = conj(delP_WT).*delB_WT;

%% Interpolation Function

function [out, lvec] = Interpolate(delA, dx, LX, lstart, i, j, k, L, ldir)

ldir = ldir/norm(ldir); % Direction along which to look
dl = dx;
total_length = LX^2 ;
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
