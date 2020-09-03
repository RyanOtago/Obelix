%% Data Directory %%%
Directory = '';
Folder    = 'Beta1Test/';
Number = 1;     % Number of file we're looking at
mnum   = ['gm0'; 'gp0'; 'gm1'; 'gp1'; 'gm2'; 'gp2']; 
filename = @(n, mn) [Directory Folder 'test00_' mnum(mn,:) '.00' sprintf('%u',n) '.vtk'];

beta = 1;

va   = 1; % Alfven Speed
vth  = 1; % Thermal Speed (I THINK I NEED TO DEFINE IN TERMS OF BETA AND VA?)
tau  = 1; % Ratio of Ion to Electron Temperatures
Z    = 1; % Ion to Electron charge ratio

sigma = 1 + tau/Z + 1/beta + sqrt((1+tau/Z)^2 + 1/(beta^2));

kappa = 1/(sigma - (2/(sigma*beta))*(1+tau/Z));     % Might not quite be kappa, would be good to check!!!!

%% Read Data from File

% Gm0
V = struct();
V = readVTKast(filename(Number,1), V, mnum(1,:));  % Convert from .vtk to .mat

NX = V.nx; NY = V.ny; NZ = V.nz;             % Number of grid-points in each direction  (Only need to do this from one file)
LX = V.x(end); LY = V.y(end); LZ = V.z(end); % Length of box in each direction
dx = LX/NX; dy = LY/NY; dz = LZ/NZ;          % Grid spacing in each direction

[i,j,k] = ndgrid((1:NX)*dx,(1:NY)*dy,(1:NZ)*dz);
[ip,jp,kp] = ndgrid((0:(NX+1))*dx,(0:(NY+1))*dy,(0:(NZ+1))*dz);
XG = permute(i, [2 1 3]); YG = permute(j, [2 1 3]); ZG = permute(k, [2 1 3]);

Gm0 = V.(mnum(1,:));

% Gp0
V = readVTKast(filename(Number,2), V, mnum(2,:));
Gp0 = V.(mnum(2,:));

% Gm1
V = readVTKast(filename(Number,3), V, mnum(3,:));
Gm1 = V.(mnum(3,:));

% Gp1
V = readVTKast(filename(Number,4), V, mnum(4,:));
Gp1 = V.(mnum(4,:));

% Gm2
V = readVTKast(filename(Number,5), V, mnum(5,:));
Gm2 = V.(mnum(5,:));

% Gp2
V = readVTKast(filename(Number,6), V, mnum(6,:));
Gp2 = V.(mnum(6,:));

%% Extracting Variables from  Data

delB     = kappa.*(sigma.*Gp0 - (1+tau/Z).*Gm0);            % Magnetic Field Fluctuations
deln     = kappa.*(sigma.*Gm0 - (tau/Z).*(2/beta).*Gp0);    % Density Fluctuations
delU     = kappa.*(vth/2).*((1+(Z/tau)+sigma)*Gm1 - (sigma + (tau/Z)*(2/beta)).*Gp1);   % Parallel Velocity Fluctuations
delPpara = kappa.*((1+(Z/tau)+sigma).*(sqrt(2).*Gm2 + Gm0) - (sigma + (tau/Z)*(2/beta)).*(sqrt(2).*Gp2 + Gp0)); %without delu terms
% delPpara = kappa.*((1+(Z/tau)+sigma).*(sqrt(2).*Gm2 - ... % Includes delU terms
%     (2*sqrt(2)/vth).*delU.*Gm1 + (1+(2*(delU).^2)/vth^2).*Gm0) - ...
%     (sigma + (tau/Z)*(2/beta)).*(sqrt(2).*Gp2 - ...
%     (2*sqrt(2)/vth).*delU.*Gp1 + (1+(2*(delU).^2)/vth^2).*Gp0));
delPperp = -(1/beta)*delB - (Z/(2*tau))*deln;               % Perpendicular Pressure Fluctuations
delP     = (1/3)*(2.*delPperp + delPpara);                  % Total Pressure

%% Linear Interpolation
lstart = [0.4 0.5 0.5];
ldir = [-0.1 sqrt(0.2) (pi^2)]; 
L = [LX, LY, LZ];
[delU_line, line_length] = Interpolate(delU, dx, LX, lstart, ip, jp, kp, L, ldir);
[delB_line,     ~] = Interpolate(delB, dx, LX, lstart, ip, jp, kp, L, ldir);
[deln_line,     ~] = Interpolate(deln, dx, LX, lstart, ip, jp, kp, L, ldir);
[delPperp_line, ~] = Interpolate(delPperp, dx, LX, lstart, ip, jp, kp, L, ldir);
[delPpara_line, ~] = Interpolate(delPpara, dx, LX, lstart, ip, jp, kp, L, ldir);
[delP_line,     ~] = Interpolate(delP, dx, LX, lstart, ip, jp, kp, L, ldir);
% 
% delB_line     = double(delB_line);      % cwt requires data to be in double format
% deln_line     = double(deln_line);
% delPperp_line = double(delPperp_line);

%% Wavelet Transform

% Plot of Wavelet Coherence
% wcoherence(delU_line,delB_line, seconds(dx), 'PhaseDisplayThreshold', 0.5)

[wcohXi , Xi_nB , period, coi] = wcoherence(deln_line,delB_line, seconds(dx));
[wcohChi, Chi_UB] = wcoherence(delU_line,delB_line, seconds(dx));
[wcohApe, alp_pe] = wcoherence(delPperp_line,delB_line, seconds(dx));
[wcohApa, alp_pa] = wcoherence(delPpara_line,delB_line, seconds(dx));
[wcohPsi, Psi_PB] = wcoherence(delP_line,delB_line, seconds(dx));

% Trim off scales larger than LX/2 to avoid measuring coherence of 'same point'
period = seconds(period);                 % Tidy this shit up
smallscale  = find(period<(LX/2));        % Only look at length scales shorter than 1/2 box-size to avoid cross-correlation
wcsUBsmall  = wcsUB(smallscale, :);
wcohUBsmall = wcohUB(smallscale, :);
periodsmall = period(smallscale);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolation Function
function [out, lvec] = Interpolate(delA, dx, LX, lstart, i, j, k, L, ldir)

ldir = ldir/norm(ldir); % Direction along which to look
dl = dx;
total_length = 100*LX^2 ;
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
