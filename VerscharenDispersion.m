%% Load Data

va    = 1;                      %input.Parameters.va;
beta_An  = logspace(-2,2,1000);    %input.Parameters.beta;
gamma = 5/3;
theta = 3*pi/7;

% s_plus  = output.sp;
% s_minus = output.sm;

%% Calculate Perturbations

% delB =  (1/(2*va))*(1+(1/beta)^2)^(-0.5)*(s_plus - s_minus);  % deltaB_parallel/B_0
% delU =  (1/(2*va))*(s_plus + s_minus);                        % deltaU_parallel/v_A
% delR = -(1/beta^2)*delB;                                      % Density
% delP = gamma*delR;                                            % Pressure

%% Find Verscharen Coefficient Magnitudes

% XiMag_grid  = abs(delR)./abs(delB);
% ChiMag_grid = abs(delU)./abs(delB);
% PsiMag_grid = abs(delP)./abs(delB);

%% Find Verscharen Coefficient Phases
% Calculate Wavelet transform of delB, delU, delR, delP for a bunch of
% scales a and positions b

% Calculate Coherence spectrum C for each quantity

% Calculate phase of coefficients by phase = atan(Im(C)/Real(C))

%% Analytic Dispersions

% !!!!!!! Need a theta? Where do I get this from?! !!!!!!!!
C_minus = sqrt(0.5*(1+(gamma*beta_An)/2) - 0.5*sqrt((1+(gamma*beta_An)/2).^2 - 2*gamma*beta_An.*(cos(theta)).^2));

Xi_An  = C_minus.^2./(C_minus.^2 - 0.5.*gamma.*beta_An.*(cos(theta))^2);
Chi_An = C_minus.*cos(theta).*0.5.*gamma.*beta_An./(C_minus.^2 - 0.5.*gamma.*beta_An.*(cos(theta))^2);
Psi_An = gamma.*beta_An.*Xi_An;

%% Plotting
figure(1)

%%% Xi  - Density
% Magnitude
subplot(2,3,1)
loglog(beta_An, abs(Xi_An), '--')    % Analytic Dispersion
% Phase
subplot(2,3,4)

%%% Chi - Velocity
% Magnitude
subplot(2,3,2)
semilogx(beta_An, abs(Chi_An), '--') % Analytic Dispersion
% Phase

%%% Psi - Pressure
% Magnitude
subplot(2,3,3)
semilogx(beta_An, abs(Psi_An), '--') % Analytic Dispersion
axis([10^(-2) 10^2 0 10])
%Phase



% For each coefficient Xi, Chi and Psi we need:
% Plot of coefficient magnitude vs beta
% Plot of coefficient phase vs beta