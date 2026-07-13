%% DESCRIPTION
%
% Terzaghi 1D consolidation problem.
%
% References:
% * Nguyen et al (2017). Modelling hydraulic fractures in porous media using flow cohesive interface elements. Eng Geol, 225:68–82.
% * Segura and Carol (2008). Coupled HM analysis using zero-thickness interface elements with double nodes—Part II: Verification and application. Int J Numer Anal Methods Geomech, 32(18):2103–23.
%
% Physics:
% * Single-phase flow hydro-mechanical (HM)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
clear; clc; % close all

%% PHYSICAL PARAMETERS

% Solid parameters
E  = 1.0e+6;       % Young's modulus (Pa)
nu = 0.3;          % Poisson's ratio (-)
n  = 0.3;          % Porosity (-)
k  = 1.15741e-12;  % Intrinsic permeability (m^2)
Ks = 1.0e+25;      % Solid bulk modulus (Pa)

% Fluid parameters
muf = 1e-3;  % Fluid dynamic viscosity (Pa*s)

% Porous medium parameters
alpha = 1.0;                % Biot coefficient
Sm    = 0.0;                % Storage coefficient
K     = E/3/(1-2*nu);       % Bulk modulus (Pa)
G     = E/2/(1+nu);         % Shear modulus (Pa)
mv    = 1 / ((4/3)*G + K);  % Confined compressibility (1/Pa)

% Consolidation coefficient
Cv = k / (muf * (Sm + alpha*alpha * mv));    

%% PROBLEM PARAMETERS

% Initial excess pore pressure (Pa)
p0 = 1.0e4;

% Total thickness of soil layer (m)
H = 1;

% Consolidation time (s)
t = 100.0;
% t = [43, 100, 262, 302.5, 500];

% Number of expansion terms
np = 1000;

% Depth intervals
dz = 0.001;   % Depth increment (m)
z  = 0:dz:H;  % Depths within soil layer
Z  = z/H;     % Normalized depth

%% TERZAGHI SOLUTION

% Initialize plot
figure;
hold on;

% Calculate pore pressure profile
for i = 1:length(t)
    % Initialize pressure vector
    p = zeros(size(z));

    % Calculate time factor
    T = Cv * t(i) / (H^2);

    % Terzaghi's solution
    for j = 1:np
        M = pi/2*(2*j-1);
        p = p + 4/pi * ((-1)^(j-1)/(2*j-1)) * cos(M*Z) * exp(-M*M*T);
    end

    % Plot pressure vector
    plot(p*p0, z/H, 'LineWidth', 2, 'DisplayName', ['Analytical: t = ',num2str(t(i)), ' s']);

    % Print pressure at bottom
    disp(['Pressure at the bottom edge: ', num2str(p0*p(1)) , ' Pa']);
end

% Formatting
xlabel('Pore Pressure (Pa)');
ylabel('Y (m)');
set(gca, 'FontSize', 18)
grid on;
box on;
legend('Location', 'best');
