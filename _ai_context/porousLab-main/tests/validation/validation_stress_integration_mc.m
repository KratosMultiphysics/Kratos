%% DESCRIPTION
%
% Validation script for mechanical constitutive models.
%
% Physics:
% * Mechanical (M)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
close all; clear; clc;

%% MODEL DEFINITION

% Material properties
rock = PorousMedia('rock');
rock.mechanical    = 'mohrCoulomb';  % Constitutive law
rock.Young         = 2.0e+4;           % Young modulus (kPa)
rock.nu            = 0.49;             % Poisson ratio
rock.cohesion      = 50.0;             % Cohesion (kPa)
rock.frictionAngle = 20.0*pi/180;      % Friction angle (rad)
rock.dilationAngle = 20.0*pi/180;      % Dilation angle
rock.stressIntAlgorithm = 'alternative';

% Material parameters vector
mat = struct('porousMedia', rock);

% Integration point
ip = IntPoint([0.0,0.0], 1.0,  Material_M(mat));
ip.initializeMechanicalAnalysisModel('PlaneStrain');

%% PLOT MC ENVELOPE

% Compression negative (engineering sign). Envelope: tau = c - sigma_n*tan(phi)
sig_n = linspace(-10*rock.cohesion, rock.cohesion/tan(rock.frictionAngle), 400);
tau_env = rock.cohesion - sig_n*tan(rock.frictionAngle);
tau_env(tau_env<0) = NaN;                % no negative shear strength

figure; grid on; box on; hold on;
plot(sig_n,  tau_env,     'k-', 'LineWidth',1.8, 'DisplayName','MC envelope (+)');
plot(sig_n, -tau_env,     'k-', 'LineWidth',1.8, 'HandleVisibility','off'); % lower branch
xlabel('$\sigma_n\;[\mathrm{kPa}]$','Interpreter','latex');
ylabel('$\tau\;[\mathrm{kPa}]$','Interpreter','latex');
legend('Location','best');
set(gca, 'FontSize',18, 'TickLabelInterpreter','latex');
axis equal

%% RUN ANALYSIS

% Initial state
ip.strainOld = [-0.0001; -0.0001; 0.0; 0.0];

% Strain increment direction
dstrain_inc = [0.0; 0.0; 0.0; 0.1];

mag  = 1.0e-02;
ninc = 50;

p  = zeros(ninc+1,1);
q  = zeros(ninc+1,1);

% NEW: path in (sigma_n, tau)
sigma_n_path = zeros(ninc+1,1);
tau_path     = zeros(ninc+1,1);

for i = 1:ninc
    ip.strain = ip.strainOld + mag * dstrain_inc;
    ip.stress = ip.mechanicalLaw();
    ip.stressOld = ip.stress;
    ip.strainOld = ip.strain;

    % Invariants
    p(i+1) = ip.constitutiveMdl.mechanical.hydrostaticStress(ip.stress);
    q(i+1) = ip.constitutiveMdl.mechanical.vonMisesStress(ip.stress);

    % NEW: principal stresses, then (sigma_n, tau) on plane of max shear
    sig = [ip.stress(1) ip.stress(4) 0;
           ip.stress(4) ip.stress(2) 0;
           0            0            ip.stress(3)];
    ev  = eig((sig+sig.')/2);           % symmetric
    s1  = max(ev);                       % largest principal
    s3  = min(ev);                       % smallest principal
    plotCircle(abs(s1-s3)/2,(s1+s3)/2)
end

%% AUXILIARY FUNCTION
function plotCircle(r,c)
    npts = 100;
    angle = linspace(0.0, 2*pi, npts);
    circle = zeros(npts,2);
    for i = 1 : npts
        circle(i,1) = c + r * cos(angle(i));
        circle(i,2) = r * sin(angle(i));
    end
    plot(circle(:,1),circle(:,2),'-r','HandleVisibility','off');
end