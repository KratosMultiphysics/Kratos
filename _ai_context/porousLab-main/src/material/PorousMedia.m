%% PorousMedia Class
% This class defines a porous media object with various physical and 
% mechanical properties. It includes methods for setting and retrieving 
% these properties, as well as performing calculations related to porous 
% media behavior.
%
%% Methods
% * *effectiveSaturationDegree*: Computes the effective saturation degree.
% * *intrinsicPermeabilityMatrix*: Returns the intrinsic permeability 
%                                  matrix.
% * *setMinLiquidRelPermeability*: Sets the minimum liquid relative 
%                                  permeability.
% * *setMinGasRelPermeability*: Sets the minimum gas relative 
%                               permeability.
% * *setUMATCapillaryPressureCurve*: Sets the user-defined capillary 
%                                    pressure curve.
% * *setUMATLiquidRelPermCurve*: Sets the user-defined liquid relative 
%                                permeability curve.
% * *setUMATGasRelPermCurve*: Sets the user-defined gas relative 
%                             permeability curve.
% * *setMechanicalConstitutiveLaw*: Sets the mechanical constitutive law.
% * *setMechanicalProperties*: Sets the mechanical properties (Young's 
%                              modulus and Poisson's ratio).
% * *setDensity*: Sets the density.
% * *getDensity*: Retrieves the density.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef PorousMedia < handle & matlab.mixin.Copyable    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id                   = '';
        Young                = [];              % Young modulus (Pa)
        nu                   = [];              % Poisson ratio
        sy0                  = [];              % Initial yield stress (Pa)
        Kp                   = [];              % Plastic modulus (Pa)
        cohesion             = [];              % Cohesion (Pa)
        frictionAngle        = [];              % Friction angle (rad)
        dilationAngle        = [];              % Dilation angle (rad)
        stressIntAlgorithm   = 'implicit';      % Stress integration algorithm
        MCmatch              = 'planestrain';   % How the DP surfaces match the MC ones
        friction             = [];              % Friction coefficient
        asympt               = [];              % Asymptotic model
        eref                 = [];              % Reference strain for asymptotic model
        sy                   = [];              % Isotropic tensile limit (Pa)
        tauy                 = [];              % Shear yield stress (Pa)
        kappa                = [];              % Ratio between the uniaxial compressive strength and the uniaxial tensile strength
        DamageThreshold      = [];              % Damage threshold
        FractureEnergyMode1  = [];              % Fracture energy associated with mode 1 (N/m)
        rho                  = [];              % Density (kg/m3)
        K                    = 0.0;             % Intrinsic permeability (m2)  
        phi                  = 0.0;             % Porosity
        biot                 = 1.0;             % Biot's coefficient
        Ks                   = 1.0e25;          % Solid bulk modulus (Pa)
        Slr                  = 0.0;             % Residual liquid saturation
        Sgr                  = 0.0;             % Residual gas saturation 
        Pb                   = 0.0;             % Gas-entry pressure (Pa)
        lambda               = 0.0;             % Curve-fitting parameter
        liqRelPermeability   = 'BrooksCorey';   % Liquid relative permeability
        gasRelPermeability   = 'BrooksCorey';   % Gas relative permeability
        capillaryPressure    = 'BrooksCorey';   % Saturation degree function
        mechanical           = 'elastic';       % Mechanical constitutive law
        b                    = [0.0;-1.0];      % Gravity force direction vector       
        SlPc_umat            = [];              % User material curve saturation law
        klr_umat             = [];              % User material curve liquid relative permeability
        kgr_umat             = [];              % User material curve gas relative permeability
        m                    = 1;               % Exponent for the polynomial relationships
    end
    properties (SetAccess = protected, GetAccess = public)
        klrmin               = 1.0e-9;          % Minimum liquid relative permeability
        kgrmin               = 1.0e-9;          % Minimum gas relative permeability
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = PorousMedia(id, permeability, porosity, ...
                biotCoefficient, solidBulkModulus, ...
                residualLiquidSaturationDegree, ...
                residualGasSaturationDegree,...
                gasEntryPressure,curveFittingParameter, ...
                liqRelPermeability, gasRelPermeability, ...
                capillaryPressure)
            if nargin == 1, this.id = id; end
            if nargin > 1
                this.id                   = id;
                this.K                    = permeability;
                this.phi                  = porosity;
                this.biot                 = biotCoefficient;
                this.Ks                   = solidBulkModulus;
                this.Slr                  = residualLiquidSaturationDegree;
                this.Sgr                  = residualGasSaturationDegree;
                this.Pb                   = gasEntryPressure;
                this.lambda               = curveFittingParameter;
                this.liqRelPermeability   = liqRelPermeability;
                this.gasRelPermeability   = gasRelPermeability;
                this.capillaryPressure    = capillaryPressure;
            end
        end
    end
    %% Public methods
    methods
        % -----------------------------------------------------------------
        % Compute the effective saturation degree
        function Se = effectiveSaturationDegree(this,Sl)
            Se = (Sl - this.Slr)/(1.0 - this.Slr - this.Sgr);
        end

        % -----------------------------------------------------------------
        % Compute the derivative of the effective saturation degree wrt the
        % liquid saturation degree
        function dSedSl = derivativeEffectiveSaturationDegree(this)
            dSedSl = 1.0/(1.0 - this.Slr - this.Sgr);
        end

        % -----------------------------------------------------------------
        % Create the intrinsic permeability matrix
        function Km = intrinsicPermeabilityMatrix(this)
            Km = this.K * eye(2);
        end

        % -----------------------------------------------------------------
        % Set the minimum liquid relative permeability
        function setMinLiquidRelPermeability(this,klrmin)
            this.klrmin = klrmin;
        end

        % -----------------------------------------------------------------
        % Set the minimum gas relative permeability
        function setMinGasRelPermeability(this,kgrmin)
            this.kgrmin = kgrmin;
        end

        % -----------------------------------------------------------------
        % Set the user defined capillary pressure curve
        function setUMATCapillaryPressureCurve(this,curve)
            curve = sortrows(curve, 2);
            this.SlPc_umat = curve;
        end

        % -----------------------------------------------------------------
        % Set the user defined liquid relative permeability curve
        function setUMATLiquidRelPermCurve(this,curve)
            curve = sortrows(curve, 1);
            this.klr_umat = curve;
        end

        % -----------------------------------------------------------------
        % Set the user defined gas relative permeability curve
        function setUMATGasRelPermCurve(this,curve)
            curve = sortrows(curve, 1);
            this.kgr_umat = curve;
        end

        % -----------------------------------------------------------------
        % Set the mecahnical constitutive law
        function setMechanicalConstitutiveLaw(this,law)
            this.mechanical = law;
        end

        % -----------------------------------------------------------------
        % Set the mechanical properties of the material
        function setMechanicalProperties(this,E,nu)
            this.Young = E;
            this.nu = nu;
        end

        % -----------------------------------------------------------------
        % Set the density of the material
        function setDensity(this,rho)
            this.rho = rho;
        end

        % -----------------------------------------------------------------
        % Get the density of the material
        function rho = getDensity(this)
            rho = this.rho;
        end
    end
end
