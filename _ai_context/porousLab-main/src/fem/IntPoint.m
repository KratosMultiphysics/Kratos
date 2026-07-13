%% IntPoint Class
% This is an abstract class that defines an integration point object.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
% 
%% Class definition
classdef IntPoint < handle    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        X                = [];   % Coordinates of the integration point in the natural coordinate system
        w                = 0.0;  % Weight associated to the integration point
        strain           = [];   % Current strain vector
        stress           = [];   % Current stress vector
        plasticstrain    = [];   % Current plastic strain vector
        statevar         = [];   % Current state variables vector
        strainOld        = [];   % Previous strain vector
        stressOld        = [];   % Previous stress vector  
        plasticstrainOld = [];   % Current plastic strain vector
        statevarOld      = [];   % Previous state variables vector
        constitutiveMdl  = [];   % Constitutive model object
        anm              = '';   % Analysis model tag
        nVar             = 4;    % Dimension of the stress and strain vectors
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = IntPoint(X,w,constitutiveMdl)
            if nargin > 0
                this.X = X;
                this.w = w;
                this.constitutiveMdl = constitutiveMdl;
            end
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Initialize analysis model (mechanical part).
        function initializeMechanicalAnalysisModel(this,anm)
            this.anm = anm;

            if strcmp(anm,'PlaneStress')
                this.nVar = 4;
            elseif strcmp(anm,'PlaneStrain')
                this.nVar = 4;
            elseif strcmp(anm,'AxisSymmetrical')
                this.nVar = 4;
            elseif strcmp(anm,'Interface')
                this.nVar = 2;
            end

            nStvar = this.constitutiveMdl.getNumberStateVar();
            this.strain      = zeros(this.nVar, 1);
            this.stress      = zeros(this.nVar, 1);
            this.statevar    = zeros(nStvar,    1);
            this.strainOld   = zeros(this.nVar, 1);
            this.stressOld   = zeros(this.nVar, 1);
            this.statevarOld = zeros(nStvar,    1);

            if this.constitutiveMdl.hasPlasticStrain()
                this.plasticstrain    = zeros(this.nVar,1);
                this.plasticstrainOld = zeros(this.nVar,1);
            end
        end

        %------------------------------------------------------------------
        % Update current strain vector.
        function updateStrainVct(this)
            this.strainOld        = this.strain;
            this.plasticstrainOld = this.plasticstrain;
        end

        %------------------------------------------------------------------
        % Update current state variable vector.
        function updateStateVar(this)
            this.statevarOld = this.statevar;
        end

        %------------------------------------------------------------------
        % Update current stress vector.
        function updateStressVct(this)
            this.stressOld = this.stress;
        end

        %------------------------------------------------------------------
        % Reset the integration point.
        function reset(this)
            nStvar = this.constitutiveMdl.getNumberStateVar();
            this.strain      = zeros(this.nVar, 1);
            this.stress      = zeros(this.nVar, 1);
            this.statevar    = zeros(nStvar,    1);
            this.strainOld   = zeros(this.nVar, 1);
            this.stressOld   = zeros(this.nVar, 1);
            this.statevarOld = zeros(nStvar,    1);
            if this.constitutiveMdl.hasPlasticStrain()
                this.plasticstrain    = zeros(this.nVar,1);
                this.plasticstrainOld = zeros(this.nVar,1);
            end
        end

        %------------------------------------------------------------------
        % Get current constitutive matrix.
        function D = getConstitutiveMtrx(this,dStrain)
            D = this.constitutiveMdl.constitutiveMtrx(dStrain,this);
        end

        %------------------------------------------------------------------
        % Compute the stress and constitutive matrix using the constitutive model
        % and update the current stress vector.
        function [stress,D] = mechanicalLaw(this)
            [stress,D] = this.constitutiveMdl.mechanicalLaw(this);
            this.stress = stress;
        end
    end
end
