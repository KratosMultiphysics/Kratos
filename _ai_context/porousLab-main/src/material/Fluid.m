%% Fluid Class
% This class defines a fluid object with properties such as density, 
% viscosity, and compressibility. It provides methods to access these 
% properties and allows initialization with specific values.
%
%% Methods
% * *Fluid*: Constructor to initialize the fluid object. If only the id is 
%           provided, other properties are set to their default values.
% * *getDensity*: Returns the density of the fluid.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef Fluid < handle  
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id   = '';
        rho  = 1000.0;      % Density (kg/m3)
        mu   = 1.0e-3;      % Viscosity (Pa*s)
        K    = 1.0e25;      % Compressibility/Bulk modulus (1/Pa)
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        % Initialization of the fluid
        function this = Fluid(id, density, viscosity, compressibility)
            if nargin == 1
                this.id = id; 
            elseif nargin > 1
                this.id   = id;
                this.rho  = density;
                this.mu   = viscosity;
                this.K    = compressibility;
            end
        end
    end
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Get the density of the fluid
        function rho = getDensity(this,~)
            rho = this.rho;
        end
        %------------------------------------------------------------------
        % Get the fluid bulk modulus based on the ideal Gas law
        function K = getBulkModulus(this,~)
            K = this.K;
        end
    end
end