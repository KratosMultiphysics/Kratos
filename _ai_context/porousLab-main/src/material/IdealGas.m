%% IdealGas Class
% This class is a subclass of the _Fluid_ class that represents an ideal 
% gas. It provides methods to calculate the density of the gas based on 
% the ideal gas law and allows setting key properties such as the 
% universal gas constant, molar mass, and temperature.
%
%% Methods
% * *getDensity*: Calculates the density of the gas using the ideal gas
%                 law considering the gas pressure
% * *setUniversalGasConstant*: Sets the value of the universal gas
%                              constant.
% * *setMolarMass*: Sets the molar mass of the gas.
% * *setTemperature*: Sets the temperature of the gas.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef IdealGas < Fluid   
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        R = 8.3144621;      % Universal gas constant (J/(mol*K)
        T = 293.15;         % Temperature (K)
        M = 0.02897;        % Molar mass (kg/mol)
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = IdealGas(id, viscosity)
            this = this@Fluid(id);
            this.rho = 0.0;
            if nargin > 1
                this.mu = viscosity;
                this.K  = this.M / (this.R * this.T);
            end
        end
    end
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Get the fluid density based on the ideal Gas law
        function rho = getDensity(this,pg)
            rho = pg * this.M / (this.R * this.T);
        end
        %------------------------------------------------------------------
        % Get the fluid bulk modulus based on the ideal Gas law
        function K = getBulkModulus(this,~)
            K = this.M / (this.R * this.T);
        end
        %------------------------------------------------------------------
        % Set the value of the universal gas constant
        function setUniversalGasConstant(this,R)
            this.R = R;
        end
        %------------------------------------------------------------------
        % Set the value of the gas molar mass
        function setMolarMass(this,M)
            this.M = M;
        end
        %------------------------------------------------------------------
        % Set the value of the temperature
        function setTemperature(this,T)
            this.T = T;
        end
    end
end