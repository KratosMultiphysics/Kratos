%% RelativePermeabilityLiakopoulosLiquid Class
% This class defines the liquid phase relative permeability model based on
% the Liakopoulos formulation. It inherits from the _RelativePermeability_
% base class and provides a specific implementation for calculating the
% relative permeability of the liquid phase.
%
%% Method
% * *calculate*: Computes the liquid phase relative permeability klr based 
%                on the liquid saturation Sl and the porous media 
%                properties.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef RelativePermeabilityLiakopoulosLiquid < RelativePermeability  
    %% Properties
    % Parameters taken from OGS-6. 
    % OGS reference:
    % Asadi, R., Ataie-Ashtiani, B. (2015): A Comparison of finite volume
    % formulations and coupling strategies for two-phase flow in deforming
    % porous media. Comput. Geosci., p. 24ff.
    properties (SetAccess = public, GetAccess = public)
        a     = 2.207;
        b     = 1.0121;
        Slmin = 0.2;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeabilityLiakopoulosLiquid()
            this = this@RelativePermeability('liakopoulosLiquid');
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability
        function klr = calculate(this, Sl, porousMedia)
            if (Sl < this.Slmin)
                klr = porousMedia.klrmin;
            elseif (Sl > 1.0)
                klr = 1.0;
            else
                klr = 1.0 - this.a * (1.0 - Sl)^this.b;
            end
        end

        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability
        function dkrdSl = derivative(this, Sl, ~)
            if (Sl < this.Slmin)
                dkrdSl = 0.0;
            elseif (Sl > 1.0)
                dkrdSl = 0.0;
            else
                Sl = min(max(Sl,this.Slmin),1.0);
                dkrdSl = this.a * this.b * (1.0 - Sl)^(this.b - 1.0);
            end
        end
        
    end
end