%% RelativePermeabilityBrooksCoreyLiquid Class
% This class implements the Brooks-Corey model for calculating the relative 
% permeability of the liquid phase in a porous medium. It inherits from the 
% _RelativePermeability_ base class.
%
%% Method
% * *calculate*: Computes the liquid phase relative permeability based on 
%                the effective saturation degree Se and the Brooks-Corey 
%                model.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef RelativePermeabilityBrooksCoreyLiquid < RelativePermeability  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeabilityBrooksCoreyLiquid()
            this = this@RelativePermeability('brooksCoreyLiquid');
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability
        function klr = calculate(~, Sl, porousMedia)
            Se = porousMedia.effectiveSaturationDegree(Sl);
            klr = Se^((2.0 + 3.0 * porousMedia.lambda)/porousMedia.lambda);
            klr = max(klr,porousMedia.klrmin);
        end

        %------------------------------------------------------------------
        % Compute the derivative of the liq phase relative permeability wrt
        % the liquid saturation degree
        function dkrdSl = derivative(~, Sl, porousMedia)
            lambda = porousMedia.lambda;
            Se = porousMedia.effectiveSaturationDegree(Sl);
            dSedSl = porousMedia.derivativeEffectiveSaturationDegree();
            dkrdSe = (Se^((2*(lambda + 1))/lambda)*(3*lambda + 2))/lambda;
            dkrdSl = dkrdSe * dSedSl;
        end
        
    end
end