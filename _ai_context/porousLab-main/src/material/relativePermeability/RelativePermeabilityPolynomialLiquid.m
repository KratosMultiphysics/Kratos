%% RelativePermeabilityPolynomialLiquid Class
% This class implements the relative permeability calculation for the gas 
% phase using a polynomial model. It inherits from the 
% _RelativePermeability_ base class and provides a specific implementation 
% for gas phase relative permeability.
% 
%% Method
% * *calculate*: Computes the gas phase relative permeability based on the 
%                liquid saturation Sl and the properties of the porous 
%                medium _porousMedia_. The calculation ensures that the 
%                relative permeability does not fall below a specified 
%                minimum value kgrmin.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef RelativePermeabilityPolynomialLiquid < RelativePermeability  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeabilityPolynomialLiquid()
            this = this@RelativePermeability('polynomialLiquid');
        end
    end

    %% Public methods
    methods
        
        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability
        function klr = calculate(~, Sl, porousMedia)
            klr = Sl*porousMedia.m;
            klr = max(klr,porousMedia.klrmin);
        end

        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability derivative wrt the
        % liquid saturation degree
        function dkrdSl = derivative(this, Sl, porousMedia)
            kr = this.calculate(Sl, porousMedia);
            if (kr==porousMedia.klrmin)
                dkrdSl = 0.0;
            else
                dkrdSl = porousMedia.m;
            end
        end
        
    end
end