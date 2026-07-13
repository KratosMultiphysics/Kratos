%% RelativePermeabilityPolynomialGas Class
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
classdef RelativePermeabilityPolynomialGas < RelativePermeability  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeabilityPolynomialGas()
            this = this@RelativePermeability('polynomialGas');
        end
    end

    %% Public methods
    methods
        
        %------------------------------------------------------------------
        % Compute the gas phase relative permeability
        function kgr = calculate(~, Sl, porousMedia)
            kgr = (1.0 - Sl)*porousMedia.m;
            kgr = max(kgr,porousMedia.kgrmin);
        end

        %------------------------------------------------------------------
        % Compute the gas phase relative permeability derivative wrt the
        % liquid saturation degree
        function dkrdSl = derivative(this, Sl, porousMedia)
            kr = this.calculate(Sl, porousMedia);
            if (kr==porousMedia.kgrmin)
                dkrdSl = 0.0;
            else
                dkrdSl = -porousMedia.m;
            end
        end
        
    end
end