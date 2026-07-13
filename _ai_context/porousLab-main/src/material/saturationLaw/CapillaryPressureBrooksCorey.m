%% CapillaryPressureBrooksCorey Class
% This class implements the Brooks-Corey capillary pressure model, which 
% is used to compute the saturation degree and its derivative for a porous
% medium. It inherits from the _CapillaryPressure_ base class.
%
% 
%% Methods
% * *saturationDegree*: Computes the liquid phase saturation degree based 
%                       on the capillary pressure and the properties of 
%                       the porous medium. 
% * *derivativeSaturationDegree*: Computes the derivative of the liquid 
%                                 phase saturation degree with respect to 
%                                 the capillary pressure.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef CapillaryPressureBrooksCorey < CapillaryPressure  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = CapillaryPressureBrooksCorey()
            this = this@CapillaryPressure('brooksCorey');
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability
        function Sl = saturationDegree(~, pc, porousMedia)
            if (pc >= porousMedia.Pb)
                Se = (porousMedia.Pb/pc)^porousMedia.lambda;
                Sl = Se * (1.0 - porousMedia.Slr - porousMedia.Sgr) + porousMedia.Slr;
            else
                Sl = 1.0;
            end
            Sl = max(min(Sl, 1.0 - porousMedia.Sgr - eps), porousMedia.Slr + eps);
        end
        
        %------------------------------------------------------------------
        % Compute the gas phase relative permeability
        function dSldpc = derivativeSaturationDegree(this, pc, porousMedia)
            % if (pc >= porousMedia.Pb)
            % 
            % else
            %     dSldpc = 1.0;
            % end
            Sl = this.saturationDegree(pc, porousMedia);
            dPcdSl = (porousMedia.Pb / (porousMedia.lambda * (porousMedia.Slr - Sl))) * ((Sl - porousMedia.Slr)/(1.0 - porousMedia.Sgr - porousMedia.Slr)) ^ (-1.0/porousMedia.lambda);
            dSldpc = 1.0 / dPcdSl;
        end
        
    end
end