%% RelativePermeabilityUMAT Class
% This class defines a relative permeability model using a user-defined 
% material (UMAT) approach. It inherits from the _RelativePermeability_ 
% base class and provides functionality to compute relative permeability 
% based on saturation curves.
%
%% Method
% * *calculate*: Computes the relative permeability kr for a given 
%                saturation Sl using linear interpolation. The result is 
%                clamped between the minimum relative permeability klrmin
%                defined in the _porousMedia_ object and 1.0.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef RelativePermeabilityUMAT < RelativePermeability  
     %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        Sl_curve = [];          % Must be sorted!!
        kr_curve = [];          % Must be sorted!!
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RelativePermeabilityUMAT(Sl_curve,kr_curve)
            this = this@RelativePermeability('umat');
            this.Sl_curve = Sl_curve;
            this.kr_curve = kr_curve;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the relative permeability
        function kr = calculate(this, Sl, porousMedia)
            kr = interp1(this.Sl_curve,this.kr_curve,Sl,'linear', 'extrap');
            kr = max(min(kr, 1.0), porousMedia.klrmin);
        end
        
    end
end