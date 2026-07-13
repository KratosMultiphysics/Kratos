%% Material_H2 class
% This class defines the material properties and behavior for a porous 
% medium involving liquid and gas phases. It includes methods for 
% computing saturation degrees, saturation derivatives, relative
% permeabilities, relative permeability derivatives, and gas-density
% derivatives for the liquid and gas phases.
%
%% Methods
% * *saturationDegree*: Computes the liquid saturation degree based on 
%                       the capillary pressure pc.
% * *derivativeSaturationDegree*: Computes the derivative of the liquid
%                                 saturation degree with respect to the
%                                 capillary pressure pc.
% * *relativePermeabilities*: Computes the liquid and gas relative
%                             permeabilities.
% * *derivativeRelPerm*: Computes the derivatives of the liquid and gas
%                        relative permeabilities with respect to liquid
%                        saturation.
% * *derivativeGasDensityWrtGasPressure*: Computes the derivative of the 
%                                         gas density with respect to the 
%                                         gas pressure Pg.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef Material_H2 < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        liqRelativePermeability = [];
        gasRelativePermeability = [];
        capillaryPressure       = [];
        liquidFluid             = Fluid();
        gasFluid                = Fluid();
        porousMedia             = PorousMedia();
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Material_H2(matData)
            this.liquidFluid = matData.liquidFluid;
            this.gasFluid    = matData.gasFluid;
            this.porousMedia = matData.porousMedia;
            % Liquid phase relative permeability function
            if strcmp('BrooksCorey',matData.porousMedia.liqRelPermeability)
                this.liqRelativePermeability = RelativePermeabilityBrooksCoreyLiquid();
            elseif strcmp('Liakopoulos',matData.porousMedia.liqRelPermeability)
                this.liqRelativePermeability = RelativePermeabilityLiakopoulosLiquid();
            elseif strcmp('PolynomialLiquid',matData.porousMedia.liqRelPermeability)
                this.liqRelativePermeability = RelativePermeabilityPolynomialLiquid();
            elseif strcmp('UMAT',matData.porousMedia.liqRelPermeability)
                curve = matData.porousMedia.klr_umat;
                this.liqRelativePermeability = RelativePermeabilityUMAT(curve(:,1),curve(:,2));
            end
            % Gas phase relative permeability function
            if strcmp('BrooksCorey',matData.porousMedia.gasRelPermeability)
                this.gasRelativePermeability = RelativePermeabilityBrooksCoreyGas();
            elseif strcmp('PolynomialGas',matData.porousMedia.gasRelPermeability)
                this.gasRelativePermeability = RelativePermeabilityPolynomialGas();
            elseif strcmp('UMAT',matData.porousMedia.gasRelPermeability)
                curve = matData.porousMedia.kgr_umat;
                this.gasRelativePermeability = RelativePermeabilityUMAT(curve(:,1),curve(:,2));
            end
            % Saturation degree function
            if strcmp('BrooksCorey',matData.porousMedia.capillaryPressure)
                this.capillaryPressure = CapillaryPressureBrooksCorey();
            elseif strcmp('UMAT',matData.porousMedia.capillaryPressure)
                curve = matData.porousMedia.SlPc_umat;
                this.capillaryPressure = CapillaryPressureUMAT(curve(:,1),curve(:,2));
            elseif strcmp('Liakopoulos',matData.porousMedia.capillaryPressure)
                this.capillaryPressure = CapillaryPressureLiakopoulos();
            elseif strcmp('Log',matData.porousMedia.capillaryPressure)
                this.capillaryPressure = CapillaryPressureLog();
            end
        end
    end
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Get the liquid saturation degree
        function Sl = saturationDegree(this,pc)
            Sl = this.capillaryPressure.saturationDegree(pc, this.porousMedia);
            if isnan(Sl) || isinf(Sl)
                error("Liquid saturation degree is NaN or Inf.");
            end
        end

        %------------------------------------------------------------------
        % Get the liquid saturation degree derivative wrt the capillary
        % pressure
        function dSldpc = derivativeSaturationDegree(this,pc)
            dSldpc = this.capillaryPressure.derivativeSaturationDegree(pc, this.porousMedia);
            if isnan(dSldpc) || isinf(dSldpc)
                error("Derivative of the liquid saturation degree is NaN or Inf.");
            end
        end

        %------------------------------------------------------------------
        % Compute the relative permeabilities
        function [klr, kgr] = relativePermeabilities(this,Sl)
            klr = this.liqRelativePermeability.calculate(Sl, this.porousMedia);
            if isnan(klr) || isinf(klr)
                error("Liquid relative permeability is NaN or Inf.");
            end
            kgr = this.gasRelativePermeability.calculate(Sl, this.porousMedia);
            if isnan(kgr) || isinf(kgr)
                error("Gas relative permeability is NaN or Inf.");
            end
        end

        %------------------------------------------------------------------
        % Compute the relative permeabilities
        function [dklrdSl, dkgrdSl] = derivativeRelPerm(this,Sl)
            dklrdSl = this.liqRelativePermeability.derivative(Sl, this.porousMedia);
            if isnan(dklrdSl) || isinf(dklrdSl)
                error("Liquid relative permeability derivative is NaN or Inf.");
            end
            dkgrdSl = this.gasRelativePermeability.derivative(Sl, this.porousMedia);
            if isnan(dkgrdSl) || isinf(dkgrdSl)
                error("Gas relative permeability derivative is NaN or Inf.");
            end
        end

        %------------------------------------------------------------------
        % Compute the derivative of the gas density wrt to the gas pressure
        function drhogdpg = derivativeGasDensityWrtGasPressure(this,rhog,pg)
            % Perturbation value
            pert = sqrt(eps);
            % Compute the gas density given a perturbation at the gas pressure
            rhogPert = this.gasFluid.getDensity(pg + pert);
            % Compute the derivative
            drhogdpg = (rhogPert - rhog)/(pert);
        end

    end
end
