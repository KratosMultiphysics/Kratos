%% CapillaryPressureUMAT class
% This class defines a capillary pressure constitutive law for porous 
% media using a user-defined material model (UMAT). It inherits from the 
% _CapillaryPressure_ base class.
%
%% Methods
% * *saturationDegree*: Computes the liquid phase saturation degree Sl for 
%                       a given capillary pressure pc and porous media 
%                       properties. The result is clamped between the 
%                       residual saturation limits (Slr and 1.0 - Sgr).
% * *derivativeSaturationDegree*: Computes the derivative of the 
%                                 saturation degree dSldpc with respect to 
%                                 capillary pressure pc using the selected 
%                                 derivative method.
% * *GetCurveDerivative*: Computes the derivative of the saturation degree 
%                         curve at a given saturation degree. The method 
%                         supports two approaches: Piecewise constant, 
%                         which corresponds to derivativeMethod = 0 and 
%                         Smoothed, which is derivativeMethod = 1, and it 
%                         uses weighted averaging of slopes for smoother 
%                         transitions.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef CapillaryPressureUMAT < CapillaryPressure 
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        pc_curve = [];          % Must be sorted!!
        Sl_curve = [];          % Must be sorted!!
        derivativeMethod = 1;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = CapillaryPressureUMAT(pc_curve,Sl_curve)
            this = this@CapillaryPressure('umat');
            this.pc_curve = pc_curve;
            this.Sl_curve = Sl_curve;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the liquid phase relative permeability
        function Sl = saturationDegree(this, pc, porousMedia)
            Sl = interp1(this.pc_curve,this.Sl_curve,pc,'linear', 'extrap');
            Sl = max(min(Sl, 1.0-porousMedia.Sgr-eps), porousMedia.Slr+eps);
        end
        
        %------------------------------------------------------------------
        % Compute the gas phase relative permeability
        function dSldpc = derivativeSaturationDegree(this, pc, porousMedia)
            % Compute the saturation degree at the perturbed values
            Sl = this.saturationDegree(pc, porousMedia);
            dSldpc = this.GetCurveDerivative(Sl);
        end

        %------------------------------------------------------------------
        % Computes the derivative of a given curve
        % Copied from OGS-5
        function derivative = GetCurveDerivative(this, xValue)
            
            % Extract the x and y values
            xPoints = this.Sl_curve;
            yPoints = this.pc_curve;
            numPoints = size(xPoints, 1);
        
            % Handle out-of-bound xValue
            if xValue < xPoints(1)
                xValue = xPoints(1);
                index = 1;
            elseif xValue > xPoints(end)
                xValue = xPoints(end);
                index = numPoints;
            else
                % Locate the interval containing xValue
                index = find(xPoints >= xValue, 1);
            end
        
            switch this.derivativeMethod
                case 0  % Piecewise constant
                    if index > 1
                        deltaX = xPoints(index) - xPoints(index - 1);
                        if abs(deltaX) > eps
                            derivative = (yPoints(index) - yPoints(index - 1)) / deltaX;
                        else
                            derivative = sign(yPoints(index) - yPoints(index - 1)) / eps;
                        end
                    else
                        derivative = 0; % Undefined derivative at the start of the curve
                    end
        
                case 1  % Smoothed
                    if index > 2 && index < numPoints
                        % Compute slopes on either side
                        slope1 = (yPoints(index) - yPoints(index - 2)) / ...
                                 (xPoints(index) - xPoints(index - 2));
                        slope2 = (yPoints(index + 1) - yPoints(index - 1)) / ...
                                 (xPoints(index + 1) - xPoints(index - 1));
                        
                        % Linear interpolation weight
                        w = (xValue - xPoints(index - 1)) / (xPoints(index) - xPoints(index - 1));
        
                        % Weighted average of slopes
                        derivative = (1 - w) * slope1 + w * slope2;
                    else
                        % Fallback to piecewise constant for boundaries
                        if index > 1
                            deltaX = xPoints(index) - xPoints(index - 1);
                            if abs(deltaX) > eps
                                derivative = (yPoints(index) - yPoints(index - 1)) / deltaX;
                            else
                                derivative = sign(yPoints(index) - yPoints(index - 1)) / eps;
                            end
                        else
                            derivative = 0; % Undefined derivative at the start of the curve
                        end
                    end
        
                otherwise
                    error('Unknown method: %d', this.derivativeMethod);
            end
            derivative = 1.0 / derivative;
        end 
    end
end