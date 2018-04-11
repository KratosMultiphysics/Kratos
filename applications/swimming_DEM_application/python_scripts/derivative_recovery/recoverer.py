# This class can be taken as an abstract template for derivation. It can be used
# for default passive behaviour.

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import weakref

class DerivativesRecoverer:
    def __init__(self, pp, model_part):
        self.pp = pp
        self.model_part = model_part
        self.cplusplus_recovery_tool = DerivativeRecoveryTool3D(model_part, pp.CFD_DEM)

class EmptyGradientRecoverer(DerivativesRecoverer):
    def __init__(self, pp, model_part):
        pass
    def RecoverGradientOfScalar(self, scalar_variable, gradient_variable):
        pass
    def RecoverGradientOfVector(self, vector_variable, gradient_variable_x, gradient_variable_y, gradient_variable_z):
        pass
    def RecoverGradientOfVelocity(self):
        pass
    def RecoverFluidFractionGradient(self):
        pass
    def RecoverPressureGradient(self):
        pass

class EmptyMaterialAccelerationRecoverer(DerivativesRecoverer):
    def __init__(self, pp, model_part):
        pass
    def RecoverMaterialAcceleration(self):
        pass
    def RecoverMaterialAccelerationFromGradient(self):
        pass

class EmptyVorticityRecoverer(DerivativesRecoverer):
    def __init__(self, pp, model_part):
        pass
    def RecoverVorticityFromGradient(self):
        pass
    def CalculateVorticityContributionOfTheGradientOfAComponent(self):
        pass

class EmptyLaplacianRecoverer(DerivativesRecoverer):
    def __init__(self, pp, model_part):
        pass
    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        pass
    def RecoverVelocityLaplacian(self):
        pass

class GradientRecoverer(EmptyGradientRecoverer):
    def __init__(self, pp, model_part):
        DerivativesRecoverer.__init__(self, pp, model_part)
    def RecoverGradientOfVelocity(self):
        self.RecoverGradientOfVector(VELOCITY, VELOCITY_X_GRADIENT, VELOCITY_Y_GRADIENT, VELOCITY_Z_GRADIENT)
    def RecoverPressureGradient(self):
        self.RecoverGradientOfScalar(PRESSURE, PRESSURE_GRADIENT)
    def RecoverFluidFractionGradient(self):
        self.RecoverGradientOfScalar(FLUID_FRACTION, FLUID_FRACTION_GRADIENT)

class MaterialAccelerationRecoverer(GradientRecoverer, EmptyMaterialAccelerationRecoverer):
    def __init__(self, pp, model_part):
        DerivativesRecoverer.__init__(self, pp, model_part)
    def RecoverMaterialAcceleration(self):
        self.RecoverMaterialAccelerationFromGradient()
    def RecoverMaterialAccelerationFromGradient(self):
        self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeFromGradient(self.model_part, VELOCITY_X_GRADIENT, VELOCITY_Y_GRADIENT, VELOCITY_Z_GRADIENT, ACCELERATION, MATERIAL_ACCELERATION)

class VorticityRecoverer(GradientRecoverer, EmptyVorticityRecoverer):
    def __init__(self, pp, model_part):
        DerivativesRecoverer.__init__(self, pp, model_part)
    def RecoverVorticityFromGradient(self):
        self.cplusplus_recovery_tool.CalculateVorticityFromGradient(self.model_part, VELOCITY_X_GRADIENT, VELOCITY_Y_GRADIENT, VELOCITY_Z_GRADIENT, VORTICITY)
    def CalculateVorticityContributionOfTheGradientOfAComponent(self):
        self.cplusplus_recovery_tool.CalculateVorticityContributionOfTheGradientOfAComponent(self.model_part, VELOCITY_COMPONENT_GRADIENT, VORTICITY)

class LaplacianRecoverer(GradientRecoverer, EmptyLaplacianRecoverer):
    def __init__(self, pp, model_part):
        DerivativesRecoverer.__init__(self, pp, model_part)
    def RecoverVelocityLaplacian(self):
        self.RecoverVectorLaplacian(VELOCITY, VELOCITY_LAPLACIAN)
