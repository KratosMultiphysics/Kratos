# This class can be taken as an abstract template for derivation. It can be used
# for default passive behaviour.

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *

class DerivativesRecoverer:
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        self.pp = pp
        self.model_part = model_part
        self.cplusplus_recovery_tool = cplusplus_recovery_tool

class GradientRecoverer(DerivativesRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        DerivativesRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)
    def RecoverGradientOfScalar(self, scalar_variable, gradient_variable):
        pass
    def RecoverGradientOfVector(self, vector_variable, gradient_variable_x, gradient_variable_y, gradient_variable_z):
        pass

class MaterialAccelerationRecoverer(DerivativesRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        DerivativesRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)
    def RecoverMaterialAcceleration(self):
        pass
    def RecoverMaterialAccelerationFromGradient(self):
        self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeFromGradient(self.model_part, VELOCITY_X_GRADIENT, VELOCITY_Y_GRADIENT, VELOCITY_Z_GRADIENT, ACCELERATION, MATERIAL_ACCELERATION)

class VorticityRecoverer(DerivativesRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        DerivativesRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)
    def RecoverVorticityFromGradient(self):
        self.cplusplus_recovery_tool.CalculateVorticityFromGradient(self.model_part, VELOCITY_X_GRADIENT, VELOCITY_Y_GRADIENT, VELOCITY_Z_GRADIENT, VORTICITY)
    def CalculateVorticityContributionOfTheGradientOfAComponent(self):
        self.cplusplus_recovery_tool.CalculateVorticityContributionOfTheGradientOfAComponent(self.model_part, VELOCITY_Z_GRADIENT, VORTICITY)

class LaplacianRecoverer(DerivativesRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        DerivativesRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)
    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        pass

class EmptyGradientRecoverer(GradientRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        pass
    def RecoverGradientOfScalar(self, scalar_variable, gradient_variable):
        pass
    def RecoverGradientOfVector(self, vector_variable, gradient_variable_x, gradient_variable_y, gradient_variable_z):
        pass

class EmptyMaterialAccelerationRecoverer(MaterialAccelerationRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        pass
    def RecoverMaterialAcceleration(self):
        pass
    def RecoverMaterialAccelerationFromGradient(self):
        pass

class EmptyVorticityRecoverer(VorticityRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        pass
    def RecoverVorticityFromGradient(self):
        pass
    def CalculateVorticityContributionOfTheGradientOfAComponent(self):
        pass

class EmptyLaplacianRecoverer(LaplacianRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        pass
    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        pass
