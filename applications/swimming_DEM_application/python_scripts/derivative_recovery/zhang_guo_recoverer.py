from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from . import recoverer

class ZhangGuoGradientRecoverer(recoverer.GradientRecoverer):
    def __init__(self, pp, model_part):
        recoverer.GradientRecoverer.__init__(self, pp, model_part)
    def RecoverGradientOfScalar(self, scalar_variable, gradient_variable):
        self.cplusplus_recovery_tool.RecoverSuperconvergentGradient(self.model_part, scalar_variable, gradient_variable)
    def RecoverGradientOfVelocity(self):
        self.cplusplus_recovery_tool.RecoverSuperconvergentGradient(self.model_part, VELOCITY_X, VELOCITY_X_GRADIENT)
        self.cplusplus_recovery_tool.RecoverSuperconvergentGradient(self.model_part, VELOCITY_Y, VELOCITY_Y_GRADIENT)
        self.cplusplus_recovery_tool.RecoverSuperconvergentGradient(self.model_part, VELOCITY_Z, VELOCITY_Z_GRADIENT)
    def RecoverPressureGradient(self):
        self.cplusplus_recovery_tool.RecoverSuperconvergentGradient(self.model_part, PRESSURE, RECOVERED_PRESSURE_GRADIENT)

class ZhangGuoMaterialAccelerationRecoverer(recoverer.MaterialAccelerationRecoverer, ZhangGuoGradientRecoverer):
    def __init__(self, pp, model_part):
        recoverer.MaterialAccelerationRecoverer.__init__(self, pp, model_part)
    def RecoverMaterialAcceleration(self):
        self.cplusplus_recovery_tool.CalculateVectorMaterialDerivative(self.model_part, VELOCITY, ACCELERATION, MATERIAL_ACCELERATION)

class ZhangGuoDirectLaplacianRecoverer(recoverer.LaplacianRecoverer):
    def __init__(self, pp, model_part):
        recoverer.LaplacianRecoverer.__init__(self, pp, model_part)
    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        self.cplusplus_recovery_tool.RecoverSuperconvergentLaplacian(self.model_part, vector_variable, laplacian_variable)

class ZhangGuoMaterialAccelerationAndLaplacianRecoverer(recoverer.LaplacianRecoverer, ZhangGuoMaterialAccelerationRecoverer):
    def __init__(self, pp, model_part):
        recoverer.LaplacianRecoverer.__init__(self, pp, model_part)
    def RecoverMaterialAcceleration(self):
        self.RecoverGradientOfVelocity()
        self.RecoverMaterialAccelerationFromGradient()
    def RecoverVelocityLaplacian(self):
        self.cplusplus_recovery_tool.RecoverSuperconvergentVelocityLaplacianFromGradient(self.model_part, VELOCITY, VELOCITY_LAPLACIAN)
