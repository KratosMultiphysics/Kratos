# importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as Fluid
from . import recoverer

import KratosMultiphysics.SwimmingDEMApplication as SDEM

class ZhangGuoGradientRecoverer(recoverer.GradientRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.GradientRecoverer.__init__(self, project_parameters, model_part)
    def RecoverGradientOfScalar(self, scalar_variable, gradient_variable):
        self.cplusplus_recovery_tool.RecoverSuperconvergentGradient(self.model_part, scalar_variable, gradient_variable, self.recovery_order)
    def RecoverGradientOfVelocity(self):
        # print(f"We calling RecoverGradientOfVelocity")
        # exit(0)
        self.cplusplus_recovery_tool.RecoverSuperconvergentGradientAlt(self.model_part, Kratos.VELOCITY_X, Kratos.VELOCITY_X_GRADIENT, self.recovery_order)
        self.cplusplus_recovery_tool.RecoverSuperconvergentGradientAlt(self.model_part, Kratos.VELOCITY_Y, Kratos.VELOCITY_Y_GRADIENT, self.recovery_order)
        self.cplusplus_recovery_tool.RecoverSuperconvergentGradientAlt(self.model_part, Kratos.VELOCITY_Z, Kratos.VELOCITY_Z_GRADIENT, self.recovery_order)

    def RecoverPressureGradient(self):
        self.cplusplus_recovery_tool.RecoverSuperconvergentGradient(self.model_part, Kratos.PRESSURE, Fluid.RECOVERED_PRESSURE_GRADIENT, self.recovery_order)
    def RecoverFluidFractionGradient(self):
        self.cplusplus_recovery_tool.RecoverSuperconvergentGradient(self.model_part, Kratos.FLUID_FRACTION, Kratos.FLUID_FRACTION_GRADIENT, self.recovery_order)

class ZhangGuoMaterialAccelerationRecoverer(recoverer.MaterialAccelerationRecoverer, ZhangGuoGradientRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.MaterialAccelerationRecoverer.__init__(self, project_parameters, model_part)
    # def RecoverMaterialAcceleration(self):
        # self.cplusplus_recovery_tool.CalculateVectorMaterialDerivative(self.model_part, Kratos.VELOCITY, Kratos.ACCELERATION, Kratos.MATERIAL_ACCELERATION)
    def RecoverMaterialAcceleration(self):
        # self.RecoverGradientOfVelocity()  # This is already called by the strategy!
        self.RecoverMaterialAccelerationFromGradient()

class ZhangGuoDirectLaplacianRecoverer(recoverer.LaplacianRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.LaplacianRecoverer.__init__(self, project_parameters, model_part)
    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        self.cplusplus_recovery_tool.RecoverSuperconvergentLaplacian(self.model_part, vector_variable, laplacian_variable)

class ZhangGuoMaterialAccelerationAndLaplacianRecoverer(recoverer.LaplacianRecoverer, ZhangGuoMaterialAccelerationRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.LaplacianRecoverer.__init__(self, project_parameters, model_part)
    def RecoverMaterialAcceleration(self):
        self.RecoverGradientOfVelocity()
        self.RecoverMaterialAccelerationFromGradient()
    def RecoverVelocityLaplacian(self):
        # self.cplusplus_recovery_tool.RecoverSuperconvergentVelocityLaplacianFromGradient(self.model_part, Kratos.VELOCITY, Kratos.VELOCITY_LAPLACIAN)
        # TODO: Uncomment this!!
        pass
