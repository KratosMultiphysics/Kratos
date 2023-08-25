# This class can be taken as an abstract template for derivation. It can be used
# for default passive behaviour.

# importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.SwimmingDEMApplication as SDEM

class DerivativesRecoverer:
    def __init__(self, project_parameters, model_part):
        self.model_part = model_part
        domain_size = model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        if domain_size == 2:
            self.cplusplus_recovery_tool = SDEM.DerivativeRecoveryTool2D(model_part, project_parameters)
        else:
            self.cplusplus_recovery_tool = SDEM.DerivativeRecoveryTool3D(model_part, project_parameters)

class EmptyGradientRecoverer(DerivativesRecoverer):
    def __init__(self, project_parameters, model_part):
        DerivativesRecoverer.__init__(self, project_parameters, model_part)
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
    def __init__(self, project_parameters, model_part):
        pass
    def RecoverMaterialAcceleration(self):
        pass
    def RecoverMaterialAccelerationFromGradient(self):
        pass

class EmptyVorticityRecoverer(DerivativesRecoverer):
    def __init__(self, project_parameters, model_part):
        pass
    def RecoverVorticityFromGradient(self):
        pass
    def CalculateVorticityContributionOfTheGradientOfAComponent(self):
        pass

class EmptyLaplacianRecoverer(DerivativesRecoverer):
    def __init__(self, project_parameters, model_part):
        pass
    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        pass
    def RecoverVelocityLaplacian(self):
        pass

class GradientRecoverer(EmptyGradientRecoverer):
    def __init__(self, project_parameters, model_part):
        DerivativesRecoverer.__init__(self, project_parameters, model_part)
    def RecoverGradientOfVelocity(self):
        self.RecoverGradientOfVector(Kratos.VELOCITY, Kratos.VELOCITY_X_GRADIENT, Kratos.VELOCITY_Y_GRADIENT, Kratos.VELOCITY_Z_GRADIENT)
    def RecoverPressureGradient(self):
        self.RecoverGradientOfScalar(Kratos.PRESSURE, Kratos.PRESSURE_GRADIENT)
    def RecoverFluidFractionGradient(self):
        self.RecoverGradientOfScalar(Kratos.FLUID_FRACTION, Kratos.FLUID_FRACTION_GRADIENT)

class MaterialAccelerationRecoverer(GradientRecoverer, EmptyMaterialAccelerationRecoverer):
    def __init__(self, project_parameters, model_part):
        GradientRecoverer.__init__(self, project_parameters, model_part)
    def RecoverMaterialAcceleration(self):
        self.RecoverMaterialAccelerationFromGradient()
    def RecoverMaterialAccelerationFromGradient(self):
        self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeFromGradient(self.model_part,
                                                                                   Kratos.VELOCITY_X_GRADIENT,
                                                                                   Kratos.VELOCITY_Y_GRADIENT,
                                                                                   Kratos.VELOCITY_Z_GRADIENT,
                                                                                   Kratos.ACCELERATION,
                                                                                   Kratos.MATERIAL_ACCELERATION)

class VorticityRecoverer(GradientRecoverer, EmptyVorticityRecoverer):
    def __init__(self, project_parameters, model_part):
        GradientRecoverer.__init__(self, project_parameters, model_part)
    def RecoverVorticityFromGradient(self):
        self.cplusplus_recovery_tool.CalculateVorticityFromGradient(self.model_part,
                                                                    Kratos.VELOCITY_X_GRADIENT,
                                                                    Kratos.VELOCITY_Y_GRADIENT,
                                                                    Kratos.VELOCITY_Z_GRADIENT,
                                                                    Kratos.VORTICITY)

    def CalculateVorticityContributionOfTheGradientOfAComponent(self):
        self.cplusplus_recovery_tool.CalculateVorticityContributionOfTheGradientOfAComponent(self.model_part,
                                                                                             Kratos.VELOCITY_COMPONENT_GRADIENT,
                                                                                             Kratos.VORTICITY)

class LaplacianRecoverer(GradientRecoverer, EmptyLaplacianRecoverer):
    def __init__(self, project_parameters, model_part):
        GradientRecoverer.__init__(self, project_parameters, model_part)
    def RecoverVelocityLaplacian(self):
        self.RecoverVectorLaplacian(Kratos.VELOCITY, Kratos.VELOCITY_LAPLACIAN)
