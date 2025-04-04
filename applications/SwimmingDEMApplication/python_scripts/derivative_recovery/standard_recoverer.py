import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as Fluid

from . import recoverer


class StandardGradientRecoverer(recoverer.GradientRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.GradientRecoverer.__init__(self, project_parameters, model_part)

    def RecoverGradientOfScalar(self, scalar_variable, gradient_variable):
        self.cplusplus_recovery_tool.CalculateGradient(self.model_part, scalar_variable, gradient_variable)

    def RecoverGradientOfVector(self, vector_variable, gradient_variable_x, gradient_variable_y, gradient_variable_z):
        self.cplusplus_recovery_tool.CalculateGradient(self.model_part, vector_variable, gradient_variable_x, gradient_variable_y, gradient_variable_z)

    def RecoverGradientOfVelocity(self):
        self.RecoverGradientOfVector(Kratos.VELOCITY, Kratos.VELOCITY_X_GRADIENT, Kratos.VELOCITY_Y_GRADIENT, Kratos.VELOCITY_Z_GRADIENT)

    def RecoverPressureGradient(self):
        self.RecoverGradientOfScalar(Kratos.PRESSURE, Fluid.RECOVERED_PRESSURE_GRADIENT)


class StandardMaterialAccelerationRecoverer(recoverer.MaterialAccelerationRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.MaterialAccelerationRecoverer.__init__(self, project_parameters, model_part)
        self.compute_exact_L2 = project_parameters['compute_exact_L2'].GetBool()

    def RecoverMaterialAcceleration(self):
        if self.compute_exact_L2:
            # self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeExactL2(self.model_part, Kratos.VELOCITY, Kratos.ACCELERATION, Kratos.MATERIAL_ACCELERATION)
            self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeExactL2Parallel(self.model_part, Kratos.VELOCITY, Kratos.ACCELERATION, Kratos.MATERIAL_ACCELERATION)
        else:
            self.cplusplus_recovery_tool.CalculateVectorMaterialDerivative(self.model_part, Kratos.VELOCITY, Kratos.ACCELERATION, Kratos.MATERIAL_ACCELERATION)


class StandardLaplacianRecoverer(recoverer.LaplacianRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.LaplacianRecoverer.__init__(self, project_parameters, model_part)

    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        self.cplusplus_recovery_tool.CalculateVectorLaplacian(self.model_part, vector_variable, laplacian_variable)
