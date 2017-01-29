from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from . import recoverer

class StandardGradientRecoverer(recoverer.GradientRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        recoverer.GradientRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)

    def RecoverGradientOfScalar(self, scalar_variable, gradient_variable):
        self.CalculateGradient(scalar_variable, gradient_variable)

    def RecoverGradientOfVector(self, scalar_variable, gradient_variable):
        self.CalculateGradient(vector_variable, gradient_variable_x, gradient_variable_y, gradient_variable_z)

class StandardMaterialAccelerationRecoverer(recoverer.MaterialAccelerationRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        recoverer.MaterialAccelerationRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)

    def RecoverMaterialAcceleration(self):
        self.cplusplus_recovery_tool.CalculateVectorMaterialDerivative(self.model_part, VELOCITY, ACCELERATION, MATERIAL_ACCELERATION)

class StandardLaplacianRecoverer(recoverer.LaplacianRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        recoverer.LaplacianRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)

    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        self.cplusplus_recovery_tool.CalculateVectorLaplacian(self.model_part, vector_variable, laplacian_variable)
