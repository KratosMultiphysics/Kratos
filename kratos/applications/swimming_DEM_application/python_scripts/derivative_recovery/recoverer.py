# This class can be taken as an abstract template for derivation. It can be used
# for default passive behaviour.

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *

class GradientRecoverer:
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        self.pp = pp
        self.model_part = model_part
        self.cplusplus_recovery_tool = cplusplus_recovery_tool

    def RecoverGradientOfScalar(self, scalar_variable, gradient_variable):
        pass

    def RecoverGradientOfVector(self, scalar_variable, gradient_variable):
        pass

class MaterialAccelerationRecoverer:
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        self.pp = pp
        self.model_part = model_part
        self.cplusplus_recovery_tool = cplusplus_recovery_tool
    def RecoverMaterialAcceleration(self):
        self.cplusplus_recovery_tool.CalculateVectorMaterialDerivative(self.model_part, VELOCITY, ACCELERATION, MATERIAL_ACCELERATION)

class LaplacianRecoverer:
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        self.pp = pp
        self.model_part = model_part
        self.cplusplus_recovery_tool = cplusplus_recovery_tool
    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        pass
