from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from . import recoverer
from . import L2_projection_recoverer

class Fortin2012GradientRecoverer(L2_projection_recoverer.L2ProjectionGradientRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        L2_projection_recoverer.L2ProjectionGradientRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)
        self.element_type = "ComputeGradientFortin20123D"
        self.condition_type = "ComputeLaplacianSimplexCondition3D"
        self.FillUpModelPart(self.element_type, self.condition_type)
        self.DOFs = (VELOCITY_Z_GRADIENT_X, VELOCITY_Z_GRADIENT_Y, VELOCITY_Z_GRADIENT_Z)
        self.AddDofs(self.DOFs)
        self.calculate_vorticity = self.pp.CFD_DEM.lift_force_type

class Fortin2012MaterialAccelerationRecoverer(Fortin2012GradientRecoverer, L2_projection_recoverer.L2ProjectionMaterialAccelerationRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        L2_projection_recoverer.L2ProjectionMaterialAccelerationRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)
        Fortin2012GradientRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)

class Fortin2012LaplacianRecoverer(L2_projection_recoverer.L2ProjectionDerivativesRecoverer, recoverer.LaplacianRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        recoverer.LaplacianRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)
        Fortin2012DerivativesRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)
        self.element_type = "ComputeLaplacianSimplex3D"
        self.condition_type = "ComputeLaplacianSimplexCondition3D"
        self.FillUpModelPart(self.element_type, self.condition_type)
        self.DOFs = (VELOCITY_LAPLACIAN_X, VELOCITY_LAPLACIAN_Y, VELOCITY_LAPLACIAN_Z)
        self.AddDofs(self.DOFs)
    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        self.SetToZero(VELOCITY_LAPLACIAN)
        self.recovery_strategy.Solve()
