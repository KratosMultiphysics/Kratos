from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from . import recoverer
from . import L2_projection_recoverer

class Pouliot2012GradientRecoverer(L2_projection_recoverer.L2ProjectionGradientRecoverer):
    def __init__(self, pp, model_part):
        L2_projection_recoverer.L2ProjectionGradientRecoverer.__init__(self, pp, model_part)
        self.element_type = "ComputeGradientPouliot20123D"
        self.condition_type = "ComputeLaplacianSimplexCondition3D"
        self.FillUpModelPart(self.element_type, self.condition_type)
        self.DOFs = (VELOCITY_COMPONENT_GRADIENT_X, VELOCITY_COMPONENT_GRADIENT_Y, VELOCITY_COMPONENT_GRADIENT_Z)
        self.AddDofs(self.DOFs)
        self.calculate_vorticity = self.pp.CFD_DEM["lift_force_type"].GetInt()

class Pouliot2012MaterialAccelerationRecoverer(Pouliot2012GradientRecoverer, L2_projection_recoverer.L2ProjectionMaterialAccelerationRecoverer):
    def __init__(self, pp, model_part, do_pre_recovery = False):
        L2_projection_recoverer.L2ProjectionMaterialAccelerationRecoverer.__init__(self, pp, model_part)
        Pouliot2012GradientRecoverer.__init__(self, pp, model_part)
        self.do_pre_recovery = do_pre_recovery

        scheme = ResidualBasedIncrementalUpdateStaticScheme()
        amgcl_smoother = AMGCLSmoother.SPAI0
        amgcl_krylov_type = AMGCLIterativeSolverType.BICGSTAB_WITH_GMRES_FALLBACK
        tolerance = 1e-12
        max_iterations = 200
        verbosity = 2 # 0->shows no information, 1->some information, 2->all the information
        gmres_size = 400

        if self.use_lumped_mass_matrix:
            linear_solver = CGSolver()
        else:
            linear_solver = AMGCLSolver(amgcl_smoother, amgcl_krylov_type, tolerance, max_iterations, verbosity,gmres_size)

        self.recovery_strategy = ResidualBasedDerivativeRecoveryStrategy(self.recovery_model_part, scheme, linear_solver, False, True, False, False)
        self.recovery_strategy.SetEchoLevel(0)

class Pouliot2012LaplacianRecoverer(L2_projection_recoverer.L2ProjectionDerivativesRecoverer, recoverer.LaplacianRecoverer):
    def __init__(self, pp, model_part):
        recoverer.LaplacianRecoverer.__init__(self, pp, model_part)
        Pouliot2012DerivativesRecoverer.__init__(self, pp, model_part)
        self.element_type = "ComputeLaplacianSimplex3D"
        self.condition_type = "ComputeLaplacianSimplexCondition3D"
        self.FillUpModelPart(self.element_type, self.condition_type)
        self.DOFs = (VELOCITY_LAPLACIAN_X, VELOCITY_LAPLACIAN_Y, VELOCITY_LAPLACIAN_Z)
        self.AddDofs(self.DOFs)
    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        self.SetToZero(VELOCITY_LAPLACIAN)
        self.recovery_strategy.Solve()

    def Solve(self):
        print("\nSolving for the fluid acceleration...")
        sys.stdout.flush()
        self.SetToZero(VELOCITY_COMPONENT_GRADIENT)
        if self.do_pre_recovery:
            self.aux_recovery_strategy.Solve()
        self.recovery_strategy.Solve()
