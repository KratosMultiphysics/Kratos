# importing the Kratos Library
import KratosMultiphysics as Kratos
from . import recoverer
from . import L2_projection_recoverer
import KratosMultiphysics.SwimmingDEMApplication.parameters_tools as PT
import sys

class Pouliot2012GradientRecoverer(L2_projection_recoverer.L2ProjectionGradientRecoverer):
    def __init__(self, project_parameters, model_part):
        L2_projection_recoverer.L2ProjectionGradientRecoverer.__init__(self, project_parameters, model_part)
        self.element_type = "ComputeGradientPouliot20123D"
        self.condition_type = "ComputeLaplacianSimplexCondition3D"
        self.FillUpModelPart(self.element_type, self.condition_type)
        self.DOFs = (Kratos.VELOCITY_COMPONENT_GRADIENT_X, Kratos.VELOCITY_COMPONENT_GRADIENT_Y, Kratos.VELOCITY_COMPONENT_GRADIENT_Z)
        self.AddDofs(self.DOFs)
        self.calculate_vorticity = (project_parameters["vorticity_calculation_type"].GetInt() > 0
                                    or PT.RecursiveFindParametersWithCondition(project_parameters["properties"],
                                                                               'vorticity_induced_lift_parameters'))

class Pouliot2012MaterialAccelerationRecoverer(Pouliot2012GradientRecoverer, L2_projection_recoverer.L2ProjectionMaterialAccelerationRecoverer):
    def __init__(self, model_part, project_parameters, do_pre_recovery = False):
        L2_projection_recoverer.L2ProjectionMaterialAccelerationRecoverer.__init__(self, project_parameters, model_part)
        Pouliot2012GradientRecoverer.__init__(self, project_parameters, model_part)
        self.do_pre_recovery = do_pre_recovery

        scheme = Kratos.ResidualBasedIncrementalUpdateStaticScheme()
        amgcl_smoother = Kratos.AMGCLSmoother.SPAI0
        amgcl_krylov_type = Kratos.AMGCLIterativeSolverType.BICGSTAB_WITH_GMRES_FALLBACK
        tolerance = 1e-12
        max_iterations = 200
        verbosity = 2 # 0->shows no information, 1->some information, 2->all the information
        gmres_size = 400

        if self.use_lumped_mass_matrix:
            linear_solver = Kratos.CGSolver()
        else:
            linear_solver = Kratos.AMGCLSolver(amgcl_smoother, amgcl_krylov_type, tolerance, max_iterations, verbosity,gmres_size)

        self.recovery_strategy = Kratos.ResidualBasedDerivativeRecoveryStrategy(self.recovery_model_part, scheme, linear_solver, False, True, False, False)
        self.recovery_strategy.SetEchoLevel(0)

class Pouliot2012LaplacianRecoverer(L2_projection_recoverer.L2ProjectionDerivativesRecoverer, recoverer.LaplacianRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.LaplacianRecoverer.__init__(self, project_parameters, model_part)
        self.element_type = "ComputeLaplacianSimplex3D"
        self.condition_type = "ComputeLaplacianSimplexCondition3D"
        self.FillUpModelPart(self.element_type, self.condition_type)
        self.DOFs = (Kratos.VELOCITY_LAPLACIAN_X, Kratos.VELOCITY_LAPLACIAN_Y, Kratos.VELOCITY_LAPLACIAN_Z)
        self.AddDofs(self.DOFs)
    def RecoverVectorLaplacian(self, vector_variable, laplacian_variable):
        self.SetToZero(Kratos.VELOCITY_LAPLACIAN)
        self.recovery_strategy.Solve()

    def Solve(self):
        Kratos.Logger.PrintInfo("SwimmingDEM", "\nSolving for the fluid acceleration...")
        sys.stdout.flush()
        self.SetToZero(Kratos.VELOCITY_COMPONENT_GRADIENT)
        if self.do_pre_recovery:
            self.recovery_strategy.Solve()
        self.recovery_strategy.Solve()
