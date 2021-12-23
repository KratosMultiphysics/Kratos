# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

from KratosMultiphysics import eigen_solver_factory
from KratosMultiphysics.StructuralMechanicsApplication import convergence_criteria_factory

def CreateSolver(main_model_part, custom_settings):
    return PrebucklingSolver(main_model_part, custom_settings)

class PrebucklingSolver(MechanicalSolver):
    """The structural mechanics prebuckling solver.

    This class creates the mechanical solvers for prebuckling analysis.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Construct the base solver.
        super().__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[PrebucklingSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "buckling_settings"     : {
                "initial_load_increment"    : 1.0,
                "small_load_increment"      : 0.0005,
                "path_following_step"       : 0.5,
                "convergence_ratio"         : 0.05,
                "make_matrices_symmetric"   : true
            },
            "eigensolver_settings" : {
                "solver_type"           : "eigen_eigensystem",
                "max_iteration"         : 1000,
                "tolerance"             : 1e-6,
                "number_of_eigenvalues" : 5,
                "echo_level"            : 1
            }
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    #### Private functions ####
    def AdvanceInTime(self, current_time):
        new_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.main_model_part.CloneTimeStep(new_time)

        return new_time

    def _CreateScheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    # Builder and Solver Eigen
    def _GetBuilderAndSolverEigen(self):
        if not hasattr(self, '_builder_and_solver_eigen'):
            self._builder_and_solver_eigen = self._CreateBuilderAndSolverEigen()
        return self._builder_and_solver_eigen

    def _CreateBuilderAndSolverEigen(self):
        linear_solver = self._GetLinearSolverEigen()
        builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        return builder_and_solver

    def _GetLinearSolverEigen(self):
        if not hasattr(self, '_linear_solver_eigen'):
            self._linear_solver_eigen = self._CreateLinearSolverEigen()
        return self._linear_solver_eigen

    def _CreateLinearSolverEigen(self):
        """Create the eigensolver"""
        return eigen_solver_factory.ConstructSolver(self.settings["eigensolver_settings"])

    # Builder and Solver Static
    def _CreateBuilderAndSolver(self):
        """This method is overridden to make sure it always uses ResidualBasedEliminationBuilderAndSolver"""
        if self.settings["builder_and_solver_settings"]["use_block_builder"].GetBool():
            warn_msg = '"Elimination Builder is required. \n'
            warn_msg += '"use_block_builder" specification will be ignored'
            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsPrebucklingAnalysis; Warning", warn_msg)
        linear_solver = self._GetLinearSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        return builder_and_solver

    # Convergence Criterion
    def _CreateConvergenceCriterion(self):
        """This method is overridden to make sure it always uses "displacement_criterion" """
        convergence_criterion_setting = self._get_convergence_criterion_settings()
        if convergence_criterion_setting["convergence_criterion"].GetString() != "displacement_criterion":
            warn_msg = 'Convergence criterion "displacement_criterion" is required. \n'
            warn_msg += '"' + convergence_criterion_setting["convergence_criterion"].GetString() + '" specification will be ignored'
            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsPrebucklingAnalysis; Warning", warn_msg)
            convergence_criterion_setting["convergence_criterion"].SetString("displacement_criterion")

        convergence_criterion = convergence_criteria_factory.convergence_criterion(convergence_criterion_setting)
        return convergence_criterion.mechanical_convergence_criterion

    def _CreateSolutionStrategy(self):
        solution_scheme = self._GetScheme()
        eigen_solver = self._GetBuilderAndSolverEigen() # The eigensolver is created here.
        builder_and_solver = self._GetBuilderAndSolver() # The linear solver is created here.
        convergence_criteria = self._GetConvergenceCriterion()
        computing_model_part = self.GetComputingModelPart()
        buckling_settings = self.settings["buckling_settings"]

        return StructuralMechanicsApplication.PrebucklingStrategy(computing_model_part,
                                                                  solution_scheme,
                                                                  eigen_solver,
                                                                  builder_and_solver,
                                                                  convergence_criteria,
                                                                  self.settings["max_iteration"].GetInt(),
                                                                  buckling_settings )
