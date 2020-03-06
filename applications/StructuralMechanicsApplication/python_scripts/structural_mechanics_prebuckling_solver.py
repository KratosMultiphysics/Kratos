from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

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
        super(PrebucklingSolver, self).__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[PrebucklingSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "buckling_settings"     : {
                "initial_load_increment"    : 1.0,
                "small_load_increment"      : 0.0005,
                "path_following_step"       : 0.5,
                "convergence_ratio"         : 0.05
            },
            "eigensolver_settings" : {
                "solver_type"           : "eigen_eigensystem",
                "max_iteration"         : 1000,
                "tolerance"             : 1e-6,
                "number_of_eigenvalues" : 5,
                "echo_level"            : 1
            }
        }""")
        this_defaults.AddMissingParameters(super(PrebucklingSolver, cls).GetDefaultSettings())
        return this_defaults

    #### Private functions ####
    def AdvanceInTime(self, current_time):
        new_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.main_model_part.CloneTimeStep(new_time)

        return new_time

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    # Builder and Solver Eigen
    def get_builder_and_solver_eigen(self):
        if not hasattr(self, '_builder_and_solver_eigen'):
            self._builder_and_solver_eigen = self._create_builder_and_solver_eigen()
        return self._builder_and_solver_eigen

    def _create_builder_and_solver_eigen(self):
        linear_solver = self.get_linear_solver_eigen()
        builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        return builder_and_solver

    def get_linear_solver_eigen(self):
        if not hasattr(self, '_linear_solver_eigen'):
            self._linear_solver_eigen = self._create_linear_solver_eigen()
        return self._linear_solver_eigen

    def _create_linear_solver_eigen(self):
        """Create the eigensolver"""
        return eigen_solver_factory.ConstructSolver(self.settings["eigensolver_settings"])

    # Builder and Solver Static
    def _create_builder_and_solver(self):
        """This method is overridden to make sure it always uses ResidualBasedEliminationBuilderAndSolver"""
        if self.settings["block_builder"].GetBool():
            warn_msg = '"Elimination Builder is required. \n'
            warn_msg += '"block_builder" specification will be ignored'
            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsPrebucklingAnalysis; Warning", warn_msg)
        linear_solver = self.get_linear_solver()
        builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        return builder_and_solver

    # Convergence Criterion
    def _create_convergence_criterion(self):
        """This method is overridden to make sure it always uses "displacement_criterion" """
        convergence_criterion_setting = self._get_convergence_criterion_settings()
        if convergence_criterion_setting["convergence_criterion"].GetString() != "displacement_criterion":
            warn_msg = 'Convergence criterion "displacement_criterion" is required. \n'
            warn_msg += '"' + convergence_criterion_setting["convergence_criterion"].GetString() + '" specification will be ignored'
            KratosMultiphysics.Logger.PrintWarning("StructuralMechanicsPrebucklingAnalysis; Warning", warn_msg)
            convergence_criterion_setting["convergence_criterion"].SetString("displacement_criterion")

        convergence_criterion = convergence_criteria_factory.convergence_criterion(convergence_criterion_setting)
        return convergence_criterion.mechanical_convergence_criterion

    def _create_mechanical_solution_strategy(self):
        solution_scheme = self.get_solution_scheme()
        eigen_solver = self.get_builder_and_solver_eigen() # The eigensolver is created here.
        builder_and_solver = self.get_builder_and_solver() # The linear solver is created here.
        convergence_criteria = self.get_convergence_criterion()
        computing_model_part = self.GetComputingModelPart()
        buckling_settings = self.settings["buckling_settings"]

        return StructuralMechanicsApplication.PrebucklingStrategy(computing_model_part,
                                                                  solution_scheme,
                                                                  eigen_solver,
                                                                  builder_and_solver,
                                                                  convergence_criteria,
                                                                  self.settings["max_iteration"].GetInt(),
                                                                  buckling_settings["initial_load_increment"].GetDouble(),
                                                                  buckling_settings["small_load_increment"].GetDouble(),
                                                                  buckling_settings["path_following_step"].GetDouble(),
                                                                  buckling_settings["convergence_ratio"].GetDouble() )
