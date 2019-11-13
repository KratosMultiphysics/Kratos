from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

from KratosMultiphysics import eigen_solver_factory
from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
import KratosMultiphysics.kratos_utilities as kratos_utils

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
        self.solverFlag = "eigen"
        KratosMultiphysics.Logger.PrintInfo("::[PrebucklingSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "buckling_settings"     : {
                "initial_step"          : 1.0,
                "small_step"            : 0.0005,
                "big_step"              : 0.5,
                "convergence_ratio"     : 0.005
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
    def SetSolverFlag( self, flag ):
        self.solverFlag = flag

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
    
    #To allow two different solvers
    def get_builder_and_solver(self):
        self._builder_and_solver = self._create_builder_and_solver()
        return self._builder_and_solver

    def get_linear_solver(self): 
        self._linear_solver = self._create_linear_solver()
        return self._linear_solver

    def _create_linear_solver(self):
        """Create the eigensolver or linear solver depending on SolverFlag"""
        if (self.solverFlag == "eigen"):
            return eigen_solver_factory.ConstructSolver(self.settings["eigensolver_settings"])
        elif( self.solverFlag == "static"):
            linear_solver_configuration = self.settings["linear_solver_settings"]
            # using a default linear solver (selecting the fastest one available)
            if kratos_utils.CheckIfApplicationsAvailable("EigenSolversApplication"):
                from KratosMultiphysics import EigenSolversApplication
            elif kratos_utils.CheckIfApplicationsAvailable("ExternalSolversApplication"):
                from KratosMultiphysics import ExternalSolversApplication

            linear_solvers_by_speed = [
                "pardiso_lu", # EigenSolversApplication (if compiled with IntAnalysisel-support)
                "sparse_lu",  # EigenSolversApplication
                "pastix",     # ExternalSolversApplication (if Pastix is included in compilation)
                "super_lu",   # ExternalSolversApplication
                "skyline_lu_factorization" # in Core, always available, but slow
            ]

            for solver_name in linear_solvers_by_speed:
                if KratosMultiphysics.LinearSolverFactory().Has(solver_name):
                    linear_solver_configuration.AddEmptyValue("solver_type").SetString(solver_name)
                    KratosMultiphysics.Logger.PrintInfo('::[MechanicalSolver]:: ',\
                        'Using "' + solver_name + '" as default linear solver')
                    return KratosMultiphysics.LinearSolverFactory().Create(linear_solver_configuration)
            
    def _create_mechanical_solution_strategy(self):
        solution_scheme = self.get_solution_scheme() 
        self.SetSolverFlag( "eigen")
        eigen_solver = self.get_builder_and_solver() # The eigensolver is created here.
        self.SetSolverFlag( "static")
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
                                                                  buckling_settings["initial_step"].GetDouble(),
                                                                  buckling_settings["small_step"].GetDouble(),
                                                                  buckling_settings["big_step"].GetDouble(),
                                                                  buckling_settings["convergence_ratio"].GetDouble() )
