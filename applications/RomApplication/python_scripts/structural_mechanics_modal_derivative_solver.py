from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.RomApplication as RomApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

# Other imports
from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory

def CreateSolver(model, custom_settings):
    return ModalDerivativeSolver(model, custom_settings)

class ModalDerivativeSolver(MechanicalSolver):
    """The structural mechanics modal derivative solver.

    This class creates the modal derivative solver for computing modal derivatives.

    See structural_mechanics_solver.py for more information.
    """
    """
    _create_solution_scheme -> implement a new scheme
    _create_convergence_criterion -> not necessary
    _create_linear_solver -> use from base class 
    _create_builder_and_solver -> do i need a new builder and solver?
    _create_mechanical_solution_strategy -> implement a new strategy

    The mechanical_solution_strategy, builder_and_solver, etc. should alway be retrieved
    using the getter functions get_mechanical_solution_strategy, get_builder_and_solver,
    etc. from this base class.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the model containing the modelpart used to construct the solver.
    settings -- Kratos parameters containing solver settings.
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super(ModalDerivativeSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[ModalDerivativeSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "derivative_type"               : "static",
            "derivative_parameter"          : "modal_coordinates",
            "sub_model_parts_list"          : [],
            "finite_difference_type"        : "forward",
            "finite_difference_step_size"   : 1e-3,
            "compute_basis_derivatives"     : true,
            "mass_orthonormalize"           : true,
            "rom_parameters_filename"       : "RomParameters.json"
        }""")
        this_defaults.AddMissingParameters(super(ModalDerivativeSolver, cls).GetDefaultSettings())
        return this_defaults

    def _create_linear_solver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        if linear_solver_configuration.Has("solver_type"): # user specified a linear solver
            return linear_solver_factory.ConstructSolver(linear_solver_configuration)
        else:
            KratosMultiphysics.Logger.PrintInfo('::[MechanicalSolver]:: No linear solver was specified, using "super_lu" as default solver')
            from KratosMultiphysics import ExternalSolversApplication
            linear_solver_configuration = KratosMultiphysics.Parameters("""{ "solver_type" : "super_lu"}""")
            return KratosMultiphysics.LinearSolverFactory().Create(linear_solver_configuration)

    def _create_solution_scheme(self):

        finite_difference_type = self.settings["finite_difference_type"].GetString()
        if finite_difference_type != "forward" and finite_difference_type != "central":
            err_msg  = '\"finite_difference_type\" can only be \"forward\" or \"central\"!'
            raise Exception(err_msg)

        derivative_parameter = self.settings["derivative_parameter"].GetString()
        if derivative_parameter == "modal_coordinates":
            derivative_parameter = RomApplication.MODAL_COORDINATE
        elif derivative_parameter == "density":
            derivative_parameter = KratosMultiphysics.DENSITY
        elif derivative_parameter == "poisson_ratio":
            derivative_parameter = KratosMultiphysics.POISSON_RATIO
        elif derivative_parameter == "young_modulus":
            derivative_parameter = KratosMultiphysics.YOUNG_MODULUS
        else:
            err_msg  = 'Given \"derivative_parameter\": ', derivative_parameter, ' is not valid!'
            raise Exception(err_msg)

        return RomApplication.ModalDerivativeScheme(derivative_parameter, self.settings)
        
    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        if self.settings["block_builder"].GetBool():
            bs_params = self.settings["builder_and_solver_settings"]
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver, bs_params)
        else:
            if self.settings["multi_point_constraints_used"].GetBool():
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolverWithConstraints(linear_solver)
            else:
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        return builder_and_solver

    def _create_mechanical_solution_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        builder_and_solver = self.get_builder_and_solver()

        return RomApplication.ModalDerivativeStrategy(computing_model_part,
                                                          mechanical_scheme,
                                                          builder_and_solver,
                                                          self.settings)
    