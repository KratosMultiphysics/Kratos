# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.RomApplication as RomApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

# Other imports

def CreateSolver(model, custom_settings):
    return ModalDerivativeSolver(model, custom_settings)

class ModalDerivativeSolver(MechanicalSolver):
    """The structural mechanics modal derivative solver.

    This class creates the modal derivative solver for computing modal derivatives.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        """Initializes the solver."""
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[ModalDerivativeSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "derivative_type"               : "static",
            "derivative_parameter"          : "modal_coordinate",
            "sub_model_parts_list"          : [],
            "finite_difference_type"        : "forward",
            "finite_difference_step_size"   : 1e-3,
            "compute_basis_derivatives"     : true,
            "mass_orthonormalize"           : false,
            "rom_parameters_filename"       : "RomParameters.json",
            "rom_settings"                  :
            {
                "nodal_unknowns" : []
            }
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def _create_solution_scheme(self):
        # Possible derivative parameters
        derivative_parameters = ["modal_coordinate", "density", "poisson_ratio", "young_modulus"]
        finite_difference_types = ["forward", "central"]
        # Check scheme related input
        derivative_parameter = self.settings["derivative_parameter"].GetString()
        if derivative_parameter not in derivative_parameters:
            err_msg  = "\"derivative_parameter\" can only be one of the following:"
            for parameter in derivative_parameters:
                err_msg += " " + parameter
            raise Exception(err_msg)

        finite_difference_type = self.settings["finite_difference_type"].GetString()
        if finite_difference_type not in finite_difference_types:
            err_msg  = "\"finite_difference_type\" can only be one of the following:"
            for fd_type in finite_difference_types:
                err_msg += " " + fd_type
            raise Exception(err_msg)

        # Create scheme
        if derivative_parameter == "modal_coordinate":
            return RomApplication.ModalDerivativeModalCoordinateScheme(self.settings)
        else:
            return RomApplication.ModalDerivativeMaterialParameterScheme(self.settings)

    def _create_builder_and_solver(self):

        if not self.settings["builder_and_solver_settings"]["use_block_builder"].GetBool():
            err_msg  = "Modal derivative analysis is only available with block builder and solver!"
            raise Exception(err_msg)

        return super()._create_builder_and_solver()

    def _create_mechanical_solution_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        modal_derivative_scheme = self.get_solution_scheme()
        builder_and_solver = self.get_builder_and_solver()

        # Check strategy related input
        sub_model_parts_list = self.settings["sub_model_parts_list"].GetStringArray()
        num_sub_model_parts = len(sub_model_parts_list)
        derivative_parameter = self.settings["derivative_parameter"].GetString()
        if derivative_parameter != "modal_coordinate" and num_sub_model_parts == 0:
            err_msg = "\"sub_model_parts_list\" is empty!"
            err_msg += " Please provide at least one SubModelPart"
            raise Exception(err_msg)
        
        for sub_model_part in sub_model_parts_list:
            if not computing_model_part.HasSubModelPart(sub_model_part):
                err_msg  = "\""+sub_model_part+"\" is not a SubModelPart of \""+computing_model_part.Name+"\""
                raise Exception(err_msg)

        derivative_types = ["static", "dynamic"]
        derivative_type = self.settings["derivative_type"].GetString()
        if derivative_type not in derivative_types:
            err_msg  = "\"derivative_type\" can only be one of the following:"
            for der_type in derivative_types:
                err_msg += " " + der_type
            raise Exception(err_msg)

        # Create strategy
        return RomApplication.ModalDerivativeStrategy(computing_model_part,
                                                          modal_derivative_scheme,
                                                          builder_and_solver,
                                                          self.settings)
    
