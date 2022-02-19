# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.RomApplication as RomApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

# Other imports
import json

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
            "sub_model_parts_list"          : [],
            "finite_difference_type"        : "forward",
            "finite_difference_step_size"   : 1e-3,
            "mass_orthonormalize"           : false,
            "rom_parameters_filename"       : "RomParameters.json"
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def _CreateScheme(self):
        # Possible derivative parameters

        finite_difference_types = ["forward", "central"]
        # Check scheme related input
        finite_difference_type = self.settings["finite_difference_type"].GetString()
        if finite_difference_type not in finite_difference_types:
            err_msg  = "\"finite_difference_type\" can only be one of the following:"
            for fd_type in finite_difference_types:
                err_msg += " " + fd_type
            raise Exception(err_msg)

        # Create scheme
        return RomApplication.ModalDerivativeScheme(self.settings)
        
    def _CreateBuilderAndSolver(self):

        if not self.settings["builder_and_solver_settings"]["use_block_builder"].GetBool():
            err_msg  = "Modal derivative analysis is only available with block builder and solver!"
            raise Exception(err_msg)

        return super()._CreateBuilderAndSolver()

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        modal_derivative_scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()

        # Check strategy related input
        derivative_types = ["static", "dynamic"]
        derivative_type = self.settings["derivative_type"].GetString()
        if derivative_type not in derivative_types:
            err_msg  = "\"derivative_type\" can only be one of the following:"
            for der_type in derivative_types:
                err_msg += " " + der_type
            raise Exception(err_msg)

        rom_parameters_filename = self.settings["rom_parameters_filename"].GetString()
        with open(rom_parameters_filename, "r") as rom_parameters_file:
            rom_parameters = json.load(rom_parameters_file)
            nodal_unknowns = rom_parameters["rom_settings"]["nodal_unknowns"]
            self.settings.AddEmptyValue("rom_settings")
            self.settings["rom_settings"].AddEmptyValue("nodal_unknowns")
            self.settings["rom_settings"]["nodal_unknowns"].SetStringArray(nodal_unknowns)

        # Create strategy
        return RomApplication.ModalDerivativeStrategy(computing_model_part,
                                                          modal_derivative_scheme,
                                                          builder_and_solver,
                                                          self.settings)