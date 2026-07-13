# Import Kratos 
import KratosMultiphysics

# Import the DeprecationManager to handle deprecated variables and provide warnings to the user
from KratosMultiphysics.deprecation_management import DeprecationManager

# Import the os module to handle file paths and directories
import os

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    """
    This function is designed to create an instance of the UnvOutputProcess class, which is responsible for handling the output of simulation results in the UNV format. It takes in a settings object and a model object, both of which are essential for configuring the output process.

    Args:
        settings (KratosMultiphysics.Parameters): A Parameters object that encapsulates a JSON string containing the configuration settings for the output process.
        model (KratosMultiphysics.Model): A Model object that contains the model parts and data necessary for the simulation.
    """
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object, encapsulating a json string")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return UnvOutputProcess(model, settings["Parameters"])

class UnvOutputProcess(KratosMultiphysics.OutputProcess):
    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        """
        This function initializes the UNV output process, setting up the model part, validating the settings, and preparing the UNV output for nodal variables and non-historical variables.

        Args:
            model (KratosMultiphysics.Model): The Kratos model containing the model parts.
            settings (KratosMultiphysics.Parameters): The settings object containing the parameters for the UNV output process.
        """
        KratosMultiphysics.OutputProcess.__init__(self)

        model_part_name = settings["model_part_name"].GetString()
        self.model = model
        self.settings = settings

        # Warning: we may be changing the parameters object here:
        self.TranslateLegacyVariablesAccordingToCurrentStandard(settings)

        # Initialize the UNV output handler with the model and settings. The default settings can be
        # found in "unv_output.cpp"; constructing the handler also validates the settings.
        self.unv_io = KratosMultiphysics.UnvOutput(self.model, self.settings)

        # Ensure the output path exists, creating it if necessary (must happen before writing the mesh)
        if self.settings["save_output_files_in_folder"].GetBool():
            output_path = self.settings["output_path"].GetString()
            if not os.path.exists(output_path):
                os.makedirs(output_path)

        # Prepare the mesh for output
        self.unv_io.InitializeMesh()
        self.unv_io.WriteMesh()

        # Initialize the output controller to manage the output intervals based on the provided settings
        self.controller = KratosMultiphysics.OutputController(model, settings)

    def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings: KratosMultiphysics.Parameters) -> None:
        """
        This function can be extended with new deprecated variables as they are generated

        The purpose of this function is to translate legacy variable names to the current standard, ensuring backward compatibility.
        It checks for the presence of deprecated variable names in the settings and replaces them with the new variable names, while also providing warnings to the user about the changes.

        Args:
            settings (KratosMultiphysics.Parameters): The settings object containing the parameters to be checked
        """
        # Defining a string to help the user understand where the warnings come from (in case any is thrown)
        context_string = type(self).__name__

        # Legacy variable name: "nodal_results"
        old_name = 'nodal_results'
        new_name = 'nodal_solution_step_data_variables'

        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

        # Check if the old variable name exists in the settings and replace it with the new variable name
        old_name = 'output_frequency'
        new_name = 'output_interval'

        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

    def PrintOutput(self) -> None:
        """
        This function prints the output for the nodal variables and nodal non-historical variables
        at the current control value, and schedules the next output.
        """
        # Print the output for the nodal variables and nodal non-historical variables at the current control value
        self.unv_io.PrintOutput()

        # Schedule next output
        self.controller.Update()

    def IsOutputStep(self) -> bool:
        """
        This function checks if the current step is an output step based on the output controller's evaluation.
        
        Returns:
            bool: True if the current step is an output step, False otherwise.
        """
        return self.controller.Evaluate()