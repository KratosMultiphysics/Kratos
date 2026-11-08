# Import necessary modules from Kratos
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    """Creates an EnSightOutputProcess instance.

    This function serves as the entry point for creating the EnSightOutputProcess
    from the project's parameters JSON file. It ensures that the input types
    are correct before instantiating the process.

    Args:
        settings (KratosMultiphysics.Parameters): The configuration parameters for the process.
        model (KratosMultiphysics.Model): The Kratos model containing the model parts.

    Raises:
        Exception: If the 'model' input is not a KratosMultiphysics.Model object.
        Exception: If the 'settings' input is not a KratosMultiphysics.Parameters object.

    Returns:
        KratosMultiphysics.OutputProcess: An instance of the EnSightOutputProcess.
    """
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object, encapsulating a json string")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return EnSightOutputProcess(model, settings["Parameters"])

class EnSightOutputProcess(KratosMultiphysics.OutputProcess):
    """Manages writing simulation results in the EnSight 6/Gold format.

    This process handles the creation of output files readable by the EnSight
    post-processing software. It uses an OutputController to determine *when*
    to write files and an EnSightOutput object to perform the actual
    I/O operations.
    """

    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        """Initializes the EnSightOutputProcess.

        Args:
            model (KratosMultiphysics.Model): The Kratos model containing the model parts.
            settings (KratosMultiphysics.Parameters): The configuration parameters for the process.
        """
        KratosMultiphysics.OutputProcess.__init__(self)

        model_part_name = settings["model_part_name"].GetString()
        self.model_part = model[model_part_name]

        # This C++ object handles the core I/O operations and validates the settings.
        # Default settings can be found in the corresponding C++ source file "ensight_output.cpp".
        self.ensight_gold_io = KratosMultiphysics.EnSightOutput(self.model_part, settings)

        # Handle folder creation for the output files.
        if settings["save_output_files_in_folder"].GetBool():
            # In a parallel run (MPI), only the main process (rank 0) should create/delete the folder
            # to avoid race conditions.
            if self.model_part.GetCommunicator().MyPID() == 0:
                output_path = settings["output_path"].GetString()
                # If the simulation is not a restart, clean up the output directory.
                if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                    kratos_utils.DeleteDirectoryIfExisting(output_path)
            # A barrier ensures that all processes wait for the main process to finish
            # directory operations before proceeding.
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

        # The OutputController determines the frequency of output (e.g., every N steps, at specific times).
        self.__controller = KratosMultiphysics.OutputController(model, settings)

    def Check(self) -> int:
        """Performs checks to ensure the process is configured correctly.

        This method delegates the check to the internal OutputController.
        It is typically called at the beginning of the simulation.

        Returns:
            int: 0 on success.
        """
        return self.__controller.Check()

    def PrintOutput(self) -> None:
        """Writes the output for the current step to an EnSight file.

        This method is called by the solver at each step where IsOutputStep() returns True.
        It triggers the I/O operation and updates the controller's internal state.
        """
        self.ensight_gold_io.PrintOutput()
        self.__controller.Update()

    def IsOutputStep(self) -> bool:
        """Determines if output should be written at the current simulation step.

        This method queries the OutputController to check if the current time or
        step meets the criteria specified in the settings.

        Returns:
            bool: True if output is required for the current step, False otherwise.
        """
        return self.__controller.Evaluate()