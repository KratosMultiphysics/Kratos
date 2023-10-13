# Importing Kratos library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SwapCoordinatesProcess(model, settings["Parameters"])

class SwapCoordinatesProcess(KM.Process):

    """ SwapCoordinatesProcess.

    This process swaps the YZ coordinates of a model part in order to merge several
    2D simulations for post-processing purpose.
    """

    def __init__(self, model, settings):
        """Constructor of the class."""
        KM.Process.__init__(self)

        default_settings = KM.Parameters("""{
                "model_part_name"        : "model_part_name"
            }""")

        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[settings["model_part_name"].GetString()]
        self.execute_initialize_solution_step_is_called = False


    def ExecuteBeforeSolutionLoop(self):
        """Perform the transformation before printing the initial mesh."""
        self.ExecuteBeforeOutputStep()


    def ExecuteInitializeSolutionStep(self):
        """Undo the transformation after the initial mesh is printed."""
        if not self.execute_initialize_solution_step_is_called:
            self.ExecuteAfterOutputStep()
            self.execute_initialize_solution_step_is_called = True


    def ExecuteBeforeOutputStep(self):
        """Swap the mesh."""
        self._SwapYZCoordinates()


    def ExecuteAfterOutputStep(self):
        """Restore the mesh swapping."""
        self._SwapYZCoordinates()


    def _SwapYZCoordinates(self):
        SW.ShallowWaterUtilities().SwapYZCoordinates(self.model_part)
        SW.ShallowWaterUtilities().SwapY0Z0Coordinates(self.model_part)
