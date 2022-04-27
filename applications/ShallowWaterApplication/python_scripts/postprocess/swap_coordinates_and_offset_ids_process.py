# Importing Kratos library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SwapCoordinatesAndOffsetIdsProcess(model, settings["Parameters"])

class SwapCoordinatesAndOffsetIdsProcess(KM.Process):
    def __init__(self, model, settings):
        """ SwapCoordinatesAndOffsetIdsProcess.

        This process provides several tools for post-processing.
        - Swap the YZ coordinates in order to make 2D simulations consistent at post process.
        - Offset the ids in order to differentiate the model parts at the post processing.
        """

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"        : "model_part_name",
                "swap_yz_coordinates"    : true,
                "nodes_ids_offset"       : 0,
                "elements_ids_offset"    : 0,
                "conditions_ids_offset"  : 0,
                "properties_ids_offset"  : 0
            }
        """)
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.swap_yz_coordinates = settings["swap_yz_coordinates"].GetBool()
        self.nodes_ids_offset = settings["nodes_ids_offset"].GetInt()
        self.elements_ids_offset = settings["elements_ids_offset"].GetInt()
        self.conditions_ids_offset = settings["conditions_ids_offset"].GetInt()
        self.properties_ids_offset = settings["properties_ids_offset"].GetInt()

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
        """Swap the mesh and offset the Ids."""
        if self.swap_yz_coordinates:
            self._SwapYZCoordinates()
        self._OffsetIds()


    def ExecuteAfterOutputStep(self):
        """Restore the mesh swapping and the Ids offset."""
        if self.swap_yz_coordinates:
            self._SwapYZCoordinates()
        self._OffsetIds(-1)


    def _SwapYZCoordinates(self):
        SW.ShallowWaterUtilities().SwapYZCoordinates(self.model_part)
        SW.ShallowWaterUtilities().SwapY0Z0Coordinates(self.model_part)


    def _OffsetIds(self, sign=1):
        SW.ShallowWaterUtilities().OffsetIds(self.model_part.Nodes, sign * self.nodes_ids_offset)
        SW.ShallowWaterUtilities().OffsetIds(self.model_part.Elements, sign * self.elements_ids_offset)
        SW.ShallowWaterUtilities().OffsetIds(self.model_part.Conditions, sign * self.conditions_ids_offset)
        SW.ShallowWaterUtilities().OffsetIds(self.model_part.Properties, sign * self.properties_ids_offset)
