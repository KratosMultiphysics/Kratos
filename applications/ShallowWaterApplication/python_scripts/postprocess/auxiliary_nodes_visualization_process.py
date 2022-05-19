import KratosMultiphysics as KM
from KratosMultiphysics.point_output_process import Interpolate

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AuxiliaryNodesVisualizationProcess(model, settings["Parameters"])

class AuxiliaryNodesVisualizationProcess(KM.Process):
    """This class creates some auxiliary nodes.

    If a reference model part is is specified, a mapping structure is created in
    order to update the solution step data from the reference the auxiliary nodes.
    The reference model part can be a specified by name or as the root model part.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return KM.Parameters('''{
            "model_part_name"           : "",
            "reference_model_part_name" : "",
            "nodes_coordinates"         : [],
            "search_configuration"      : "initial",
            "search_tolerance"          : 1e-6
        }''')

    _search_configurations = {
        "initial" : KM.Configuration.Initial,
        "current" : KM.Configuration.Current
    }

    def __init__(self, model, settings):
        """Constructor of the process."""
        super().__init__()

        # Validate settings
        settings.AddMissingParameters(self.GetDefaultParameters())
        self.settings = settings

        # Retrieve the model part
        self.model_part = model.CreateModelPart(self.settings["model_part_name"].GetString())
        self.reference_model_part = None
        if self.settings["reference_model_part_name"].GetString():
            self.reference_model_part = model.CreateModelPart(self.settings["reference_model_part_name"].GetString())
            self.model_part.ProcessInfo = self.reference_model_part.ProcessInfo
            KM.MergeVariableListsUtility().Merge(self.reference_model_part, self.model_part)
        elif self.model_part.IsSubModelPart():
            self.reference_model_part = self.model_part.GetRootModelPart()

    def ExecuteInitialize(self):
        """Create the auxiliary nodes and set up the mapping structure."""
        # Initialize the variables
        if self.reference_model_part is not None:
            self.node_id = self.reference_model_part.NumberOfNodes()
            variables_names = self.model_part.GetHistoricalVariablesNames()
            self.variables = [KM.KratosGlobals.GetVariable(name) for name in variables_names]
        else:
            self.node_id = 0

        # Create the auxiliary geometries
        for coord in self.settings["nodes_coordinates"]:
            self._AddNode(coord)

        # Initialize the mapping
        if self.reference_model_part is not None:
            self._SearchNodes()

    def ExecuteBeforeOutputStep(self):
        for node, entity, area_coords in zip(self.found_positions, self.entities, self.area_coords):
            for variable in self.variables:
                value = Interpolate(variable, entity, area_coords, historical_value=True)
                node.SetSolutionStepValue(variable, value)

    def _AddNode(self, params):
        self.node_id += 1
        coordinates = params.GetVector()
        self.model_part.CreateNewNode(self.node_id, coordinates[0], coordinates[1], coordinates[2])

    def _SearchNodes(self):
        search_configuration = self._search_configurations[self.settings["search_configuration"].GetString()]
        search_tolerance = self.settings["search_tolerance"].GetDouble()

        self.entities = []
        self.area_coords = []
        self.found_positions = []

        for node in self.model_part.Nodes:
            sf_values = KM.Vector()
            found_id = KM.BruteForcePointLocator(self.reference_model_part).FindElement(node, sf_values, search_configuration, search_tolerance)
            if found_id > -1:
                self.entities.append(self.reference_model_part.Elements[found_id])
                self.area_coords.append(sf_values)
                self.found_positions.append(node)
