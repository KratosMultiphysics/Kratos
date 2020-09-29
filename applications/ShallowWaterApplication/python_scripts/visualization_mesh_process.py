import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VisualizationMeshProcess(Model, settings["Parameters"])

## This process adds a topography model part and mark the wet and dry domain
class VisualizationMeshProcess(KM.Process):

    def __init__(self, Model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"              : "model_part_name",
                "topographic_model_part_name"  : "topographic_model_part",
                "create_topographic_model_part": true,
                "update_topography"            : false,
                "update_free_surface"          : false,
                "topography_variable"          : "BATHYMETRY",
                "free_surface_variable"        : "FREE_SURFACE_ELEVATION",
                "nodal_variables_to_transfer"  : [],
                "nonhistorical_variables_to_transfer" : []
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)

        self.computing_model_part = Model[settings["model_part_name"].GetString()]

        # Creating the utility for the topographic model part management
        self.use_topographic_model_part = settings["create_topographic_model_part"].GetBool()
        self.update_topography = False
        self.update_free_surface = False
        if self.use_topographic_model_part:
            self.topographic_model_part = Model.CreateModelPart(settings["topographic_model_part_name"].GetString())
            self.topography_utility = SW.ReplicateModelPartUtility(self.computing_model_part, self.topographic_model_part)
            self.nodal_variables = self._GenerateVariableListFromInput(settings["nodal_variables_to_transfer"])
            self.nonhistorical_variables = self._GenerateVariableListFromInput(settings["nonhistorical_variables_to_transfer"])
            self.topography_variable = KM.KratosGlobals.GetVariable(settings["topography_variable"].GetString())
            self.free_surface_variable = KM.KratosGlobals.GetVariable(settings["free_surface_variable"].GetString())
            self.update_topography = settings["update_topography"].GetBool()
            self.update_free_surface = settings["update_free_surface"].GetBool()

        # The DefineAuxiliaryProperties method duplicates the current number of properties:
        # For each property, it creates another one, which means dry state.
        # It should be called only once, otherwise, the number of properties will increase without limits
        self.properties_utility = SW.PostProcessUtilities(self.computing_model_part)
        self.properties_utility.DefineAuxiliaryProperties()

        KM.VariableUtils().SetNonHistoricalVariableToZero(SW.WATER_HEIGHT, self.computing_model_part.Nodes)
        KM.VariableUtils().SetNonHistoricalVariableToZero(SW.WATER_SURFACE, self.computing_model_part.Nodes)

    def ExecuteInitialize(self):
        if self.use_topographic_model_part:
            self.topography_utility.Replicate()

    def ExecuteBeforeSolutionLoop(self):
        SW.ShallowWaterUtilities().IdentifyWetDomain(self.computing_model_part, KM.FLUID, 0.0)
        SW.ShallowWaterUtilities().ComputeVisualizationWaterSurface(self.computing_model_part)
        SW.ShallowWaterUtilities().ComputeVisualizationWaterHeight(self.computing_model_part, KM.FLUID, 0.0)
        if self.use_topographic_model_part:
            for variable in self.nodal_variables:
                self.topography_utility.TransferVariable(variable)
            for variable in self.nonhistorical_variables:
                self.topography_utility.TransferNonHistoricalVariable(variable)
            self.topography_utility.SetDestinationMeshZCoordinate(self.topography_variable)
            if self.update_free_surface:
                self.topography_utility.SetOriginMeshZCoordinate(self.free_surface_variable)
            else:
                self.topography_utility.SetOriginMeshZCoordinate()

    def ExecuteBeforeOutputStep(self):
        # The elements should be active to be included in the GidOutputProcess
        KM.VariableUtils().SetFlag(KM.ACTIVE, True, self.computing_model_part.Elements)
        SW.ShallowWaterUtilities().ComputeFreeSurfaceElevation(self.computing_model_part)
        SW.ShallowWaterUtilities().ComputeVisualizationWaterSurface(self.computing_model_part)
        SW.ShallowWaterUtilities().ComputeVisualizationWaterHeight(self.computing_model_part, KM.FLUID, 0.0)
        self.properties_utility.AssignDryWetProperties(KM.FLUID)

        if self.use_topographic_model_part:
            for variable in self.nodal_variables:
                self.topography_utility.TransferVariable(variable)
            for variable in self.nonhistorical_variables:
                self.topography_utility.TransferNonHistoricalVariable(variable)
            if self.update_topography:
                self.topography_utility.SetDestinationMeshZCoordinate(self.topography_variable)
            if self.update_free_surface:
                self.topography_utility.SetOriginMeshZCoordinate(self.free_surface_variable)

    def ExecuteAfterOutputStep(self):
        self.properties_utility.RestoreDryWetProperties()
        if self.use_topographic_model_part and self.update_free_surface:
            self.topography_utility.SetOriginMeshZCoordinate()

    def _GenerateVariableListFromInput(self, variables_array):
        '''Parse a list of variables from input.'''
        # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
        return [ KM.KratosGlobals.GetVariable(variables_array[i].GetString()) for i in range(variables_array.size()) ]
