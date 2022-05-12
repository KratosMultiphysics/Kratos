# Importing Kratos library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

# Importing useful utilities
from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VisualizationMeshProcess(model, settings["Parameters"])

class VisualizationMeshProcess(KM.Process):
    """ VisualizationMeshProcess.

    This process provides several tools for post-processing.
    - Generation of an auxiliary model part for the topography visualization as a separate file.
    - Setting the TOPOGRAPHY and FREE_SURFACE_ELEVATION into DISPLACEMENT_Z or Z-coordinate in order to view the mesh deformation.
    - Saving the HEIGHT and FREE_SURFACE_ELEVATION as non-historical only on wet nodes. Dry nodes are set with the no-data value of GiD.
    """

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""
            {
                "model_part_name"                : "model_part_name",
                "topographic_model_part_name"    : "topographic_model_part",
                "create_topographic_model_part"  : true,
                "free_surface_deformation_mode"  : "nodal_displacement",
                "topography_deformation_mode"    : "z_coordinate",
                "mean_water_level"               : 0.0,
                "nodal_variables_to_transfer"    : ["TOPOGRAPHY"],
                "nonhistorical_variables_to_transfer" : []
            }
            """)

    _mesh_deformation_modes = {
        "z_coordinate"       : True,
        "nodal_displacement" : False
    }

    def __init__(self, model, settings):
        """Constructor with Model and Parameters."""

        KM.Process.__init__(self)
        self.model = model
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.computing_model_part = self.model[settings["model_part_name"].GetString()]

        # Get the deformation mode options
        self.deform_free_surface = self._GetDeformMeshFlag(settings["free_surface_deformation_mode"])
        self.deform_topography = self._GetDeformMeshFlag(settings["topography_deformation_mode"])
        self.mean_water_level = settings["mean_water_level"].GetDouble()

        # Creating the topographic model part if specified
        self.topographic_model_part = None
        if settings["create_topographic_model_part"].GetBool():
            self.topographic_model_part = self.model.CreateModelPart(settings["topographic_model_part_name"].GetString())
        else:
            if self.model.HasModelPart(settings["topographic_model_part_name"].GetString()):
                self.topographic_model_part = self.model.GetModelPart(settings["topographic_model_part_name"].GetString())

        # Creating the variables list
        self.nodal_variables = GenerateVariableListFromInput(settings["nodal_variables_to_transfer"])
        self.nonhistorical_variables = GenerateVariableListFromInput(settings["nonhistorical_variables_to_transfer"])


    def ExecuteInitialize(self):
        """Generate the topographic model part if specified or it already exists."""
        if self.topographic_model_part is not None:
            self._DuplicateModelPart()
        KM.VariableUtils().SetNonHistoricalVariable(KM.TEMPERATURE, 2, self.computing_model_part.Nodes)


    def ExecuteBeforeSolutionLoop(self):
        """Initialize the visualization operations."""

        # Deform the mesh according to the input options
        if not self.deform_free_surface:
            self._FlattenMeshCoordinates(self.computing_model_part)
            self._InitializeDisplacement(self.computing_model_part)

        # Deform the topography and transfer the nodal variables
        if self.topographic_model_part is not None:
            if not self.deform_topography:
                self._FlattenMeshCoordinates(self.topographic_model_part)
                self._InitializeDisplacement(self.topographic_model_part)

        self.ExecuteBeforeOutputStep()


    def ExecuteBeforeOutputStep(self):
        """Perform the visualization operations."""
        # Set the results only over the wet domain as non-historical
        self._StoreNonHistoricalVariablesGiDNoDataIfDry()

        # Deform the mesh according to the input options
        if self.deform_free_surface:
            self._DeformMesh(self.computing_model_part, SW.FREE_SURFACE_ELEVATION)
        else:
            self._UpdateDisplacement(self.computing_model_part, SW.FREE_SURFACE_ELEVATION)

        # Deform the topography and transfer the nodal variables
        if self.topographic_model_part is not None:
            self._TransferVariables()

            if self.deform_topography:
                self._DeformMesh(self.topographic_model_part, SW.TOPOGRAPHY)
            else:
                self._UpdateDisplacement(self.topographic_model_part, SW.TOPOGRAPHY)


    def ExecuteAfterOutputStep(self):
        """Restore the mesh deformation."""
        if self.deform_free_surface:
            self._RestoreMesh(self.computing_model_part)
        
        if self.topographic_model_part is not None:
            if self.deform_topography:
                self._RestoreMesh(self.topographic_model_part)


    def _StoreNonHistoricalVariablesGiDNoDataIfDry(self):
        SW.ShallowWaterUtilities().StoreNonHistoricalGiDNoDataIfDry(self.computing_model_part, SW.HEIGHT)
        SW.ShallowWaterUtilities().StoreNonHistoricalGiDNoDataIfDry(self.computing_model_part, SW.FREE_SURFACE_ELEVATION)


    def _DuplicateModelPart(self):
        KM.MergeVariableListsUtility().Merge(self.computing_model_part, self.topographic_model_part)
        self.topographic_model_part.ProcessInfo = self.computing_model_part.ProcessInfo
        element_num_nodes = len(self.computing_model_part.Elements.__iter__().__next__().GetNodes())
        condition_num_nodes = len(self.computing_model_part.Conditions.__iter__().__next__().GetNodes())
        reference_element = "Element2D{}N".format(element_num_nodes)
        reference_condition = "LineCondition2D{}N".format(condition_num_nodes)
        KM.DuplicateMeshModeler(self.computing_model_part).GenerateMesh(
            self.topographic_model_part, reference_element, reference_condition)
        KM.CopyPropertiesModeler(self.computing_model_part, self.topographic_model_part).SetupModelPart()
        SW.ShallowWaterUtilities().OffsetIds(self.topographic_model_part.Nodes)
        SW.ShallowWaterUtilities().OffsetIds(self.topographic_model_part.Elements)
        SW.ShallowWaterUtilities().OffsetIds(self.topographic_model_part.Conditions)
        SW.ShallowWaterUtilities().OffsetIds(self.topographic_model_part.Properties)


    def _TransferVariables(self):
        for variable in self.nodal_variables:
            KM.VariableUtils().CopyModelPartNodalVar(
                variable,
                self.computing_model_part,
                self.topographic_model_part,
                0)
        for variable in self.nonhistorical_variables:
            KM.VariableUtils().CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar(
                variable, variable,
                self.computing_model_part,
                self.topographic_model_part,
                KM.Flags(), False)


    def _FlattenMeshCoordinates(self, model_part):
        SW.ShallowWaterUtilities().SetMeshZCoordinateToZero(model_part)
        SW.ShallowWaterUtilities().SetMeshZ0CoordinateToZero(model_part)
        SW.ShallowWaterUtilities().OffsetMeshZCoordinate(model_part, self.mean_water_level)


    def _DeformMesh(self, model_part, variable):
        SW.ShallowWaterUtilities().SetMeshZCoordinate(model_part, variable)
        SW.ShallowWaterUtilities().OffsetMeshZCoordinate(model_part, self.mean_water_level)


    def _RestoreMesh(self, model_part):
        SW.ShallowWaterUtilities().SetMeshZCoordinateToZero(model_part)


    @staticmethod
    def _InitializeDisplacement(model_part):
        KM.VariableUtils().SetNonHistoricalVariableToZero(KM.DISPLACEMENT, model_part.Nodes)


    @staticmethod
    def _UpdateDisplacement(model_part, variable):
        KM.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(variable, KM.DISPLACEMENT_Z,
                                                                   model_part, model_part, 0)


    def _GetDeformMeshFlag(self, mesh_deformation_mode):
        try:
            value = self._mesh_deformation_modes[mesh_deformation_mode.GetString()]
        except KeyError:
            msg = "VisualizationMeshProcess. Unknown deformation mode '{}'. The possible options are: \n".format(mesh_deformation_mode.GetString())
            for key in self._mesh_deformation_modes.keys():
                msg += " - {}\n".format(key)
            raise Exception(msg)
        return value
