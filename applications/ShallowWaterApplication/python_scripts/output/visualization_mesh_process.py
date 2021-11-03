# Importing Kratos library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

# Importing useful utilities
from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VisualizationMeshProcess(Model, settings["Parameters"])

class VisualizationMeshProcess(KM.Process):
    def __init__(self, Model, settings):
        """ VisualizationMeshProcess.

        This process provides several tools for post-processing.
        - Generation of an auxiliary model part for the topography visualization as a separate file.
        - Setting the TOPOGRAPHY and FREE_SURFACE_ELEVATION into DISPLACEMENT_Z or Z-coordinate to view the mesh deformation.
        - Duplication of the properties to use them as a dry-wet flag.
        """

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"                : "model_part_name",
                "topographic_model_part_name"    : "topographic_model_part",
                "create_topographic_model_part"  : true,
                "mesh_deformation_mode"          : "use_nodal_displacement",
                "topography_variable"            : "TOPOGRAPHY",
                "free_surface_variable"          : "FREE_SURFACE_ELEVATION",
                "nodal_variables_to_transfer"    : ["TOPOGRAPHY"],
                "nonhistorical_variables_to_transfer" : []
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)

        self.computing_model_part = Model[settings["model_part_name"].GetString()]

        # Getting the deformation mode options and storing the variable for the free surface
        mesh_deformation_mode = settings["mesh_deformation_mode"].GetString()
        if mesh_deformation_mode == "use_z_coordinate":
            self.deform_mesh = True
        elif mesh_deformation_mode =="use_nodal_displacement":
            self.deform_mesh = False
        else:
            msg = """VisualizationMeshProcess.
            Unkown 'mesh_deformation_mode'. The possible options are:
             - use_z_coordinate
             - use_nodal_displacement
            The input is {}
            """
            raise Exception(msg.format(mesh_deformation_mode))
        self.free_surface_variable = KM.KratosGlobals.GetVariable(settings["free_surface_variable"].GetString())

        # Creating the utility for the topographic model part management
        self.use_topographic_model_part = settings["create_topographic_model_part"].GetBool()
        if self.use_topographic_model_part:
            self.topographic_model_part = Model.CreateModelPart(settings["topographic_model_part_name"].GetString())
            self.topography_utility = SW.ReplicateModelPartUtility(self.computing_model_part, self.topographic_model_part)
            self.nodal_variables = GenerateVariableListFromInput(settings["nodal_variables_to_transfer"])
            self.nonhistorical_variables = GenerateVariableListFromInput(settings["nonhistorical_variables_to_transfer"])
            self.topography_variable = KM.KratosGlobals.GetVariable(settings["topography_variable"].GetString())

    def ExecuteInitialize(self):
        if self.use_topographic_model_part:
            self.topography_utility.Replicate()

    def ExecuteBeforeSolutionLoop(self):
        # Set the results only over the wet domain as non-historical
        self._StoreNonHistoricalVariablesGiDNoDataIfDry()

        # Transferring the nodal variables
        if self.use_topographic_model_part:
            for variable in self.nodal_variables:
                self.topography_utility.TransferVariable(variable)
            for variable in self.nonhistorical_variables:
                self.topography_utility.TransferNonHistoricalVariable(variable)

        # Moving the mesh according to the specified options
        if self.deform_mesh:
            SW.ShallowWaterUtilities().SetMeshZCoordinate(self.computing_model_part, self.free_surface_variable)
            if self.use_topographic_model_part:
                SW.ShallowWaterUtilities().SetMeshZCoordinate(self.topographic_model_part, self.topography_variable)
        else:
            SW.ShallowWaterUtilities().SetMeshZCoordinateToZero(self.computing_model_part)
            KM.VariableUtils().SetNonHistoricalVariableToZero(KM.DISPLACEMENT, self.computing_model_part.Nodes)
            if self.use_topographic_model_part:
                SW.ShallowWaterUtilities().SetMeshZCoordinateToZero(self.topographic_model_part)
                KM.VariableUtils().SetNonHistoricalVariableToZero(KM.DISPLACEMENT, self.topographic_model_part.Nodes)

    def ExecuteBeforeOutputStep(self):
        # Set the results only over the wet domain as non-historical
        self._StoreNonHistoricalVariablesGiDNoDataIfDry()

        # Transferring the nodal variables
        if self.use_topographic_model_part:
            for variable in self.nodal_variables:
                self.topography_utility.TransferVariable(variable)
            for variable in self.nonhistorical_variables:
                self.topography_utility.TransferNonHistoricalVariable(variable)

        # Moving the mesh according to the specified options
        if self.deform_mesh:
            self.topography_utility.SetOriginMeshZCoordinate(self.free_surface_variable)
            if self.use_topographic_model_part:
                self.topography_utility.SetDestinationMeshZCoordinate(self.topography_variable)
        else:
            KM.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(
                self.free_surface_variable,
                KM.DISPLACEMENT_Z,
                self.computing_model_part,
                self.computing_model_part,
                0)
            if self.use_topographic_model_part:
                KM.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(
                    self.topography_variable,
                    KM.DISPLACEMENT_Z,
                    self.topographic_model_part,
                    self.topographic_model_part,
                    0)

    def ExecuteAfterOutputStep(self):
        # Restoring the mesh
        if self.deform_mesh:
            SW.ShallowWaterUtilities().SetMeshZCoordinateToZero(self.computing_model_part)

    def _StoreNonHistoricalVariablesGiDNoDataIfDry(self):
        SW.ShallowWaterUtilities().StoreNonHistoricalGiDNoDataIfDry(self.computing_model_part, SW.HEIGHT)
        SW.ShallowWaterUtilities().StoreNonHistoricalGiDNoDataIfDry(self.computing_model_part, SW.FREE_SURFACE_ELEVATION)
