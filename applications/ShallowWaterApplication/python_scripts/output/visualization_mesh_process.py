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
    def __init__(self, model, settings):
        """ VisualizationMeshProcess.

        This process provides several tools for post-processing.
        - Generation of an auxiliary model part for the topography visualization as a separate file.
        - Setting the TOPOGRAPHY and FREE_SURFACE_ELEVATION into DISPLACEMENT_Z or Z-coordinate in order to view the mesh deformation.
        - Saving the HEIGHT and FREE_SURFACE_ELEVATION as non-historical only on wet nodes. Dry nodes are set with the no-data value of GiD.
        """

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"                : "model_part_name",
                "topographic_model_part_name"    : "topographic_model_part",
                "create_topographic_model_part"  : true,
                "mesh_deformation_mode"          : "use_nodal_displacement",
                "nodal_variables_to_transfer"    : ["TOPOGRAPHY"],
                "nonhistorical_variables_to_transfer" : []
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)

        self.computing_model_part = model[settings["model_part_name"].GetString()]

        # Getting the deformation mode options
        self.deform_mesh = self._GetDeformMeshFlag(settings["mesh_deformation_mode"].GetString())

        # Creating the topographic model part if specified
        self.topographic_model_part = None
        if settings["create_topographic_model_part"].GetBool():
            self.topographic_model_part = model.CreateModelPart(settings["topographic_model_part_name"].GetString())
        else:
            if model.Has(settings["topographic_model_part_name"].GetString()):
                self.topographic_model_part = model.GetModelPart(settings["topographic_model_part_name"].GetString())

        # Creating the variables list
        self.nodal_variables = GenerateVariableListFromInput(settings["nodal_variables_to_transfer"])
        self.nonhistorical_variables = GenerateVariableListFromInput(settings["nonhistorical_variables_to_transfer"])


    def ExecuteInitialize(self):
        if self.topographic_model_part is not None:
            self._DuplicateModelPart()


    def ExecuteBeforeSolutionLoop(self):
        # Set the results only over the wet domain as non-historical
        self._StoreNonHistoricalVariablesGiDNoDataIfDry()

        # Transferring the nodal variables
        if self.topographic_model_part is not None:
            self._TransferVariables()

        # Deform the mesh according to the input options
        if self.deform_mesh:
            self._DeformMesh()
        else:
            self._InitializeDisplacement()
            self._UpdateDisplacement()


    def ExecuteBeforeOutputStep(self):
        # Set the results only over the wet domain as non-historical
        self._StoreNonHistoricalVariablesGiDNoDataIfDry()

        # Transferring the nodal variables
        if self.topographic_model_part is not None:
            self._TransferVariables()

        # Deform the mesh according to the input options
        if self.deform_mesh:
            self._DeformMesh()
        else:
            self._UpdateDisplacement()


    def ExecuteAfterOutputStep(self):
        # Restoring the mesh
        if self.deform_mesh:
            self._RestoreMesh()


    def _StoreNonHistoricalVariablesGiDNoDataIfDry(self):
        SW.ShallowWaterUtilities().StoreNonHistoricalGiDNoDataIfDry(self.computing_model_part, SW.HEIGHT)
        SW.ShallowWaterUtilities().StoreNonHistoricalGiDNoDataIfDry(self.computing_model_part, SW.FREE_SURFACE_ELEVATION)


    def _DuplicateModelPart(self):
        KM.MergeVariableListsUtility().Merge(self.computing_model_part, self.topographic_model_part)
        element_num_nodes = len(self.computing_model_part.Elements.__iter__().__next__().GetNodes())
        condition_num_nodes = len(self.computing_model_part.Conditions.__iter__().__next__().GetNodes())
        reference_element = "Element2D{}N".format(element_num_nodes)
        reference_condition = "LineCondition2D{}N".format(condition_num_nodes)
        KM.DuplicateMeshModeler(self.computing_model_part).GenerateMesh(
            self.topographic_model_part, reference_element, reference_condition)
        SW.ShallowWaterUtilities().OffsetIds(self.topographic_model_part.Nodes)
        SW.ShallowWaterUtilities().OffsetIds(self.topographic_model_part.Elements)
        SW.ShallowWaterUtilities().OffsetIds(self.topographic_model_part.Conditions)


    def _TransferVariables(self):
        for variable in self.nodal_variables:
            KM.VariableUtils().CopyModelPartNodalVar(
                variable,
                self.computing_model_part,
                self.topographic_model_part,
                0)
        for variable in self.nonhistorical_variables:
            KM.VariableUtils().CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar(
                variable,
                self.computing_model_part,
                self.topographic_model_part,
                KM.Flags())


    def _DeformMesh(self):
        SW.ShallowWaterUtilities().SetMeshZCoordinate(
            self.computing_model_part,
            SW.FREE_SURFACE_ELEVATION)
        if self.topographic_model_part is not None:
            SW.ShallowWaterUtilities().SetMeshZCoordinate(
                self.topographic_model_part,
                SW.TOPOGRAPHY)


    def _RestoreMesh(self):
        SW.ShallowWaterUtilities().SetMeshZCoordinateToZero(
            self.computing_model_part)
        if self.topographic_model_part is not None:
            SW.ShallowWaterUtilities().SetMeshZCoordinateToZero(
                self.topographic_model_part)


    def _InitializeDisplacement(self):
        KM.VariableUtils().SetNonHistoricalVariableToZero(
            KM.DISPLACEMENT,
            self.computing_model_part.Nodes)
        if self.topographic_model_part is not None:
            KM.VariableUtils().SetNonHistoricalVariableToZero(
                KM.DISPLACEMENT,
                self.topographic_model_part.Nodes)


    def _UpdateDisplacement(self):
        KM.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(
            SW.FREE_SURFACE_ELEVATION,
            KM.DISPLACEMENT_Z,
            self.computing_model_part,
            self.computing_model_part,
            0)
        if self.topographic_model_part is not None:
            KM.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(
                SW.TOPOGRAPHY,
                KM.DISPLACEMENT_Z,
                self.topographic_model_part,
                self.topographic_model_part,
                0)


    def _GetDeformMeshFlag(self, mesh_deformation_mode):
        if mesh_deformation_mode == "use_z_coordinate":
            return True
        elif mesh_deformation_mode == "use_nodal_displacement":
            return False
        else:
            msg = """VisualizationMeshProcess.
            Unknown 'mesh_deformation_mode' = '{}'. The possible options are:
             - use_z_coordinate
             - use_nodal_displacement
            """
            raise Exception(msg.format(mesh_deformation_mode))
