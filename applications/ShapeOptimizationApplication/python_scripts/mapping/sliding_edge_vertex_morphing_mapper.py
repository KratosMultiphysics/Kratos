# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Schmoelz David, https://github.com/dschmoelz
#
# ==============================================================================

import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
try:
    import KratosMultiphysics.MappingApplication as KMA
    SpacialMapperFactory = KMA.MapperFactory
except ImportError:
    SpacialMapperFactory = None
from ..custom_ios.wrl_io import WrlIO

class SlidingEdgeVertexMorphingMapper():
    """
    The InPlaneVertexMorphingMapper extends the standard Vertex Morphing approach
    by restricting the shape update of nodes to an in-plane motion only. The nodes
    are only allowed to float on a predefined background mesh.
    The background mesh can be the initial mesh of the design surface, or another mesh
    describing the same geometry (ideally also filling holes in the surface).
    This is especially important if the design surface extends during the optimization.

    Limitations:
    - Damping can only be used if all cartesian directions (x,y,z) are damped.
    - The projection of gradients on the surface normals has to be deactivated.
    """

    def __init__(self, origin_model_part, destination_model_part, settings):
        from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
        if not SpacialMapperFactory:
            raise Exception("SlidingEdgeVertexMorphingMapper: MappingApplication is required!")

        self.settings = settings
        self.origin_model_part = origin_model_part
        self.destination_model_part = destination_model_part

        extracted_vm_settings = settings.Clone()
        extracted_vm_settings["sliding_edge_morphing"].SetBool(False)
        extracted_vm_settings.RemoveValue("sliding_edge_morphing_settings")
        self.vm_mapper = mapper_factory.CreateMapper(origin_model_part, destination_model_part, extracted_vm_settings)

        in_plane_settings = self.settings["sliding_edge_morphing_settings"]
        in_plane_settings.RecursivelyValidateAndAssignDefaults(self.GetDefaultInPlaneSettings())

        self._background_model = KM.Model()
        if in_plane_settings["model_import_settings"]["input_type"].GetString() in ["mdpa", "vrml", "wrl"]:
            background_main_mesh = self._background_model.CreateModelPart("background_mesh")
            background_main_mesh.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            background_main_mesh.AddNodalSolutionStepVariable(KM.NORMAL)
            background_main_mesh.AddNodalSolutionStepVariable(KSO.NORMALIZED_SURFACE_NORMAL)
            background_main_mesh.AddNodalSolutionStepVariable(KSO.BACKGROUND_COORDINATE)
        else:
            raise Exception("Other model part input options are not yet implemented.")

        self.background_mesh = None  # created in Initialize
        self.spacial_mapper = None  # created in Initialize

    @classmethod
    def GetDefaultInPlaneSettings(cls):
        return KM.Parameters("""{
            "model_import_settings" : {
                "input_type" : "mdpa",
                "input_filename" : "UNKNOWN_NAME"
            },
            "sliding_edge_sub_model_part_name" : "",
            "background_sub_model_part_name" : "",
            "spacial_mapper_settings" : {
                "mapper_type": "nearest_element",
                "echo_level" : 0
            }
        }""")

    def GetStructuralSimilaritySetttings(self):
        return KM.Parameters("""{
            "
        }""")

    def Initialize(self):
        self.vm_mapper.Initialize()

        sliding_edge_settings = self.settings["sliding_edge_morphing_settings"]

        background_main_mesh = self._background_model.GetModelPart("background_mesh")
        if sliding_edge_settings["model_import_settings"]["input_type"].GetString() == "mdpa":
            model_part_io = KM.ModelPartIO(sliding_edge_settings["model_import_settings"]["input_filename"].GetString())
            model_part_io.ReadModelPart(background_main_mesh)
        elif sliding_edge_settings["model_import_settings"]["input_type"].GetString() in ["vrml", "wrl"]:
            model_part_io = WrlIO(sliding_edge_settings["model_import_settings"]["input_filename"].GetString())
            model_part_io.ReadModelPart(background_main_mesh)

        background_sub_model_part_name = sliding_edge_settings["background_sub_model_part_name"].GetString()
        if background_sub_model_part_name == "":
            self.background_mesh = background_main_mesh
        else:
            self.background_mesh = background_main_mesh.GetSubModelPart(background_sub_model_part_name)

        sliding_edge_sub_model_part_name = sliding_edge_settings["sliding_edge_sub_model_part_name"].GetString()
        root_model_part = self.destination_model_part.GetRootModelPart()
        self.sliding_edge_mesh = root_model_part.GetSubModelPart(sliding_edge_sub_model_part_name)

        background_geometry_utilities = KSO.GeometryUtilities(self.background_mesh)
        background_geometry_utilities.ComputeUnitSurfaceNormals()
        self.spacial_mapper = SpacialMapperFactory.CreateMapper(
            self.background_mesh, self.sliding_edge_mesh, sliding_edge_settings["spacial_mapper_settings"])

        KSO.MeshControllerUtilities(self.background_mesh).WriteCoordinatesToVariable(KSO.BACKGROUND_COORDINATE)

        #self._structural_mechanics_analysis = StructuralMechanicsAnalysis(model, self.MeshSolverSettings)

    def Update(self):
        self.vm_mapper.Update()
        self.spacial_mapper.UpdateInterface()
        self.spacial_mapper.Map(KSO.NORMALIZED_SURFACE_NORMAL, KSO.BACKGROUND_NORMAL)

        # Model kopieren als Structural model
        
        structural_mechanics_analysis_setings = self.GetStructuralSimilaritySetttings()
        self._structural_mechanics_analysis = StructuralMechanicsAnalysis(model, structural_mechanics_analysis_setings)

    def Map(self, origin_variable, destination_variable):
        self.vm_mapper.Map(origin_variable, destination_variable)

        geometry_utilities = KSO.GeometryUtilities(self.sliding_edge_mesh)
        geometry_utilities.ProjectNodalVariableOnTangentPlane(
            destination_variable, KSO.BACKGROUND_NORMAL)

        # non-linear correction
        self._CorrectOutOfPlanePart(destination_variable)

    def InverseMap(self, destination_variable, origin_variable):
        geometry_utilities = KSO.GeometryUtilities(self.sliding_edge_mesh)
        geometry_utilities.ProjectNodalVariableOnTangentPlane(
            destination_variable, KSO.BACKGROUND_NORMAL)

        self.vm_mapper.InverseMap(destination_variable, origin_variable)

    def _CorrectOutOfPlanePart(self, destination_variable):
        geometry_utilities = KSO.GeometryUtilities(self.sliding_edge_mesh)
        mesh_utilities = KSO.MeshControllerUtilities(self.sliding_edge_mesh)

        mesh_utilities.UpdateMeshAccordingInputVariable(destination_variable)
        mesh_utilities.SetReferenceMeshToMesh()

        self.spacial_mapper.UpdateInterface()
        self.spacial_mapper.Map(KSO.BACKGROUND_COORDINATE, KSO.BACKGROUND_COORDINATE)
        mesh_utilities.SubtractCoordinatesFromVariable(KSO.BACKGROUND_COORDINATE, KSO.OUT_OF_PLANE_DELTA)
        geometry_utilities.ProjectNodalVariableOnDirection(KSO.OUT_OF_PLANE_DELTA, KSO.BACKGROUND_NORMAL)

        mesh_utilities.RevertMeshUpdateAccordingInputVariable(destination_variable)
        mesh_utilities.SetReferenceMeshToMesh()

        mesh_utilities.AddFirstVariableToSecondVariable(KSO.OUT_OF_PLANE_DELTA, destination_variable)

    def _CreateStructuralModel(self):
        
        self._structural_model = KM.Model()
        self._structural_model.CreateSubModelPart("MeshMotionCondtions")
        self._structural_model.CreateSubModelPart("FixedConditions")
