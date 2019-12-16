# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
try:
    import KratosMultiphysics.MappingApplication as KMA
    SpacialMapperFactory = KMA.MapperFactory
except ImportError:
    SpacialMapperFactory = None

# TODO temporary aliases for some variables specific for in plane mapping
# on background mesh
_BACKGROUND_COORDINATE = KM.LOCAL_AXIS_1
# on destination
_MAPPED_BACKGROUND_COORDINATE = KM.LOCAL_AXIS_1
_MAPPED_BACKGROUND_NORMAL = KM.LOCAL_AXIS_2
_OUT_OF_PLANE_DELTA = KM.LOCAL_AXIS_3

class InPlaneVertexMorphingMapper():

    def __init__(self, origin_model_part, destination_model_part, settings, vm_mapper):
        if not SpacialMapperFactory:
            raise Exception("InPlaneVertexMorphingMapper: MappingApplication is required!")

        self.settings = settings
        self.origin_model_part = origin_model_part
        self.destination_model_part = destination_model_part

        # TODO destination does already exist can not add here anymore - non historical??
        # self.destination_model_part.AddNodalSolutionStepVariable(_MAPPED_BACKGROUND_COORDINATE)
        # self.destination_model_part.AddNodalSolutionStepVariable(_MAPPED_BACKGROUND_NORMAL)
        # self.destination_model_part.AddNodalSolutionStepVariable(_OUT_OF_PLANE_DELTA)

        self.vm_mapper = vm_mapper

        self._background_model = KM.Model()
        if self.settings["in_plane_settings"]["model_import_settings"]["input_type"].GetString() in ["mdpa", "vrml", "wrl"]:
            self.background_mesh = self._background_model.CreateModelPart("background_mesh")
            self.background_mesh.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            self.background_mesh.AddNodalSolutionStepVariable(KM.NORMAL)
            self.background_mesh.AddNodalSolutionStepVariable(KSO.NORMALIZED_SURFACE_NORMAL)
            self.background_mesh.AddNodalSolutionStepVariable(_BACKGROUND_COORDINATE)
        else:
            raise Exception("Other model part input options are not yet implemented.")
        self.spacial_mapper = None

    def Initialize(self):
        self.vm_mapper.Initialize()

        if self.settings["in_plane_settings"]["model_import_settings"]["input_type"].GetString() == "mdpa":
            model_part_io = KM.ModelPartIO(self.settings["in_plane_settings"]["model_import_settings"]["input_filename"].GetString())
            model_part_io.ReadModelPart(self.background_mesh)
        elif self.settings["in_plane_settings"]["model_import_settings"]["input_type"].GetString() in ["vrml", "wrl"]:
            model_part_io = WrlIO(self.settings["in_plane_settings"]["model_import_settings"]["input_filename"].GetString())
            model_part_io.ReadModelPart(self.background_mesh)

        background_geometry_utilities = KSO.GeometryUtilities(self.background_mesh)
        background_geometry_utilities.ComputeUnitSurfaceNormals()
        self.spacial_mapper = SpacialMapperFactory.CreateMapper(
            self.background_mesh, self.destination_model_part, self.settings["in_plane_settings"]["spacial_mapper_settings"])
        KSO.MeshControllerUtilities(self.background_mesh).WriteCoordinatesToVariable(_BACKGROUND_COORDINATE)

    def Update(self):
        self.vm_mapper.Update()
        self.spacial_mapper.UpdateInterface()
        self.spacial_mapper.Map(KSO.NORMALIZED_SURFACE_NORMAL, _MAPPED_BACKGROUND_NORMAL)

    def Map(self, origin_variable, destination_variable):
        self.vm_mapper.Map(origin_variable, destination_variable)

        # project on tangent plane
        geometry_utilities = KSO.GeometryUtilities(self.destination_model_part)
        geometry_utilities.ProjectNodalVariableToTangentOfDirection(
            destination_variable, _MAPPED_BACKGROUND_NORMAL)

        self._CorrectOutOfPlanePart(destination_variable)

        # TODO possible improvement if backgroundmesh too small: damp based on distance from background mesh?

    def InverseMap(self, destination_variable, origin_variable):
        # TODO possible improvement if backgroundmesh too small: damp based on distance from background mesh?

        # project on tangent plane
        geometry_utilities = KSO.GeometryUtilities(self.destination_model_part)
        geometry_utilities.ProjectNodalVariableToTangentOfDirection(
            destination_variable, _MAPPED_BACKGROUND_NORMAL)

        self.vm_mapper.InverseMap(destination_variable, origin_variable)

    def _CorrectOutOfPlanePart(self, destination_variable):
        geometry_utilities = KSO.GeometryUtilities(self.destination_model_part)
        mesh_utilities = KSO.MeshControllerUtilities(self.destination_model_part)

        # temp update of destination geometry
        mesh_utilities.UpdateMeshAccordingInputVariable(destination_variable)
        mesh_utilities.SetReferenceMeshToMesh()

        # calculate out of plane delta
        self.spacial_mapper.UpdateInterface()
        self.spacial_mapper.Map(_BACKGROUND_COORDINATE, _MAPPED_BACKGROUND_COORDINATE)
        mesh_utilities.CalculateDistanceToVariable(_MAPPED_BACKGROUND_COORDINATE, _OUT_OF_PLANE_DELTA)
        geometry_utilities.ProjectNodalVariableToDirection(_OUT_OF_PLANE_DELTA, _MAPPED_BACKGROUND_NORMAL)

        mesh_utilities.RevertMeshUpdateAccordingInputVariable(destination_variable)
        mesh_utilities.SetReferenceMeshToMesh()

        mesh_utilities.AddFirstVariableToSecondVariable(_OUT_OF_PLANE_DELTA, destination_variable)

        #TODO damping should only be used if activated for x, y, z - where to check??!!
