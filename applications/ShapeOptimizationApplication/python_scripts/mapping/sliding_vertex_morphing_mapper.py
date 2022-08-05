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
try:
    import KratosMultiphysics.StructuralMechanicsApplication as KSM
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
except ImportError:
    KSM = None
    StructuralMechanicsAnalysis = None
SpacialMapperFactory = KM.MapperFactory

from ..custom_ios.wrl_io import WrlIO

class SlidingVertexMorphingMapper():
    """
    The SlidingVertexMorphingMapper extends the standard Vertex Morphing approach
    by restricting the nodal shape update of a specified model part. The nodes of this model
    part have to lie on a predefined background mesh and are only allowed to float/slide on it.
    SlidingVertexMorphingMapper uses similar techniques as direction damping and
    InPlaneVertexMorphingMapper.

    Limitations:
    - Damping can only be used if all cartesian directions (x,y,z) are damped.
    - The projection of gradients on the surface normals has to be deactivated.
    """

    def __init__(self, origin_model_part, destination_model_part, settings):
        from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
        if not SpacialMapperFactory:
            raise Exception("SlidingVertexMorphingMapper: MappingApplication is required!")
        if not KSM:
            raise Exception("SlidingVertexMorphingMapper: StructuralMechanicsApplication is required!")

        self.settings = settings
        self.origin_model_part = origin_model_part
        self.destination_model_part = destination_model_part

        extracted_vm_settings = settings.Clone()
        extracted_vm_settings["sliding_morphing"].SetBool(False)
        extracted_vm_settings.RemoveValue("sliding_morphing_settings")
        self.vm_mapper = mapper_factory.CreateMapper(origin_model_part, destination_model_part, extracted_vm_settings)

        in_plane_settings = self.settings["sliding_morphing_settings"]
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
        self.spacial_mapper = None   # created in Initialize
        self.sliding_mesh = None     # created in Initialize

    @classmethod
    def GetDefaultInPlaneSettings(cls):
        return KM.Parameters("""{
            "model_import_settings" : {
                "input_type" : "mdpa",
                "input_filename" : "UNKNOWN_NAME"
            },
            "sliding_sub_model_part_name" : "",
            "background_sub_model_part_name" : "",
            "spacial_mapper_settings" : {
                "mapper_type": "nearest_element",
                "echo_level" : 0
            }
        }""")

    def GetStructuralSimilaritySettings(self):
        return KM.Parameters(
        """{
            "problem_data"    : {
                "problem_name"    : "structural_similarity",
                "parallel_type"   : "OpenMP",
                "start_time"      : 0.0,
                "end_time"        : 1.0,
                "echo_level"      : 0
            },
            "solver_settings" : {
                "solver_type"     : "static",
                "analysis_type"   : "linear",
                "model_part_name" : "main",
                "domain_size"     : 3,
                "time_stepping" : {
                    "time_step"       : 1.0
                },
                "model_import_settings"              : {
                    "input_type"     : "use_input_model_part"
                },
                "rotation_dofs"                      : true
            },
            "processes" : {
                "constraints_process_list" : [
                {
                    "python_module" : "assign_vector_variable_process",
                    "kratos_module" : "KratosMultiphysics",
                    "help"          : "This process fixes the selected components of a given vector variable",
                    "process_name"  : "AssignVectorVariableProcess",
                    "Parameters"    : {
                        "mesh_id"         : 0,
                        "model_part_name" : "main.fixed_bc",
                        "variable_name"   : "DISPLACEMENT",
                        "value"           : [0.0,0.0,0.0]
                    }
                }
                ]
            }
        }""")

    def Initialize(self):
        self.vm_mapper.Initialize()

        sliding_settings = self.settings["sliding_morphing_settings"]

        background_main_mesh = self._background_model.GetModelPart("background_mesh")
        if sliding_settings["model_import_settings"]["input_type"].GetString() == "mdpa":
            model_part_io = KM.ModelPartIO(sliding_settings["model_import_settings"]["input_filename"].GetString())
            model_part_io.ReadModelPart(background_main_mesh)
        elif sliding_settings["model_import_settings"]["input_type"].GetString() in ["vrml", "wrl"]:
            model_part_io = WrlIO(sliding_settings["model_import_settings"]["input_filename"].GetString())
            model_part_io.ReadModelPart(background_main_mesh)

        background_sub_model_part_name = sliding_settings["background_sub_model_part_name"].GetString()
        if background_sub_model_part_name == "":
            self.background_mesh = background_main_mesh
        else:
            self.background_mesh = background_main_mesh.GetSubModelPart(background_sub_model_part_name)

        sliding_sub_model_part_name = sliding_settings["sliding_sub_model_part_name"].GetString()
        root_model_part = self.destination_model_part.GetRootModelPart()
        self.sliding_mesh = root_model_part.GetSubModelPart(sliding_sub_model_part_name)

        background_geometry_utilities = KSO.GeometryUtilities(self.background_mesh)
        background_geometry_utilities.ComputeUnitSurfaceNormals()
        self.spacial_mapper = SpacialMapperFactory.CreateMapper(
            self.background_mesh, self.sliding_mesh, sliding_settings["spacial_mapper_settings"])

        KSO.MeshControllerUtilities(self.background_mesh).WriteCoordinatesToVariable(KSO.BACKGROUND_COORDINATE)

    def Update(self):
        self.vm_mapper.Update()
        self.spacial_mapper.UpdateInterface()
        self.spacial_mapper.Map(KSO.NORMALIZED_SURFACE_NORMAL, KSO.BACKGROUND_NORMAL)

        # create structural similarity model and analysis
        self._structural_model = KM.Model()
        self._CreateStructuralModel(self._structural_model)

        self._structural_model_correction = KM.Model()
        self._CreateStructuralModel(self._structural_model_correction)

        structural_mechanics_analysis_setings = self.GetStructuralSimilaritySettings()
        # two structural similarity analysis for damping vectors and non-linear correction
        self._structural_mechanics_analysis = StructuralMechanicsAnalysis(self._structural_model, structural_mechanics_analysis_setings)
        self._structural_mechanics_analysis_correction = StructuralMechanicsAnalysis(self._structural_model_correction, structural_mechanics_analysis_setings)

        # propagate background normal field from sliding mesh to destination
        self._PerformStructuralSimilarityAnalysis(self._structural_mechanics_analysis, self._structural_model, KSO.BACKGROUND_NORMAL)

    def Map(self, origin_variable, destination_variable):
        self.vm_mapper.Map(origin_variable, destination_variable)

        geometry_utilities = KSO.GeometryUtilities(self.destination_model_part)

        geometry_utilities.ProjectNodalVariableOnTangentPlane(
            destination_variable, KSO.BACKGROUND_NORMAL)

        # non-linear correction on sliding mesh
        self._CorrectOutOfPlanePart(destination_variable)

        mesh_utilities_struct =  KSO.MeshControllerUtilities(self._structural_model_correction.GetModelPart("main"))

        # map linear shape update to structural model
        ne_mapper_settings = KM.Parameters("""{
                                "mapper_type": "nearest_element",
                                "echo_level" : 0
                            }""")
        model_mapper = SpacialMapperFactory.CreateMapper(
            self.destination_model_part, self._structural_model_correction.GetModelPart("main"),
            ne_mapper_settings)

        model_mapper.Map(destination_variable, destination_variable)

        # update structural model
        mesh_utilities_struct.UpdateMeshAccordingInputVariable(destination_variable)
        mesh_utilities_struct.SetReferenceMeshToMesh()

        # propagate non-linear correction from sliding mesh to destination
        self._PerformStructuralSimilarityAnalysis(self._structural_mechanics_analysis_correction, self._structural_model_correction, KSO.OUT_OF_PLANE_DELTA)

        # revert structural model
        mesh_utilities_struct.RevertMeshUpdateAccordingInputVariable(destination_variable)
        mesh_utilities_struct.SetReferenceMeshToMesh()

        # add non-linear correction to shape update
        mesh_utilities = KSO.MeshControllerUtilities(self.destination_model_part)
        mesh_utilities.AddFirstVariableToSecondVariable(KSO.OUT_OF_PLANE_DELTA, destination_variable)

    def InverseMap(self, destination_variable, origin_variable):
        geometry_utilities = KSO.GeometryUtilities(self.destination_model_part)

        geometry_utilities.ProjectNodalVariableOnTangentPlane(
            destination_variable, KSO.BACKGROUND_NORMAL)

        self.vm_mapper.InverseMap(destination_variable, origin_variable)

    def _CorrectOutOfPlanePart(self, destination_variable):
        geometry_utilities = KSO.GeometryUtilities(self.sliding_mesh)
        mesh_utilities = KSO.MeshControllerUtilities(self.sliding_mesh)

        mesh_utilities.UpdateMeshAccordingInputVariable(destination_variable)
        mesh_utilities.SetReferenceMeshToMesh()

        self.spacial_mapper.UpdateInterface()
        self.spacial_mapper.Map(KSO.BACKGROUND_COORDINATE, KSO.BACKGROUND_COORDINATE)
        mesh_utilities.SubtractCoordinatesFromVariable(KSO.BACKGROUND_COORDINATE, KSO.OUT_OF_PLANE_DELTA)
        geometry_utilities.ProjectNodalVariableOnDirection(KSO.OUT_OF_PLANE_DELTA, KSO.BACKGROUND_NORMAL)

        mesh_utilities.RevertMeshUpdateAccordingInputVariable(destination_variable)
        mesh_utilities.SetReferenceMeshToMesh()

    def _PerformStructuralSimilarityAnalysis(self, structural_analysis, structural_model, input_variable):

        mesh_utilities = KSO.MeshControllerUtilities(structural_model.GetModelPart("main"))
        mesh_prescribed_part = structural_model.GetModelPart("main.prescribed_bc")

        # nearest element mapping from destination to sliding mesh
        ne_mapper_settings = KM.Parameters("""{
                                "mapper_type": "nearest_element",
                                "echo_level" : 0
                            }""")
        prescribed_mapper = SpacialMapperFactory.CreateMapper(
            self.destination_model_part, mesh_prescribed_part,
            ne_mapper_settings)

        prescribed_mapper.Map(input_variable, KM.DISPLACEMENT)

        # fix dbc on sliding mesh
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_X, True, mesh_prescribed_part.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_Y, True, mesh_prescribed_part.Nodes)
        KM.VariableUtils().ApplyFixity(KM.DISPLACEMENT_Z, True, mesh_prescribed_part.Nodes)

        structural_analysis.Initialize()
        structural_analysis.RunSolutionLoop()
        structural_analysis.Finalize()

        # nearest element mapping from structural model to destination
        mesh_utilities.RevertMeshUpdateAccordingInputVariable(KM.DISPLACEMENT)
        mesh_utilities.SetReferenceMeshToMesh()

        ne_mapper_settings = KM.Parameters("""{
                                "mapper_type": "nearest_element",
                                "echo_level" : 0
                            }""")
        main_mapper = SpacialMapperFactory.CreateMapper(
            structural_model.GetModelPart("main"), self.destination_model_part,
            ne_mapper_settings)
        main_mapper.Map(KM.DISPLACEMENT, input_variable)

    def _CreateStructuralModel(self, structural_model):

        main_part = structural_model.CreateModelPart("main")
        prescribed_part = main_part.CreateSubModelPart("prescribed_bc")
        fixed_part = main_part.CreateSubModelPart("fixed_bc")

        # add solution variables necessary for structural mechanics analysis
        main_part.AddNodalSolutionStepVariable(KSO.SHAPE_UPDATE)
        main_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        main_part.AddNodalSolutionStepVariable(KM.ROTATION)
        main_part.AddNodalSolutionStepVariable(KM.REACTION)
        main_part.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_PRESSURE)
        main_part.AddNodalSolutionStepVariable(KM.NEGATIVE_FACE_PRESSURE)
        main_part.AddNodalSolutionStepVariable(KSM.POINT_LOAD)
        main_part.AddNodalSolutionStepVariable(KSM.LINE_LOAD)
        main_part.AddNodalSolutionStepVariable(KSM.SURFACE_LOAD)
        main_part.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)
        main_part.AddNodalSolutionStepVariable(KM.ROTATION)
        main_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        main_part.AddNodalSolutionStepVariable(KSM.POINT_MOMENT)
        main_part.AddNodalSolutionStepVariable(KSO.NORMALIZED_SURFACE_NORMAL)

        material_properties = main_part.GetProperties()[0]
        self._ApplyMaterialProperties(material_properties)

        # create nodes of prescribed bc
        prescribed_nodes_added = set()
        for node in self.sliding_mesh.Nodes:
            prescribed_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
            prescribed_nodes_added.add(node.Id)

        # Identify mesh motion area / nodes to be damped
        KM.VariableUtils().SetFlag(KM.STRUCTURE, False, main_part.Nodes)
        radius = self.settings["filter_radius"].GetDouble()
        search_based_functions = KSO.SearchBasedFunctions(self.destination_model_part)
        search_based_functions.FlagNodesInRadius(prescribed_part.Nodes, KM.STRUCTURE, 2*radius)

        # create nodes for main_part
        main_nodes_added = set()
        for node in self.destination_model_part.Nodes:
            if node.Is(KM.STRUCTURE):
                main_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
                main_nodes_added.add(node.Id)

        # create elements for main_part
        ele_id = 0
        node_ids_to_remove = []
        fixed_node_ids = []
        for surface_condition in self.destination_model_part.Conditions:
            node_counter = 0
            number_of_nodes = len(surface_condition.GetNodes())
            ele_node_ids = []
            for node in surface_condition.GetNodes():
                if node.Id in main_nodes_added:
                    ele_node_ids.append(node.Id)
                    node_counter += 1
            if node_counter == number_of_nodes:
                ele_id += 1
                condition_dim = surface_condition.GetGeometry().LocalSpaceDimension()
                if condition_dim != 2:
                    raise Exception("SlidingVertexMorphingMapper: Design model part can only be a surface!")
                if number_of_nodes == 3:
                    main_part.CreateNewElement("ShellThinElement3D3N", ele_id, ele_node_ids, material_properties)
                elif number_of_nodes == 4:
                    main_part.CreateNewElement("ShellThinElement3D4N", ele_id, ele_node_ids, material_properties)
                elif number_of_nodes == 6:
                    edge_nodes = ele_node_ids[:3]
                    node_ids_to_remove.extend(ele_node_ids[3:])
                    main_part.CreateNewElement("ShellThinElement3D3N", ele_id, edge_nodes, material_properties)
                elif number_of_nodes == 8 or number_of_nodes == 9:
                    edge_nodes = ele_node_ids[:4]
                    node_ids_to_remove.extend(ele_node_ids[4:])
                    main_part.CreateNewElement("ShellThinElement3D4N", ele_id, edge_nodes, material_properties)
            else:
                fixed_node_ids.extend(ele_node_ids)

        # remove nodes if higher order surface conditions have been used
        node_ids_to_remove = list(set(node_ids_to_remove))
        nodes_to_remove = main_part.CreateSubModelPart("nodes_to_remove")
        nodes_to_remove.AddNodes(node_ids_to_remove)

        KM.VariableUtils().SetFlag(KM.TO_ERASE, True, nodes_to_remove.Nodes)
        main_part.RemoveNodes(KM.TO_ERASE)

        # add dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, main_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, main_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, main_part)

        # add fixed nodes to fixed_part
        fixed_node_ids = list(set(fixed_node_ids) - set(node_ids_to_remove))
        fixed_part.AddNodes(fixed_node_ids)

    def _ApplyMaterialProperties(self, material_properties):
        # define properties
        # artificially high bending stiffness via large thickness
        radius = self.settings["filter_radius"].GetDouble()
        material_properties.SetValue(KM.YOUNG_MODULUS,1.0E+8)
        material_properties.SetValue(KM.POISSON_RATIO,0.0)
        material_properties.SetValue(KM.THICKNESS,2*radius)
        material_properties.SetValue(KM.DENSITY,0.0)

        constitutive_law = KSM.LinearElasticPlaneStress2DLaw()

        material_properties.SetValue(KM.CONSTITUTIVE_LAW,constitutive_law)
