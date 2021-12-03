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
import KratosMultiphysics.StructuralMechanicsApplication as KSM
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
                "echo_level"      : 0,
                "analysis_type"   : "linear",
                "model_part_name" : "Main",
                "domain_size"     : 3,
                "time_stepping" : {
                    "time_step"       : 1.0
                },
                "model_import_settings"              : {
                    "input_type"     : "use_input_model_part"
                },
                "line_search"                        : false,
                "convergence_criterion"              : "residual_criterion",
                "displacement_relative_tolerance"    : 0.0001,
                "displacement_absolute_tolerance"    : 1e-9,
                "residual_relative_tolerance"        : 0.0001,
                "residual_absolute_tolerance"        : 1e-9,
                "max_iteration"                      : 10,
                "linear_solver_settings"             : {
                    "solver_type" : "LinearSolversApplication.sparse_lu",
                    "scaling"     : false,
                    "verbosity"   : 0
                },
                "rotation_dofs"                      : true
            },
            "processes" : {
                "constraints_process_list" : [
                {
                    "python_module": "mpc_normal_support_process",
                    "process_name": "MPCNormalSupportProcess",
                    "Parameters": {
                        "model_part_name": "Main.without_edges"
                    }
                },
                {
                    "python_module" : "assign_vector_variable_process",
                    "kratos_module" : "KratosMultiphysics",
                    "help"          : "This process fixes the selected components of a given vector variable",
                    "process_name"  : "AssignVectorVariableProcess",
                    "Parameters"    : {
                        "mesh_id"         : 0,
                        "model_part_name" : "Main.FixedConditions",
                        "variable_name"   : "DISPLACEMENT",
                        "value"           : [0.0,0.0,0.0]
                    }
                }
                ],
                "loads_process_list"       : []
            }
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

        self._output_boundary_geometry()
        
        #self._structural_mechanics_analysis = StructuralMechanicsAnalysis(model, self.MeshSolverSettings)

    def Update(self):
        self.vm_mapper.Update()
        self.spacial_mapper.UpdateInterface()
        self.spacial_mapper.Map(KSO.NORMALIZED_SURFACE_NORMAL, KSO.BACKGROUND_NORMAL)

        # Model kopieren als Structural model
        self._structural_model = KM.Model()

        structural_mechanics_analysis_setings = self.GetStructuralSimilaritySetttings()
        self._structural_mechanics_analysis = StructuralMechanicsAnalysis(self._structural_model, structural_mechanics_analysis_setings)
        self._CreateStructuralModel()

        # print(self._structural_model)

    def Map(self, origin_variable, destination_variable):
        self.vm_mapper.Map(origin_variable, destination_variable)

        geometry_utilities = KSO.GeometryUtilities(self.sliding_edge_mesh)
        geometry_utilities.ProjectNodalVariableOnTangentPlane(
            destination_variable, KSO.BACKGROUND_NORMAL)

        # non-linear correction
        for node in self.sliding_edge_mesh.Nodes:
            if node.Id == 1:
                print("node {} coordinates bevor _CorrectOutOfPlanePart [{},{}.{}]".format(node.Id,node.X, node.Y, node.Z))

        self._CorrectOutOfPlanePart(destination_variable)

        # self._structural_model = None
        # self._structural_mechanics_analysis = None

    def InverseMap(self, destination_variable, origin_variable):
        #geometry_utilities = KSO.GeometryUtilities(self.sliding_edge_mesh)
        geometry_utilities = KSO.GeometryUtilities(self.destination_model_part)

        for node in self.sliding_edge_mesh.Nodes:
            if node.Id == 1:
                print("node {} coordinates bevor _PerformStructuralSimilarityAnalysis [{},{}.{}]".format(node.Id,node.X, node.Y, node.Z))


        self._PerformStructuralSimilarityAnalysis(KSO.BACKGROUND_NORMAL)

        geometry_utilities.ProjectNodalVariableOnTangentPlane(
            destination_variable, KSO.BACKGROUND_NORMAL)

        self.vm_mapper.InverseMap(destination_variable, origin_variable)

    def _CorrectOutOfPlanePart(self, destination_variable):
        geometry_utilities = KSO.GeometryUtilities(self.sliding_edge_mesh)
        mesh_utilities = KSO.MeshControllerUtilities(self.sliding_edge_mesh)

        for node in self.destination_model_part.Nodes:
            node.SetValue(KM.NODE_PROPERTY_ID, node.Id)

        for node in self.sliding_edge_mesh.Nodes:
            if node.Id == 1:
                print("node {} coordinates bevor [{},{}.{}]".format(node.Id,node.X, node.Y, node.Z))
        mesh_utilities.UpdateMeshAccordingInputVariable(destination_variable)
        for node in self.sliding_edge_mesh.Nodes:
            if node.Id == 1:
                print("node {} coordinates danach [{},{}.{}]".format(node.Id,node.X, node.Y, node.Z))

        mesh_utilities.SetReferenceMeshToMesh()
        # for node in self.sliding_edge_mesh.Nodes:
        #     print("node {} coordinates nach SetReference [{},{}.{}]".format(node.Id,node.X, node.Y, node.Z))

        self.spacial_mapper.UpdateInterface()
        self.spacial_mapper.Map(KSO.BACKGROUND_COORDINATE, KSO.BACKGROUND_COORDINATE)
        mesh_utilities.SubtractCoordinatesFromVariable(KSO.BACKGROUND_COORDINATE, KSO.OUT_OF_PLANE_DELTA)
        geometry_utilities.ProjectNodalVariableOnDirection(KSO.OUT_OF_PLANE_DELTA, KSO.BACKGROUND_NORMAL)

        mesh_utilities.RevertMeshUpdateAccordingInputVariable(destination_variable)
        mesh_utilities.SetReferenceMeshToMesh()

        mesh_utilities.AddFirstVariableToSecondVariable(KSO.OUT_OF_PLANE_DELTA, destination_variable)

    def _PerformStructuralSimilarityAnalysis(self, input_variable):
        
        mesh_utilities = KSO.MeshControllerUtilities(self._structural_model.GetModelPart("Main"))

        for node in self._structural_model.GetModelPart("Main.MeshMotionConditions").Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT, node.GetSolutionStepValue(input_variable))
            node.Fix(KM.DISPLACEMENT_X)
            node.Fix(KM.DISPLACEMENT_Y)
            node.Fix(KM.DISPLACEMENT_Z)

        #self._structural_mechanics_analysis.Run()
        self._structural_mechanics_analysis.Initialize()
        self._structural_mechanics_analysis.RunSolutionLoop()
        #self._structural_mechanics_analysis.Finalize()


        # for node in self._structural_model.GetModelPart("Main.without_edges").Nodes:
        #     node.SetSolutionStepValue(input_variable, node.GetSolutionStepValue(KM.DISPLACEMENT))   

        for node in self.destination_model_part.Nodes:
            node.SetSolutionStepValue(input_variable, node.GetSolutionStepValue(KM.DISPLACEMENT))

        # for node in self._structural_model.GetModelPart("Main.MeshMotionConditions").Nodes:
        #     node.SetSolutionStepValue(KM.DISPLACEMENT, [0.0,0.0,0.0])


        self._structural_mechanics_analysis.Clear()
        self._structural_mechanics_analysis.Finalize()

        for node in self._structural_model.GetModelPart("Main.MeshMotionConditions").Nodes:
            node.Free(KM.DISPLACEMENT_X)
            node.Free(KM.DISPLACEMENT_Y)
            node.Free(KM.DISPLACEMENT_Z)

        mesh_utilities.RevertMeshUpdateAccordingInputVariable(KM.DISPLACEMENT)
        mesh_utilities.SetReferenceMeshToMesh()
        

        #for node in 

        #self._output_structural_model()


    def _CreateStructuralModel(self):
        
        
        #main_part = self._structural_model.CreateModelPart("Main")
        main_part = self._structural_model.GetModelPart("Main")
        #main_part.AddNodalSolutionStepVariable(KSO.BACKGROUND_NORMAL)
        mesh_motion_conditions = main_part.CreateSubModelPart("MeshMotionConditions")
        fixed_nodes = main_part.CreateSubModelPart("FixedConditions")
        #main_part.AddNodalSolutionStepVariable(KM.STRUCTURE)

        material_properties = main_part.GetProperties()[0]
        self._apply_material_properties(material_properties)
        #main_part.AddProperties(material_properties)
        #print(material_properties)
        for node in self.sliding_edge_mesh.Nodes:
            # neue Knoten erstellen! können selbe IDs haben
            #mesh_motion_conditions.CreateNewNode(node.Id, node.X, node.Y, node.Z)
            mesh_motion_conditions.AddNode(node, 0)

            #print(node.GetValue(KSO.BACKGROUND_NORMAL))
        
        # Identify mesh motion area / nodes to be damped
        KM.VariableUtils().SetFlag(KM.STRUCTURE, False, main_part.Nodes)
        radius = self.settings["filter_radius"].GetDouble()
        search_based_functions = KSO.SearchBasedFunctions(self.destination_model_part)
        search_based_functions.FlagNodesInRadius(mesh_motion_conditions.Nodes, KM.STRUCTURE, 2*radius)

        for node in self.destination_model_part.Nodes:
            if node.Is(KM.STRUCTURE):
                #main_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
                main_part.AddNode(node, 0)
                node.SetValue(KM.NODE_PROPERTY_ID, node.Id)

        ele_id = 0
        for surface_condition in self.destination_model_part.Conditions:
            node_counter = 0
            number_of_nodes = len(surface_condition.GetNodes())
            ele_node_ids = []
            for node in surface_condition.GetNodes():
                if main_part.HasNode(node.Id):
                    ele_node_ids.append(node.Id)
                    node_counter += 1
                else:
                    break

            if node_counter == number_of_nodes:
                    ele_id += 1
                    # check LocalDimension of condition = 2? wenn 1 oder 3 => Fehler
                    # wieviele knoten hat condition? => erstelle Element 3D3N, 3D4N; erstelle für 3D6N, 3D8N jeweils 2 element
                    # dafür brauche ich die eckknoten (indices sind immer gleich => geometrie aus z.b. triangle_3d_6.h), lösche die mittelkantenknoten
                    main_part.CreateNewElement("ShellThinElement3D3N", ele_id, ele_node_ids, material_properties)
                    main_part.AddCondition(surface_condition)

        
        # find edges
        # main_part.AddNodalSolutionStepVariable(NEIGHBOUR_NODES)
        # main_part.AddNodalSolutionStepVariable(NEIGHBOUR_ELEMENTS)
        # neighbor_search = KM.FindNodalNeighboursProcess(main_part)
        # neighbor_search.ClearNeighbours()
        # neighbor_search.Execute()

        neighbor_node_search = KM.FindGlobalNodalNeighboursProcess(main_part)
        neighbor_node_search.ClearNeighbours()
        neighbor_node_search.Execute()
        neighbour_element_search = KM.FindGlobalNodalElementalNeighboursProcess(main_part)
        neighbour_element_search.ClearNeighbours()
        neighbour_element_search.Execute()

        KSO.GeometryUtilities(main_part).ExtractEdgeNodes("FixedConditions")
        #edge_part = main_part.GetSubModelPart("auto_surface_nodes")
       
        node_ids_to_remove = []
        for node in fixed_nodes.Nodes:
            if mesh_motion_conditions.HasNode(node.Id):
                node_ids_to_remove.append(node.Id)

        for node_id in node_ids_to_remove:
            fixed_nodes.RemoveNode(node_id)

        #self._output_edge_model(main_part)
        # for node in main_part.Nodes:
        #     neighbour_node_ids = node.GetValue(KM.NEIGHBOUR_NODES)
        #     for neighbour_node_id in neighbour_node_ids:
        #         neighbour_node = main_part.GetNode(neighbour_node_id)
        #         neighbour_element_ids = node.GetValue(KM.NEIGHBOUR_ELEMENTS) 

        main_part.Check(main_part.ProcessInfo)
        

        without_edges = main_part.CreateSubModelPart("without_edges")
        for node in main_part.Nodes:
            without_edges.AddNode(node, 0)

        node_ids_to_remove = []
        for node in fixed_nodes.Nodes:
            node_ids_to_remove.append(node.Id)
        for node in mesh_motion_conditions.Nodes:
            node_ids_to_remove.append(node.Id)
        for node_id in node_ids_to_remove:
            without_edges.RemoveNode(node_id)

        ele_id = 0
        for surface_condition in main_part.Conditions:
            node_counter = 0
            number_of_nodes = len(surface_condition.GetNodes())
            ele_node_ids = []
            for node in surface_condition.GetNodes():
                if without_edges.HasNode(node.Id):
                    ele_node_ids.append(node.Id)
                    node_counter += 1
                else:
                    break
            if node_counter == number_of_nodes:
                    ele_id += 1
                    without_edges.AddCondition(surface_condition)

        # neighbor_node_search = KM.FindGlobalNodalNeighboursProcess(main_part)
        # neighbor_node_search.ClearNeighbours()
        # neighbor_node_search.Execute()
        # neighbour_node_id_map = neighbor_node_search.GetNeighbourIds(main_part.Nodes)
        # print(neighbour_node_id_map)


        #self._output_structural_model()

        # neighbour_element_search = KM.FindGlobalNodalElementalNeighboursProcess(main_part)
        # neighbour_element_search.ClearNeighbours()
        # neighbour_element_search.Execute()

        # for node in main_part.Nodes:
        #     print(node.GetValue(KM.NUMBER_OF_NEIGHBOUR_ELEMENTS))
        # neighbour_element_id_map = neighbour_element_search.GetNeighbourIds(main_part.Nodes)        
        # print(neighbour_element_id_map)

        # for node in main_part.Nodes:
        #     for neighbour_node_id in neighbour_node_id_map[node.Id]:
        #         neighbour_node = main_part.GetNode(neighbour_node_id)
        #         neighbour_element_ids = node.GetValue(KM.NEIGHBOUR_ELEMENTS) 


    # def _apply_material_properties(self, model):
    #     #define properties
    #     model.GetProperties()[1].SetValue(KM.YOUNG_MODULUS,2.06900E+8)
    #     model.GetProperties()[1].SetValue(KM.POISSON_RATIO,2.90000E-01)
    #     model.GetProperties()[1].SetValue(KM.THICKNESS,0.1)
    #     model.GetProperties()[1].SetValue(KM.DENSITY,7.85)

    #     constitutive_law = KSM.LinearElasticPlaneStress2DLaw()

    #     model.GetProperties()[1].SetValue(KM.CONSTITUTIVE_LAW,constitutive_law)

    def _output_structural_model(self):
        from KratosMultiphysics.vtk_output_process import VtkOutputProcess

        self.vtk_output_process = VtkOutputProcess(self._structural_model,
                                    KM.Parameters("""{
                                            "model_part_name"                    : "Main",
                                            "nodal_solution_step_data_variables" : ["DISPLACEMENT"],
                                            "nodal_data_value_variables"         : ["NODE_PROPERTY_ID", "BACKGROUND_NORMAL"],
                                            "nodal_flags"                        : ["STRUCTURE"],
                                            "write_ids"                          : true
                                        }
                                        """)
                                    )

        self.vtk_output_process.ExecuteInitialize()
        self.vtk_output_process.ExecuteBeforeSolutionLoop()
        self.vtk_output_process.ExecuteInitializeSolutionStep()
        self.vtk_output_process.PrintOutput()
        self.vtk_output_process.ExecuteFinalizeSolutionStep()
        self.vtk_output_process.ExecuteFinalize()

    def _output_edge_model(self, Main):
        from KratosMultiphysics.vtk_output_process import VtkOutputProcess

        self.vtk_output_process = VtkOutputProcess(self._structural_model,
                                    KM.Parameters("""{
                                            "model_part_name"                    : "Main.FixedConditions",
                                            "nodal_solution_step_data_variables" : [],
                                            "nodal_data_value_variables"         : ["NODE_PROPERTY_ID"],
                                            "nodal_flags"                        : ["STRUCTURE"],
                                            "write_ids"                          : true
                                        }
                                        """)
                                    )

        self.vtk_output_process.ExecuteInitialize()
        self.vtk_output_process.ExecuteBeforeSolutionLoop()
        self.vtk_output_process.ExecuteInitializeSolutionStep()
        self.vtk_output_process.PrintOutput()
        self.vtk_output_process.ExecuteFinalizeSolutionStep()
        self.vtk_output_process.ExecuteFinalize()


    def _apply_material_properties(self, material_properties):
        #define properties
        material_properties.SetValue(KM.YOUNG_MODULUS,2.06900E+8)
        material_properties.SetValue(KM.POISSON_RATIO,2.90000E-01)
        material_properties.SetValue(KM.THICKNESS,0.1)
        material_properties.SetValue(KM.DENSITY,7.85)

        constitutive_law = KSM.LinearElasticPlaneStress2DLaw()

        material_properties.SetValue(KM.CONSTITUTIVE_LAW,constitutive_law)


    def _GetMaterialProperties(self):
        return KM.Parameters("""{
            "model_part_name": "Main",
            "properties_id": 0,
            "Material": {
                    "name": "Shell_Material",
                    "constitutive_law": {
                            "name": "LinearElasticPlaneStress2DLaw"
                    },
                    "Variables": {
                            "DENSITY": 7.85,
                            "YOUNG_MODULUS": 2.06900E+8,
                            "POISSON_RATIO": 2.90000E-01,
                            "THICKNESS": 1.00000E-01
                    },
                    "Tables": {}
            }
        }""")

    def _output_boundary_geometry(self):
        from KratosMultiphysics.vtk_output_process import VtkOutputProcess

        self.vtk_output_process = VtkOutputProcess(self._background_model,
                                    KM.Parameters("""{
                                            "model_part_name"                    : "background_mesh",
                                            "nodal_solution_step_data_variables" : [],
                                            "nodal_data_value_variables"         : ["NODE_PROPERTY_ID"],
                                            "nodal_flags"                        : [],
                                            "write_ids"                          : true
                                        }
                                        """)
                                    )

        self.vtk_output_process.ExecuteInitialize()
        self.vtk_output_process.ExecuteBeforeSolutionLoop()
        self.vtk_output_process.ExecuteInitializeSolutionStep()
        self.vtk_output_process.PrintOutput()
        self.vtk_output_process.ExecuteFinalizeSolutionStep()
        self.vtk_output_process.ExecuteFinalize()

    def _output_destination_geometry(self):
        from KratosMultiphysics.vtk_output_process import VtkOutputProcess

        self.vtk_output_process = VtkOutputProcess(self.destination_model_part,
                                    KM.Parameters("""{
                                            "model_part_name"                    : "background_mesh",
                                            "nodal_solution_step_data_variables" : [],
                                            "nodal_data_value_variables"         : ["NODE_PROPERTY_ID"],
                                            "nodal_flags"                        : [],
                                            "write_ids"                          : true
                                        }
                                        """)
                                    )

        self.vtk_output_process.ExecuteInitialize()
        self.vtk_output_process.ExecuteBeforeSolutionLoop()
        self.vtk_output_process.ExecuteInitializeSolutionStep()
        self.vtk_output_process.PrintOutput()
        self.vtk_output_process.ExecuteFinalizeSolutionStep()
        self.vtk_output_process.ExecuteFinalize()

