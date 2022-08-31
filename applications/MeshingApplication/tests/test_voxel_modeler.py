from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as core
import KratosMultiphysics.AdditiveManufacturingApplication as additive
from KratosMultiphysics.vtk_output_process import VtkOutputProcess
import math


def run_modelers(current_model, modelers_list) :

    from KratosMultiphysics.modeler_factory import KratosModelerFactory
    factory = KratosModelerFactory()
    list_of_modelers = factory.ConstructListOfModelers(
        current_model,
        modelers_list)

    for modeler in list_of_modelers:
        modeler.SetupGeometryModel()

    for modeler in list_of_modelers:
        modeler.PrepareGeometryModel()

    for modeler in list_of_modelers:
        modeler.SetupModelPart()


class TestVoxelMeshModeler(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(core.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerNoMdpa(self):

        ## create model
        model = core.Model()

        modelers_list = core.Parameters("""
        [{
            "modeler_name" : "VoxelMeshGeneratorModeler",
            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
            "Parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.01, 0.01, 0.02],
                        "min_point" : [ 0.02, 0.02, -0.02],
                        "max_point" : [ 0.04, 0.04, 0.04]
                    }
                },
                "coloring_settings_list": [
				{
                    "type" : "outer_faces_of_cells_with_color",
					"color": -2,
                    "cell_color": 1
				}
                ],
                "entities_generator_list": [
                {
                    "type" : "elements_with_cell_color",
					"model_part_name": "main_model_part.workpiece",
					"color": 1,
                    "properties_id": 1
				},
                {
                    "type" : "conditions_with_face_color",
                    "model_part_name": "main_model_part.workpiece_boundaries",
                    "color": -2,
                    "properties_id": 2
                }

                ]
            }
        }]""")

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that elements, nodes and conditions were generated correctly
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 12 )
        self.assertEqual( len(main_model_part.Nodes), 36 )
        self.assertEqual( len(main_model_part.Conditions), 32 )

        for element in main_model_part.Elements:
            self.assertGreaterEqual(element.GetGeometry().Volume(), 0.00)

    @KratosUnittest.skipIf(core.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerElementInput(self):

        ## create model
        model = core.Model()

        modelers_list = core.Parameters("""
        [{
            "modeler_name" : "VoxelMeshGeneratorModeler",
            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
            "Parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "data/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]                    }
                },
                "coloring_settings_list": [
				{
                    "type" : "cells_with_inside_center",
					"model_part_name": "skin_model_part.workpiece",
					"color": -1
				}
                ],
                "entities_generator_list": [
                {
                    "type" : "elements_with_cell_color",
					"model_part_name": "main_model_part.workpiece",
					"color": -1,
                    "properties_id": 1
				} 
                ]
            }
        }]""")

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        skin_model_part = model[modelers_list[0]["Parameters"]["input_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(skin_model_part.Elements), 12 )
        self.assertEqual( len(skin_model_part.Nodes), 8 )

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 8 )
        self.assertEqual( len(main_model_part.Nodes), 27 )

        for element in main_model_part.Elements:
            self.assertGreaterEqual(element.GetGeometry().Volume(), 0.00)


    @KratosUnittest.skipIf(core.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerTetrahedralElementInput(self):

        ## create model
        model = core.Model()

        modelers_list = core.Parameters("""
        [{
            "modeler_name" : "VoxelMeshGeneratorModeler",
            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
            "Parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "data/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]                    }
                },
                "coloring_settings_list": [
				{
                    "type" : "cells_with_inside_center",
					"model_part_name": "skin_model_part.workpiece",
					"color": -1
				}
                ],
                "entities_generator_list": [
                {
                    "type" : "tetrahedral_elements_with_cell_color",
					"model_part_name": "main_model_part.workpiece",
					"color": -1,
                    "properties_id": 1
				} 
                ]
            }
        }]""")

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        skin_model_part = model[modelers_list[0]["Parameters"]["input_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(skin_model_part.Elements), 12 )
        self.assertEqual( len(skin_model_part.Nodes), 8 )

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 48 ) # 8 cells with 6 tetrahedra each
        self.assertEqual( len(main_model_part.Nodes), 27 )

        for element in main_model_part.Elements:
            self.assertGreaterEqual(element.GetGeometry().Volume(), 0.00)


    @KratosUnittest.skipIf(core.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerConditionInput(self):

        ## create model
        model = core.Model()

        modelers_list = core.Parameters("""
        [{
            "modeler_name" : "VoxelMeshGeneratorModeler",
            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
            "Parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "data/cube_skin_mesh_conditions",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]                    }
                },
                "coloring_settings_list": [
				{
                    "type" : "cells_with_inside_center",
					"model_part_name": "skin_model_part.workpiece",
					"color": -1,
                    "input_entities": "conditions"
				}
                ],
                "entities_generator_list": [
                {
                    "type" : "elements_with_cell_color",
					"model_part_name": "main_model_part.workpiece",
					"color": -1,
                    "properties_id": 1
				} 
                ]
            }
        }]""")

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 8 )
        self.assertEqual( len(main_model_part.Nodes), 27 )

    @KratosUnittest.skipIf(core.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerConnectedElements(self):

        ## create model
        model = core.Model()

        modelers_list = core.Parameters("""
        [{
            "modeler_name" : "VoxelMeshGeneratorModeler",
            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
            "Parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "data/un_connected_parts",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" :  [0.00, 0.00, 0.00],
                        "max_point" : [ 0.09, 0.09, 0.09]                    
                        }
                },
                "coloring_settings_list": [
				{
                    "type" : "cells_in_touch",
					"model_part_name": "skin_model_part.workpiece",
					"color": -1,
                    "input_entities": "conditions"
				},
				{
                    "type" : "connected_cells_in_touch",
					"model_part_name": "skin_model_part.base",
					"color": -2,
                    "cell_color": -1,
                    "input_entities": "conditions"
				}
                ],
                "entities_generator_list": [
                {
                    "type" : "elements_with_cell_color",
					"model_part_name": "main_model_part.workpiece",
					"color": -2,
                    "properties_id": 1
				}, 
                {
                    "type" : "elements_with_cell_color",
					"model_part_name": "main_model_part.excluded",
					"color": -1,
                    "properties_id": 2
				} 
                ]
            }
        }]""")

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 6 )
        self.assertEqual( len(main_model_part.Nodes), 33 )

        workpiec_model_part = model["main_model_part.workpiece"]
        # Check nb of node and elements
        self.assertEqual( len(workpiec_model_part.Elements), 4 )
        self.assertEqual( len(workpiec_model_part.Nodes), 24 )

        excluded_model_part = model["main_model_part.excluded"]
        # Check nb of node and elements
        self.assertEqual( len(excluded_model_part.Elements), 2 )
        self.assertEqual( len(excluded_model_part.Nodes), 12 ) 


    @KratosUnittest.skipIf(core.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerInRadius(self):

        modelers_list = core.Parameters("""
        [{
            "modeler_name" : "VoxelMeshGeneratorModeler",
            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
            "Parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "data/two_particles",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]                    }
                },
                "coloring_settings_list": [
				{
                    "type" : "cells_with_center_in_sphere_arround_nodes",
					"model_part_name": "skin_model_part",
					"color": -1
				}
                ],
                "entities_generator_list": [
                {
                    "type" : "elements_with_cell_color",
					"model_part_name": "main_model_part.workpiece",
					"color": 1,
                    "properties_id": 1
				} 
                ]
            }
        }]""")

        ## create model
        model = core.Model()

        skin_model_part_name = modelers_list[0]["Parameters"]["input_model_part_name"].GetString()
        skin_model_part = model.CreateModelPart(skin_model_part_name)

        skin_model_part.AddNodalSolutionStepVariable(core.RADIUS)

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]


        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 12 )
        self.assertEqual( len(main_model_part.Nodes), 52 )

        node1 = skin_model_part.GetNode(1)
        node2 = skin_model_part.GetNode(2)
        
        for element in main_model_part.Elements:
            center = element.GetGeometry().Center()
            for i in range(1,3):
                node = skin_model_part.GetNode(i)
                distance_vector = node - center
                distance = math.sqrt(sum([x**2 for x in distance_vector]))
                radius = node.GetSolutionStepValue(core.RADIUS)
                self.assertGreater( distance, radius )

class TestOpenStructureVoxelizer(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(core.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestOpenStructureVoxelizerElementInput(self):

        modelers_list = core.Parameters("""
        [{
            "modeler_name" : "VoxelMeshGeneratorModeler",
            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
            "Parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "data/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]                    }
                },
                "coloring_settings_list": [
				{
                    "type" : "cells_in_touch",
					"model_part_name": "skin_model_part.workpiece",
					"color": -1
				}
                ],
                "entities_generator_list": [
                {
                    "type" : "elements_with_cell_color",
					"model_part_name": "main_model_part.workpiece",
					"color": -1,
                    "properties_id": 1
				} 
                ]
            }
        }]""")

        ## create model
        model = core.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]

        # skin_model_part = model[modelers_list[0]["Parameters"]["input_model_part_name"].GetString()]
        # vtk_output_configuration = core.Parameters("""{
        #         "model_part_name"        : \""""+main_model_part.Name+"""\",
        #         "output_sub_model_parts" : true,
        #         "nodal_solution_step_data_variables" : []
        #     }""")
            
        # vtk_output = VtkOutputProcess(model, vtk_output_configuration)
        # vtk_output.ExecuteInitialize()
        # vtk_output.ExecuteBeforeSolutionLoop()
        # vtk_output.ExecuteInitializeSolutionStep()
        # vtk_output.PrintOutput()
        # vtk_output.ExecuteFinalizeSolutionStep()
        # vtk_output.ExecuteFinalize()

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 18 )
        self.assertEqual( len(main_model_part.Nodes), 48 )

    @KratosUnittest.skipIf(core.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestOpenStructureVoxelizerConditionInput(self):

        modelers_list = core.Parameters("""
        [{
            "modeler_name" : "VoxelMeshGeneratorModeler",
            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
            "Parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "data/cube_skin_mesh_conditions",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]                    }
                },
                "coloring_settings_list": [
				{
                    "type" : "cells_in_touch",
					"model_part_name": "skin_model_part.workpiece",
					"color": -1,
                    "input_entities": "conditions"
				}
                ],
                "entities_generator_list": [
                {
                    "type" : "elements_with_cell_color",
					"model_part_name": "main_model_part.workpiece",
					"color": -1,
                    "properties_id": 1
				} 
                ]
            }
        }]""")

        ## create model
        model = core.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 18 )
        self.assertEqual( len(main_model_part.Nodes), 48 )

    @KratosUnittest.skipIf(core.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestFirstCellTouchVoxelizer(self):

        modelers_list = core.Parameters("""
        [{
            "modeler_name" : "VoxelMeshGeneratorModeler",
            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
            "Parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "data/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.05,-0.02, 0.04],
                        "max_point" : [ 0.01, 0.01, 0.07]                    }
                },
                "coloring_settings_list": [
				{
                    "type" : "cells_in_touch",
					"model_part_name": "skin_model_part.workpiece",
					"color": -1
				}
                ],
                "entities_generator_list": [
                {
                    "type" : "elements_with_cell_color",
					"model_part_name": "main_model_part.workpiece",
					"color": -1,
                    "properties_id": 1
				} 
                ]
            }
        }]""")

        ## create model
        model = core.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]

        # skin_model_part = model[modelers_list[0]["Parameters"]["input_model_part_name"].GetString()]
        # vtk_output_configuration = core.Parameters("""{
        #         "model_part_name"        : \""""+main_model_part.Name+"""\",
        #         "output_sub_model_parts" : true,
        #         "nodal_solution_step_data_variables" : []
        #     }""")
            
        # vtk_output = VtkOutputProcess(model, vtk_output_configuration)
        # vtk_output.ExecuteInitialize()
        # vtk_output.ExecuteBeforeSolutionLoop()
        # vtk_output.ExecuteInitializeSolutionStep()
        # vtk_output.PrintOutput()
        # vtk_output.ExecuteFinalizeSolutionStep()
        # vtk_output.ExecuteFinalize()

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 1 )
        self.assertEqual( len(main_model_part.Nodes), 8 )

class TestBoundaryConditionVoxelizer(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(core.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestWallConditionVoxelizer(self):

        modelers_list = core.Parameters("""
        [{
            "modeler_name" : "VoxelMeshGeneratorModeler",
            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
            "Parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "data/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]                    }
                },
                "coloring_settings_list": [
				{
                    "type" : "cells_with_inside_center",
					"model_part_name": "skin_model_part.workpiece",
					"color": -1
				},
				{
                    "type" : "cells_faces",
					"model_part_name": "skin_model_part.workpiece",
					"color": -1,
                    "cell_color": -1
				}
                ],
                "entities_generator_list": [
                {
                    "type" : "conditions_with_face_color",
					"model_part_name": "main_model_part.workpiece",
					"color": -1,
                    "properties_id": 1
				} 
                ]
            }
        }]""")

        ## create model
        model = core.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Conditions), 24 )
        self.assertEqual( len(main_model_part.Nodes), 26 )

        for condition in main_model_part.Conditions:
            self.assertGreaterEqual(condition.GetGeometry().Area(), 0.00)

    def test_TestOuterWallConditionVoxelizer(self):

        modelers_list = core.Parameters("""
        [{
            "modeler_name" : "VoxelMeshGeneratorModeler",
            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
            "Parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "data/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]                    }
                },
                "coloring_settings_list": [
				{
                    "type" : "cells_with_inside_center",
					"model_part_name": "skin_model_part.workpiece",
					"color": -1
				},
				{
                    "type" : "outer_faces_of_cells_with_color",
					"color": -1,
                    "cell_color": -1
				}
                ],
                "entities_generator_list": [
                {
                    "type" : "conditions_with_face_color",
					"model_part_name": "main_model_part.workpiece",
					"color": -1,
                    "properties_id": 1
				} 
                ]
            }
        }]""")

        ## create model
        model = core.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Conditions), 24 )
        self.assertEqual( len(main_model_part.Nodes), 26 )

        for condition in main_model_part.Conditions:
            self.assertGreaterEqual(condition.GetGeometry().Area(), 0.00)

    @KratosUnittest.skip("Missing mdpa")
    def test_TestDiagonalCubeWallConditionVoxelizer(self):

        modelers_list = core.Parameters("""
        [{
            "modeler_name" : "VoxelMeshGeneratorModeler",
            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
            "Parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "data/diagonal_cube_skin_mesh",
                "key_plane_generator": {
                            "Parameters": {
                                "max_point": [
                                    0.12,
                                    0.12,
                                    0.12
                                ],
                                "min_point": [
                                    -0.04,
                                    -0.04,
                                    0.02
                                ],
                                "voxel_sizes": [
                                    0.004,
                                    0.004,
                                    0.004
                                ]
                            }
                },
                "coloring_settings_list": [
				{
                    "type" : "cells_with_inside_center",
					"model_part_name": "workpiece",
					"color": -1
				},
				{
                    "type" : "cells_faces",
					"model_part_name": "workpiece",
					"color": -1,
                    "cell_color": -1
				}
                ],
                "entities_generator_list": [
                {
                    "type" : "conditions_with_face_color",
					"model_part_name": "workpiece",
					"color": -1,
                    "properties_id": 1
				} 
                ]
            }
        }]""")

        ## create model
        model = core.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        # self.assertEqual( len(main_model_part.Conditions), 24 )
        # self.assertEqual( len(main_model_part.Nodes), 27 )

        for condition in main_model_part.Conditions:
            self.assertGreaterEqual(condition.GetGeometry().Area(), 0.00)


    @KratosUnittest.skip("Missing mdpa")
    def test_TestCompleteModel(self):

        modelers_list = core.Parameters("""
            [
				{
		            "modeler_name" : "VoxelMeshGeneratorModeler",
		            "kratos_module" : "KratosMultiphysics.AdditiveManufacturingApplication",
		            "Parameters" : {
		                "output_model_part_name" : "main_model_part",
		                "input_model_part_name" : "input_model_part",
		                "mdpa_file_name" : "data/complete_model",
		                "key_plane_generator": {
		                    "Parameters" : {
								"voxel_sizes" : [0.005, 0.005, 0.005],
		                        "min_point" : [ -0.02,-0.03, 0.0],
		                        "max_point" : [ 0.078656,  0.02, 0.113746]
								}
		                },
		                "coloring_settings_list": [
							{
			                    "type" : "cells_in_touch",
								"model_part_name": "input_model_part.support.support_body",
								"color": -1,
						        "input_entities": "elements"

							},
							{
			                    "type" : "cells_with_inside_center",
								"model_part_name": "input_model_part.printpart.printpart_lateral_bc",
								"color": -2,
								"input_entities": "conditions"
							},
							{
			                    "type" : "cells_faces",
								"model_part_name": "input_model_part.support.support_base_bc",
								"color": -3,
			                    "cell_color": -1,
								"input_entities": "conditions"
							},
							{
			                    "type" : "cells_faces",
								"model_part_name": "input_model_part.printpart.printpart_base_bc",
								"color": -4,
			                    "cell_color": -2,
								"input_entities": "conditions"
							}						
		                ],
		                "entities_generator_list": [
							{
			                    "type" : "elements_with_cell_color",
								"model_part_name": "main_model_part.support_output.supports_volume_output",
								"color": -1,
			                    "properties_id": 1
							},
			                {
			                    "type" : "elements_with_cell_color",
								"model_part_name": "main_model_part.printpart_output.printpart_volume_output",
								"color": -2,
			                    "properties_id": 2
							},
							{
			                    "type" : "conditions_with_face_color",
								"model_part_name": "main_model_part.support_output.support_base_bc_output",
								"color": -3,
			                    "properties_id": 3
							},
							{
			                    "type" : "conditions_with_face_color",
								"model_part_name": "main_model_part.printpart_output.printpart_base_bc_output",
								"color": -4,
			                    "properties_id": 4
							}
		                ]
		            }
		        }
			 ]
        """)

        ## create model
        model = core.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["Parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        # self.assertEqual( len(main_model_part.Conditions), 24 )
        # self.assertEqual( len(main_model_part.Nodes), 27 )

        for condition in main_model_part.Conditions:
            self.assertGreaterEqual(condition.GetGeometry().Area(), 0.00)



if __name__ == '__main__':
    
    tester = TestVoxelMeshModeler()
    tester.test_TestVoxelMeshModelerNoMdpa()
    tester.test_TestVoxelMeshModelerElementInput()
    tester.test_TestVoxelMeshModelerConditionInput()
    tester.test_TestVoxelMeshModelerTetrahedralElementInput()
    tester.test_TestVoxelMeshModelerConnectedElements()
    tester.test_TestVoxelMeshModelerInRadius()

    tester = TestOpenStructureVoxelizer()
    tester.test_TestOpenStructureVoxelizerElementInput()
    tester.test_TestOpenStructureVoxelizerConditionInput()
    tester.test_TestFirstCellTouchVoxelizer()

    tester = TestBoundaryConditionVoxelizer()
    tester.test_TestWallConditionVoxelizer()
    tester.test_TestOuterWallConditionVoxelizer()
