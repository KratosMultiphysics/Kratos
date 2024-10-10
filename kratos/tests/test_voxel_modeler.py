import math

import KratosMultiphysics as kratos
from KratosMultiphysics import KratosUnittest
from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory


def run_modelers(model, settings):
    factory = KratosModelParametersFactory(model)
    list_of_modelers = factory.ConstructListOfItems(settings)

    for modeler in list_of_modelers:
        modeler.SetupGeometryModel()

    for modeler in list_of_modelers:
        modeler.PrepareGeometryModel()

    for modeler in list_of_modelers:
        modeler.SetupModelPart()


class TestVoxelMeshModeler(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerNoMdpa(self):

        ## create model
        model = kratos.Model()

        modelers_list = kratos.Parameters("""[{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
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
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 12 )
        self.assertEqual( len(main_model_part.Nodes), 36 )
        self.assertEqual( len(main_model_part.Conditions), 32 )

        for element in main_model_part.Elements:
            self.assertGreaterEqual(element.GetGeometry().Volume(), 0.00)


    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerElementInput(self):

        ## create model
        model = kratos.Model()

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]
                    }
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
        skin_model_part = model[modelers_list[0]["parameters"]["input_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(skin_model_part.Elements), 12 )
        self.assertEqual( len(skin_model_part.Nodes), 8 )

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 8 )
        self.assertEqual( len(main_model_part.Nodes), 27 )

        for element in main_model_part.Elements:
            self.assertGreaterEqual(element.GetGeometry().Volume(), 0.00)


    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerQuadraticElementInput(self):

        ## create model
        model = kratos.Model()

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]
                    }
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
                    "type" : "quadratic_elements_with_cell_color",
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
        skin_model_part = model[modelers_list[0]["parameters"]["input_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(skin_model_part.Elements), 12 )
        self.assertEqual( len(skin_model_part.Nodes), 8 )

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 8 )
        self.assertEqual( len(main_model_part.Nodes), 125 )

        expected_volume = 0.03**3
        for element in main_model_part.Elements:
            self.assertEqualTolerance(element.GetGeometry().Volume(), expected_volume, 1e-12)


    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerElementInputAutoSkin(self):

        ## create model
        model = kratos.Model()

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/cube_skin_mesh",
                "key_plane_generator": {
                    "type" : "outer_shell",
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03]
                    }
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
        skin_model_part = model[modelers_list[0]["parameters"]["input_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(skin_model_part.Elements), 12 )
        self.assertEqual( len(skin_model_part.Nodes), 8 )

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 8 )
        self.assertEqual( len(main_model_part.Nodes), 27 )

        for element in main_model_part.Elements:
            self.assertGreaterEqual(element.GetGeometry().Volume(), 0.00)


    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerTetrahedralElementInput(self):

        ## create model
        model = kratos.Model()

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]
                    }
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
        skin_model_part = model[modelers_list[0]["parameters"]["input_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(skin_model_part.Elements), 12 )
        self.assertEqual( len(skin_model_part.Nodes), 8 )

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 48 ) # 8 cells with 6 tetrahedra each
        self.assertEqual( len(main_model_part.Nodes), 27 )

        for element in main_model_part.Elements:
            self.assertGreaterEqual(element.GetGeometry().Volume(), 0.00)


    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerConditionInput(self):

        ## create model
        model = kratos.Model()

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/cube_skin_mesh_conditions",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]
                    }
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
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 8 )
        self.assertEqual( len(main_model_part.Nodes), 27 )


    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerConnectedElements(self):

        ## create model
        model = kratos.Model()

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/un_connected_parts",
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
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

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


    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestVoxelMeshModelerInRadius(self):

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/two_particles",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]
                    }
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
        model = kratos.Model()

        skin_model_part_name = modelers_list[0]["parameters"]["input_model_part_name"].GetString()
        skin_model_part = model.CreateModelPart(skin_model_part_name)

        skin_model_part.AddNodalSolutionStepVariable(kratos.RADIUS)

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]


        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 12 )
        self.assertEqual( len(main_model_part.Nodes), 52 )

        for element in main_model_part.Elements:
            center = element.GetGeometry().Center()
            for i in range(1,3):
                node = skin_model_part.GetNode(i)
                distance_vector = node - center
                distance = math.sqrt(sum([x**2 for x in distance_vector]))
                radius = node.GetSolutionStepValue(kratos.RADIUS)
                self.assertGreater( distance, radius )


class TestOpenStructureVoxelizer(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestOpenStructureVoxelizerElementInput(self):

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/cube_skin_mesh",
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
        model = kratos.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 18 )
        self.assertEqual( len(main_model_part.Nodes), 48 )


    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestOpenStructureVoxelizerConditionInput(self):

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/cube_skin_mesh_conditions",
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
        model = kratos.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 18 )
        self.assertEqual( len(main_model_part.Nodes), 48 )


    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestFirstCellTouchVoxelizer(self):

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.05,-0.02, 0.04],
                        "max_point" : [ 0.01, 0.01, 0.07]
                    }
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
        model = kratos.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Elements), 1 )
        self.assertEqual( len(main_model_part.Nodes), 8 )


class TestBoundaryConditionVoxelizer(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestWallConditionVoxelizer(self):

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]
                    }
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
        model = kratos.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Conditions), 24 )
        self.assertEqual( len(main_model_part.Nodes), 26 )

        for condition in main_model_part.Conditions:
            self.assertGreaterEqual(condition.GetGeometry().Area(), 0.00)


    @KratosUnittest.skipIf(kratos.IsDistributedRun(), "This test is designed for serial runs only.")
    def test_TestFacesBetweenColors(self):

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]
                    }
                },
                "coloring_settings_list": [
                {
                    "type" : "cells_with_inside_center",
                    "model_part_name": "skin_model_part.workpiece",
                    "color": -1
                },
                {
                    "type" : "cells_faces_between_colors",
                    "outside_color" : 1.0,
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
        model = kratos.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Conditions), 8 )
        self.assertEqual( len(main_model_part.Nodes), 15 )

        for condition in main_model_part.Conditions:
            self.assertGreaterEqual(condition.GetGeometry().Area(), 0.00)


    def test_TestOuterWallConditionVoxelizer(self):

        modelers_list = kratos.Parameters("""
        [{
            "name" : "KratosMultiphysics.VoxelMeshGeneratorModeler",
            "parameters" : {
                "output_model_part_name" : "main_model_part",
                "input_model_part_name" : "skin_model_part",
                "mdpa_file_name" : "test_files/voxel_mesh_modeler/cube_skin_mesh",
                "key_plane_generator": {
                    "Parameters" : {
                        "voxel_sizes" : [0.03, 0.03, 0.03],
                        "min_point" : [-0.03,-0.03, 0.03],
                        "max_point" : [ 0.09, 0.09, 0.09]
                    }
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
        model = kratos.Model()

        # create modeler obj
        run_modelers(model, modelers_list)

        # Check that all element have deactivated in ExecuteBeforeSolutionLoop()
        main_model_part = model[modelers_list[0]["parameters"]["output_model_part_name"].GetString()]

        # Check nb of node and elements
        self.assertEqual( len(main_model_part.Conditions), 24 )
        self.assertEqual( len(main_model_part.Nodes), 26 )

        for condition in main_model_part.Conditions:
            self.assertGreaterEqual(condition.GetGeometry().Area(), 0.00)


if __name__ == '__main__':
    KratosUnittest.main()
