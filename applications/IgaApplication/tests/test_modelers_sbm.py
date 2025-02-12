import KratosMultiphysics
import KratosMultiphysics.IgaApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

def run_modelers(current_model, modelers_list):
    from KratosMultiphysics.modeler_factory import KratosModelerFactory
    factory = KratosModelerFactory()
    list_of_modelers = factory.ConstructListOfModelers(current_model, modelers_list)

    for modeler in list_of_modelers:
        modeler.SetupGeometryModel()

    for modeler in list_of_modelers:
        modeler.PrepareGeometryModel()

    for modeler in list_of_modelers:
        modeler.SetupModelPart()

class TestModelersSbm(KratosUnittest.TestCase):
    
    # test the call to the nurbs_geometry_modeler_sbm to create a rectangle + the breps
    def test_nurbs_geometry_2d_modeler_control_points(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsGeometryModelerSbm",
            "Parameters": {
                "model_part_name" : "IgaModelPart",
                "lower_point_xyz": [0.0,0.0,0.0],
                "upper_point_xyz": [1.0,1.0,0.0],
                "lower_point_uvw": [0.0,0.0,0.0],
                "upper_point_uvw": [1.0,1.0,0.0],
                "polynomial_order" : [4, 1],
                "number_of_knot_spans" : [3,2]
            }
        }]
        """)
        run_modelers(current_model, modeler_settings)

        model_part = current_model.GetModelPart("IgaModelPart")
        # Check control points
        nodes_mp = model_part.NodesArray(0)

        for i in range(3):
            self.assertAlmostEqual(nodes_mp[i*7+0].X,0.0)
            self.assertAlmostEqual(nodes_mp[i*7+1].X,0.083333333334)
            self.assertAlmostEqual(nodes_mp[i*7+2].X,0.25)
            self.assertAlmostEqual(nodes_mp[i*7+3].X,0.5)
            self.assertAlmostEqual(nodes_mp[i*7+4].X,0.75)
            self.assertAlmostEqual(nodes_mp[i*7+5].X,0.916666666667)
            self.assertAlmostEqual(nodes_mp[i*7+6].X,1.0)
        for i in range(7):
            self.assertAlmostEqual(nodes_mp[i].Y,0.0)
            self.assertAlmostEqual(nodes_mp[i].Z,0.0)
        for i in range(7,14):
            self.assertAlmostEqual(nodes_mp[i].Y,0.5)
            self.assertAlmostEqual(nodes_mp[i].Z,0.0)
        for i in range(14,21):
            self.assertAlmostEqual(nodes_mp[i].Y,1.0)
            self.assertAlmostEqual(nodes_mp[i].Z,0.0)

        geometry = model_part.GetGeometry(1)

        # Check if geometry holds same nodes as model part
        for i, node_geo in enumerate(geometry):
            self.assertEqual(node_geo.Id, nodes_mp[i].Id )
            self.assertEqual(node_geo.X, nodes_mp[i].X )
            self.assertEqual(node_geo.Y, nodes_mp[i].Y )
            self.assertEqual(node_geo.Z, nodes_mp[i].Z )

        # test the creation of the breps
        self.assertEqual(current_model["IgaModelPart"].NumberOfGeometries(), 5)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().X, 0.5)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().Y, 0.5)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().Z, 0.0)
        

    # test for SBM
    def test_nurbs_geometry_2d_modeler_sbm_outer_skin_boundary(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsGeometryModelerSbm",
            "Parameters": {
                "model_part_name" : "IgaModelPart",
                "lower_point_xyz": [0.0,0.0,0.0],
                "upper_point_xyz": [2.0,2.0,0.0],
                "lower_point_uvw": [0.0,0.0,0.0],
                "upper_point_uvw": [2.0,2.0,0.0],
                "polynomial_order" : [1, 1],
                "number_of_knot_spans" : [5,5],
                "sbm_parameters": {
                    "lambda_outer": 0.5,
                    "number_of_inner_loops": 0
                },
                "skin_model_part_outer_initial_name": "skinModelPart_outer_initial",
                "skin_model_part_name": "skinModelPart"
            }
        }]
        """)
       # New model
        current_model = KratosMultiphysics.Model()
        skin_model_part_outer_initial = current_model.CreateModelPart("skinModelPart_outer_initial")

        skin_model_part_outer_initial.CreateNewProperties(1)

        skin_model_part_outer_initial.CreateNewNode(1, 0.0, 0.0, 0.0)
        skin_model_part_outer_initial.CreateNewNode(2, 2.0, 0.0, 0.0)
        skin_model_part_outer_initial.CreateNewNode(3, 2.0, 2.0, 0.0)
        skin_model_part_outer_initial.CreateNewNode(4, 0.0, 2.0, 0.0)

        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 1, [1, 2], skin_model_part_outer_initial.GetProperties()[1])
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 2, [2, 3], skin_model_part_outer_initial.GetProperties()[1])
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 3, [3, 4], skin_model_part_outer_initial.GetProperties()[1])
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 4, [4, 1], skin_model_part_outer_initial.GetProperties()[1])

        run_modelers(current_model, modeler_settings)
        
        # test surrogate nodes
        self.assertEqual(current_model["IgaModelPart.surrogate_outer"].GetNodes()[43].X, 0.4)
        self.assertEqual(current_model["IgaModelPart.surrogate_outer"].GetNodes()[43].Y, 2.0)
        self.assertEqual(current_model["IgaModelPart.surrogate_outer"].GetNodes()[43].Z, 0.0)

        self.assertEqual(current_model["IgaModelPart.surrogate_outer"].GetNodes()[51].X, 2.0)
        self.assertEqual(current_model["IgaModelPart.surrogate_outer"].GetNodes()[51].Y, 0.4)
        self.assertEqual(current_model["IgaModelPart.surrogate_outer"].GetNodes()[51].Z, 0.0)

        # test the creation of the breps
        self.assertEqual(current_model["IgaModelPart"].NumberOfGeometries(), 21)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().X, 1.0)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().Y, 1.0)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().Z, 0.0)

    
    # inner boundary single inner loop
    def test_nurbs_geometry_2d_modeler_sbm_inner_skin_boundary(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsGeometryModelerSbm",
            "Parameters": {
                "model_part_name" : "IgaModelPart",
                "lower_point_xyz": [0.0,0.0,0.0],
                "upper_point_xyz": [4.0,6.0,0.0],
                "lower_point_uvw": [0.0,0.0,0.0],
                "upper_point_uvw": [4.0,6.0,0.0],
                "polynomial_order" : [2, 2],
                "number_of_knot_spans" : [6,9],
                "sbm_parameters": {
                    "lambda_inner": 0.5,
                    "number_of_inner_loops": 1
                },
                "skin_model_part_inner_initial_name": "skinModelPart_inner_initial",
                "skin_model_part_name": "skinModelPart"
            }
        }]
        """)
       # New model
        current_model = KratosMultiphysics.Model()
        skin_model_part_inner_initial = current_model.CreateModelPart("skinModelPart_inner_initial")

        skin_model_part_inner_initial.CreateNewProperties(1)

        skin_model_part_inner_initial.CreateNewNode(1, 0.5, 0.5, 0.0)
        skin_model_part_inner_initial.CreateNewNode(2, 2.0, 0.5, 0.0)
        skin_model_part_inner_initial.CreateNewNode(3, 2.0, 3.5, 0.0)
        skin_model_part_inner_initial.CreateNewNode(4, 0.5, 3.5, 0.0)

        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 1, [1, 2], skin_model_part_inner_initial.GetProperties()[1])
        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 2, [2, 3], skin_model_part_inner_initial.GetProperties()[1])
        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 3, [3, 4], skin_model_part_inner_initial.GetProperties()[1])
        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 4, [4, 1], skin_model_part_inner_initial.GetProperties()[1])

        run_modelers(current_model, modeler_settings)

        # test surrogate nodes
        self.assertEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[92].X, 2/3)
        self.assertEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[92].Y, 8/3)
        self.assertEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[92].Z, 0.0)

        self.assertEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[98].X, 2.0)
        self.assertEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[98].Y, 4/3)
        self.assertEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[98].Z, 0.0)

        # test the creation of the breps
        self.assertEqual(current_model["IgaModelPart"].NumberOfGeometries(), 17)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().X, 2.0)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().Y, 3.0)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().Z, 0.0)


    # inner boundary double inner loops
    def test_nurbs_geometry_2d_modeler_sbm_two_inner_skin_boundary(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsGeometryModelerSbm",
            "Parameters": {
                "model_part_name" : "IgaModelPart",
                "lower_point_xyz": [0.0,0.0,0.0],
                "upper_point_xyz": [4.0,6.0,0.0],
                "lower_point_uvw": [0.0,0.0,0.0],
                "upper_point_uvw": [4.0,6.0,0.0],
                "polynomial_order" : [2, 2],
                "number_of_knot_spans" : [10,15],
                "sbm_parameters": {
                    "lambda_inner": 0.5,
                    "number_of_inner_loops": 2
                },
                "skin_model_part_inner_initial_name": "skinModelPart_inner_initial",
                "skin_model_part_name": "skinModelPart"
            }
        }]
        """)

        current_model = KratosMultiphysics.Model()

        skin_model_part_inner_initial = current_model.CreateModelPart("skinModelPart_inner_initial")
        skin_model_part_inner_initial.CreateNewProperties(1)

        # object 1
        skin_model_part_inner_initial.CreateNewNode(1, 0.5, 0.5, 0.0)
        skin_model_part_inner_initial.CreateNewNode(2, 2.0, 0.5, 0.0)
        skin_model_part_inner_initial.CreateNewNode(3, 2.0, 3.5, 0.0)
        skin_model_part_inner_initial.CreateNewNode(4, 0.5, 3.5, 0.0)

        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 1, [1, 2], skin_model_part_inner_initial.GetProperties()[1])
        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 2, [2, 3], skin_model_part_inner_initial.GetProperties()[1])
        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 3, [3, 4], skin_model_part_inner_initial.GetProperties()[1])
        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 4, [4, 1], skin_model_part_inner_initial.GetProperties()[1])

        # object 2
        skin_model_part_inner_initial.CreateNewNode(5, 1.5, 4.6, 0.0)
        skin_model_part_inner_initial.CreateNewNode(6, 3.7, 4.7, 0.0)
        skin_model_part_inner_initial.CreateNewNode(7, 3.0, 5.5, 0.0)

        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 5, [5, 6], skin_model_part_inner_initial.GetProperties()[1])
        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 6, [6, 7], skin_model_part_inner_initial.GetProperties()[1])
        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 7, [7, 5], skin_model_part_inner_initial.GetProperties()[1])

        run_modelers(current_model, modeler_settings)

        # test surrogate nodes
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[226].X, 1.6)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[226].Y, 0.4)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[226].Z, 0.0)

        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[235].X, 3.2)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[235].Y, 5.2)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[235].Z, 0.0)

        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[213].X, 0.4)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[213].Y, 3.6)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[213].Z, 0.0)

        # test the creation of the breps
        self.assertEqual(current_model["IgaModelPart"].NumberOfGeometries(), 41)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().X, 2.0)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().Y, 3.0)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().Z, 0.0)

    
    # inner boundary single inner loop + outer
    def test_nurbs_geometry_2d_modeler_sbm_inner_outer_skin_boundary(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsGeometryModelerSbm",
            "Parameters": {
                "model_part_name" : "IgaModelPart",
                "lower_point_xyz": [0.0,0.0,0.0],
                "upper_point_xyz": [4.0,6.0,0.0],
                "lower_point_uvw": [0.0,0.0,0.0],
                "upper_point_uvw": [4.0,6.0,0.0],
                "polynomial_order" : [2, 2],
                "number_of_knot_spans" : [20,10],
                "sbm_parameters": {
                    "lambda_inner": 1.0,
                    "lambda_outer": 0.5,
                    "number_of_inner_loops": 1
                },
                "skin_model_part_inner_initial_name": "skinModelPart_inner_initial",
                "skin_model_part_outer_initial_name": "skinModelPart_outer_initial",
                "skin_model_part_name": "skinModelPart"
            }
        }]
        """)
        current_model = KratosMultiphysics.Model()

        # Create Inner skin boundary

        skin_model_part_inner_initial = current_model.CreateModelPart("skinModelPart_inner_initial")
        skin_model_part_inner_initial.CreateNewProperties(1)
        skin_model_part_inner_initial.CreateNewNode(1, 1.5, 0.5, 0.0)
        skin_model_part_inner_initial.CreateNewNode(2, 3.0, 0.5, 0.0)
        skin_model_part_inner_initial.CreateNewNode(3, 3.0, 3.5, 0.0)
        skin_model_part_inner_initial.CreateNewNode(4, 1.5, 3.5, 0.0)
        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 1, [1, 2], skin_model_part_inner_initial.GetProperties()[1])
        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 2, [2, 3], skin_model_part_inner_initial.GetProperties()[1])
        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 3, [3, 4], skin_model_part_inner_initial.GetProperties()[1])
        skin_model_part_inner_initial.CreateNewCondition("LineCondition2D2N", 4, [4, 1], skin_model_part_inner_initial.GetProperties()[1])

        # Create Outer skin boundary
        skin_model_part_outer_initial = current_model.CreateModelPart("skinModelPart_outer_initial")
        skin_model_part_outer_initial.CreateNewProperties(1)
        skin_model_part_outer_initial.CreateNewNode(1, 0.4, 0.0, 0.0)
        skin_model_part_outer_initial.CreateNewNode(2, 3.5, 0.0, 0.0)
        skin_model_part_outer_initial.CreateNewNode(3, 3.2, 5.0, 0.0)
        skin_model_part_outer_initial.CreateNewNode(4, 1.0, 5.5, 0.0)
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 1, [1, 2], skin_model_part_outer_initial.GetProperties()[1])
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 2, [2, 3], skin_model_part_outer_initial.GetProperties()[1])
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 3, [3, 4], skin_model_part_outer_initial.GetProperties()[1])
        skin_model_part_outer_initial.CreateNewCondition("LineCondition2D2N", 4, [4, 1], skin_model_part_outer_initial.GetProperties()[1])

        run_modelers(current_model, modeler_settings)

        # test surrogate nodes
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[281].X, 2.8)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[281].Y, 0.6)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[281].Z, 0.0)

        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[269].X, 1.6)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[269].Y, 3.0)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_inner"].GetNodes()[269].Z, 0.0)

        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_outer"].GetNodes()[295].X, 0.8)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_outer"].GetNodes()[295].Y, 3.6)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_outer"].GetNodes()[295].Z, 0.0)

        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_outer"].GetNodes()[320].X, 3.6)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_outer"].GetNodes()[320].Y, 0.6)
        self.assertAlmostEqual(current_model["IgaModelPart.surrogate_outer"].GetNodes()[320].Z, 0.0)

        # test the creation of the breps
        self.assertEqual(current_model["IgaModelPart"].NumberOfGeometries(), 73)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().X, 2.0)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().Y, 3.0)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().Z, 0.0)


if __name__ == '__main__':
    KratosUnittest.main()
