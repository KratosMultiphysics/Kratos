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

    def test_iga_modeler_outer_sbm(self):
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
        
        modeler_settings = KratosMultiphysics.Parameters("""
        [
            {
            "modeler_name": "NurbsGeometryModelerSbm",
            "Parameters": {
                    "model_part_name" : "IgaModelPart",
                    "lower_point_xyz": [0.0,0.0,0.0],
                    "upper_point_xyz": [2.0,2.0,0.0],
                    "lower_point_uvw": [0.0,0.0,0.0],
                    "upper_point_uvw": [2.0,2.0,0.0],
                    "polynomial_order" : [1, 1],
                    "number_of_knot_spans" : [5,5],
                    "lambda_outer": 0.5,
                    "number_of_inner_loops": 0,
                    "skin_model_part_outer_initial_name": "skinModelPart_outer_initial",           
                    "skin_model_part_name": "skin_model_part"
                }
            },
            {
                "modeler_name": "IgaModelerSbm",
                "Parameters": {
                    "echo_level": 0,
                    "skin_model_part_name": "skin_model_part",
                    "analysis_model_part_name": "IgaModelPart",
                    "element_condition_list": [
                        {
                            "geometry_type": "GeometrySurface",
                            "iga_model_part": "ComputationalDomain",
                            "type": "element",
                            "name": "LaplacianIGAElement",
                            "shape_function_derivatives_order": 3,
                            "variables": [
                            {
                                "variable_name": "BODY_FORCE",
                                "value": ["0.0", "0.0", "0.0"]
                            }]
                        },
                        {
                            "geometry_type": "SurfaceEdge",
                            "iga_model_part": "SBM_Support_outer",
                            "type": "condition",
                            "name": "SbmLaplacianConditionDirichlet",
                            "shape_function_derivatives_order": 3, 
                            "sbm_parameters": {
                                "is_inner" : false
                            }
                        }
                    ] // element condition list
                }
            }] // iga modeler
            """)
        
        iga_model_part = current_model.CreateModelPart("IgaModelPart")
        iga_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        
        run_modelers(current_model, modeler_settings)

        support_model_part = current_model.GetModelPart("IgaModelPart.SBM_Support_outer")
        computational_model_part = current_model.GetModelPart("IgaModelPart.ComputationalDomain")

        # # Check if all needed node are within the model parts
        self.assertEqual(support_model_part.NumberOfNodes(), 240)
        self.assertEqual(support_model_part.NumberOfConditions(), 60)
        self.assertEqual(computational_model_part.NumberOfNodes(), 400)
        self.assertEqual(computational_model_part.NumberOfConditions(), 0)
        
        self.assertEqual(support_model_part.GetNodes()[6].X, 2.0)
        self.assertEqual(support_model_part.GetNodes()[6].Y, 0.0)
        self.assertEqual(support_model_part.GetNodes()[12].X, 2.0)
        self.assertEqual(support_model_part.GetNodes()[12].Y, 0.4)

        self.assertEqual(support_model_part.GetConditions()[21].Info(), "\"SbmLaplacianConditionDirichlet\" #21")
        self.assertEqual(support_model_part.GetConditions()[80].Info(), "\"SbmLaplacianConditionDirichlet\" #80")

        self.assertEqual(computational_model_part.GetElements()[13].Info(), "LaplacianIgaElement #13")
        self.assertEqual(computational_model_part.GetElements()[40].Info(), "LaplacianIgaElement #40")

    
    def test_iga_modeler_inner_outer_sbm(self):
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

        
        modeler_settings = KratosMultiphysics.Parameters("""
        [
            {
                "modeler_name": "NurbsGeometryModelerSbm",
                "Parameters": {
                        "model_part_name" : "IgaModelPart",
                        "lower_point_xyz": [0.0,0.0,0.0],
                        "upper_point_xyz": [4.0,6.0,0.0],
                        "lower_point_uvw": [0.0,0.0,0.0],
                        "upper_point_uvw": [4.0,6.0,0.0],
                        "polynomial_order" : [2, 2],
                        "number_of_knot_spans" : [20,10],
                        "lambda_inner": 1.0,
                        "lambda_outer": 0.5,
                        "number_of_inner_loops": 1,
                        "skin_model_part_inner_initial_name": "skinModelPart_inner_initial",
                        "skin_model_part_outer_initial_name": "skinModelPart_outer_initial",
                        "skin_model_part_name": "skinModelPart"
                }
            },
            {
                "modeler_name": "IgaModelerSbm",
                "Parameters": {
                    "echo_level": 0,
                    "skin_model_part_name": "skinModelPart",
                    "analysis_model_part_name": "IgaModelPart",
                    "element_condition_list": [
                        {
                            "geometry_type": "GeometrySurface",
                            "iga_model_part": "ComputationalDomain",
                            "type": "element",
                            "name": "LaplacianIGAElement",
                            "shape_function_derivatives_order": 3,
                            "variables": [
                            {
                                "variable_name": "BODY_FORCE",
                                "value": ["0.0", "0.0", "0.0"]
                            }]
                        },
                        {
                            "geometry_type": "SurfaceEdge",
                            "iga_model_part": "SBM_Support_outer",
                            "type": "condition",
                            "name": "SbmLaplacianConditionNeumann",
                            "shape_function_derivatives_order": 3, 
                            "sbm_parameters": {
                                "is_inner" : false
                            }
                        },
                        {
                            "geometry_type": "SurfaceEdge",
                            "iga_model_part": "SBM_Support_inner",
                            "type": "condition",
                            "name": "SbmLaplacianConditionDirichlet",
                            "shape_function_derivatives_order": 4, 
                            "sbm_parameters": {
                                "is_inner" : true
                            }
                        }
                    ] // element condition list
                }
            }
        ] // iga modeler
        """)

        iga_model_part = current_model.CreateModelPart("IgaModelPart")
        iga_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        
        run_modelers(current_model, modeler_settings)

        support_model_part_outer = current_model.GetModelPart("IgaModelPart.SBM_Support_outer")
        support_model_part_inner = current_model.GetModelPart("IgaModelPart.SBM_Support_inner")
        computational_model_part = current_model.GetModelPart("IgaModelPart.ComputationalDomain")

        # # Check if all needed node are within the model parts
        self.assertEqual(support_model_part_inner.NumberOfNodes(), 990)
        self.assertEqual(support_model_part_inner.NumberOfConditions(), 110)
        self.assertEqual(support_model_part_outer.NumberOfNodes(), 2250)
        self.assertEqual(support_model_part_outer.NumberOfConditions(), 250)
        self.assertEqual(computational_model_part.NumberOfNodes(), 7371)
        self.assertEqual(computational_model_part.NumberOfConditions(), 0)
        self.assertEqual(computational_model_part.NumberOfElements(), 819)

        self.assertEqual(support_model_part_inner.GetConditions()[323].Info(), "\"SbmLaplacianConditionDirichlet\" #323")
        self.assertEqual(support_model_part_inner.GetConditions()[432].Info(), "\"SbmLaplacianConditionDirichlet\" #432")
        self.assertEqual(support_model_part_outer.GetConditions()[73].Info(), "\"SbmLaplacianConditionNeumann\" #73")
        self.assertEqual(support_model_part_outer.GetConditions()[322].Info(), "\"SbmLaplacianConditionNeumann\" #322")
        self.assertEqual(computational_model_part.GetElements()[13].Info(), "LaplacianIgaElement #13")
        self.assertEqual(computational_model_part.GetElements()[40].Info(), "LaplacianIgaElement #40")

    
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
                "lambda_outer": 0.5,
                "number_of_inner_loops": 0,
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

        expected_nodes = [
            [0, 0.4, 0],
            [2, 0.4, 0],
            [0, 0.8, 0],
            [2, 0.8, 0],
            [0, 1.2, 0],
            [2, 1.2, 0],
            [0, 1.6, 0],
            [2, 1.6, 0],
            [0, 2, 0],
            [2, 2, 0],
            [0.4, 0, 0],
            [0.4, 2, 0],
            [0.8, 0, 0],
            [0.8, 2, 0],
            [1.2, 0, 0],
            [1.2, 2, 0],
            [1.6, 0, 0],
            [1.6, 2, 0],
            [2, 0, 0],
            [2, 2, 0]
        ]

        # test surrogate nodes
        i = 0
        for cond in current_model["IgaModelPart.surrogate_outer"].Conditions:
            expected_boundary = (i%2==0) #should be alternate between true and false
            self.assertAlmostEqual(cond.GetNodes()[1].X, expected_nodes[i][0])
            self.assertAlmostEqual(cond.GetNodes()[1].Y, expected_nodes[i][1])
            self.assertAlmostEqual(cond.GetNodes()[1].Z, expected_nodes[i][2])
            self.assertEqual(cond.Is(KratosMultiphysics.BOUNDARY), expected_boundary)
            i += 1
        
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
                "lambda_inner": 0.5,
                "number_of_inner_loops": 1,
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

        expected_nodes = [
            [2/3, 4/3, 0],
            [2, 4/3, 0],
            [2/3, 2, 0],
            [2, 2, 0],
            [2/3, 8/3, 0],
            [2, 8/3, 0],
            [2/3, 10/3, 0],
            [2, 10/3, 0],
            [4/3, 2/3, 0],
            [4/3, 10/3, 0],
            [2, 2/3, 0],
            [2, 10/3, 0]
        ]
        # test surrogate nodes
        i = 0
        for cond in current_model["IgaModelPart.surrogate_inner"].Conditions:
            expected_boundary = (i%2==0) #should be alternate between true and false
            self.assertAlmostEqual(cond.GetNodes()[1].X, expected_nodes[i][0])
            self.assertAlmostEqual(cond.GetNodes()[1].Y, expected_nodes[i][1])
            self.assertAlmostEqual(cond.GetNodes()[1].Z, expected_nodes[i][2])
            self.assertEqual(cond.Is(KratosMultiphysics.BOUNDARY), expected_boundary)
            i += 1

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
                "lambda_inner": 0.5,
                "number_of_inner_loops": 2,
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
        
        expected_nodes = [
            [0.4, 0.8, 0],
            [2, 0.8, 0],
            [0.4, 1.2, 0],
            [2, 1.2, 0],
            [0.4, 1.6, 0],
            [2, 1.6, 0],
            [0.4, 2, 0],
            [2, 2, 0],
            [0.4, 2.4, 0],
            [2, 2.4, 0],
            [0.4, 2.8, 0],
            [2, 2.8, 0],
            [0.4, 3.2, 0],
            [2, 3.2, 0],
            [0.4, 3.6, 0],
            [2, 3.6, 0],
            [0.8, 0.4, 0],
            [0.8, 3.6, 0],
            [1.2, 0.4, 0],
            [1.2, 3.6, 0],
            [1.6, 0.4, 0],
            [1.6, 3.6, 0],
            [2, 0.4, 0],
            [2, 3.6, 0],
            [2, 5.2, 0],
            [3.6, 5.2, 0],
            [2.8, 5.6, 0],
            [3.2, 5.6, 0],
            [2.4, 4.8, 0],
            [2.4, 5.2, 0],
            [2.8, 4.8, 0],
            [2.8, 5.2, 0],
            [3.2, 4.8, 0],
            [3.2, 5.6, 0],
            [3.6, 4.8, 0],
            [3.6, 5.2, 0]
        ]

        # test surrogate nodes
        i = 0
        for cond in current_model["IgaModelPart.surrogate_inner"].Conditions:
            expected_boundary = (i%2==0) #should be alternate between true and false
            self.assertAlmostEqual(cond.GetNodes()[1].X, expected_nodes[i][0])
            self.assertAlmostEqual(cond.GetNodes()[1].Y, expected_nodes[i][1])
            self.assertAlmostEqual(cond.GetNodes()[1].Z, expected_nodes[i][2])
            self.assertEqual(cond.Is(KratosMultiphysics.BOUNDARY), expected_boundary)
            i += 1

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
                "lambda_inner": 1.0,
                "lambda_outer": 0.5,
                "number_of_inner_loops": 1,
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

        expected_nodes_inner = [
            [1.6, 1.2, 0],
            [3, 1.2, 0],
            [1.6, 1.8, 0],
            [3, 1.8, 0],
            [1.6, 2.4, 0],
            [3, 2.4, 0],
            [1.6, 3, 0],
            [3, 3, 0],
            [1.8, 0.6, 0],
            [1.8, 3, 0],
            [2, 0.6, 0],
            [2, 3, 0],
            [2.2, 0.6, 0],
            [2.2, 3, 0],
            [2.4, 0.6, 0],
            [2.4, 3, 0],
            [2.6, 0.6, 0],
            [2.6, 3, 0],
            [2.8, 0.6, 0],
            [2.8, 3, 0],
            [3, 0.6, 0],
            [3, 3, 0]
        ]

        expected_nodes_outer = [
            [0.4, 0.6, 0],
            [3.6, 0.6, 0],
            [0.4, 1.2, 0],
            [3.4, 1.2, 0],
            [0.6, 1.8, 0],
            [3.4, 1.8, 0],
            [0.6, 2.4, 0],
            [3.4, 2.4, 0],
            [0.6, 3, 0],
            [3.4, 3, 0],
            [0.8, 3.6, 0],
            [3.4, 3.6, 0],
            [0.8, 4.2, 0],
            [3.2, 4.2, 0],
            [0.8, 4.8, 0],
            [3.2, 4.8, 0],
            [1, 5.4, 0],
            [2.8, 5.4, 0],
            [0.6, 0, 0],
            [0.6, 1.2, 0],
            [0.8, 0, 0],
            [0.8, 3, 0],
            [1, 0, 0],
            [1, 4.8, 0],
            [1.2, 0, 0],
            [1.2, 5.4, 0],
            [1.4, 0, 0],
            [1.4, 5.4, 0],
            [1.6, 0, 0],
            [1.6, 5.4, 0],
            [1.8, 0, 0],
            [1.8, 5.4, 0],
            [2, 0, 0],
            [2, 5.4, 0],
            [2.2, 0, 0],
            [2.2, 5.4, 0],
            [2.4, 0, 0],
            [2.4, 5.4, 0],
            [2.6, 0, 0],
            [2.6, 5.4, 0],
            [2.8, 0, 0],
            [2.8, 5.4, 0],
            [3, 0, 0],
            [3, 4.8, 0],
            [3.2, 0, 0],
            [3.2, 4.8, 0],
            [3.4, 0, 0],
            [3.4, 3.6, 0],
            [3.6, 0, 0],
            [3.6, 0.6, 0]
        ]

        # test surrogate nodes
        i = 0
        for cond in current_model["IgaModelPart.surrogate_inner"].Conditions:
            expected_boundary = (i%2==0) #should be alternate between true and false
            self.assertAlmostEqual(cond.GetNodes()[1].X, expected_nodes_inner[i][0])
            self.assertAlmostEqual(cond.GetNodes()[1].Y, expected_nodes_inner[i][1])
            self.assertAlmostEqual(cond.GetNodes()[1].Z, expected_nodes_inner[i][2])
            self.assertEqual(cond.Is(KratosMultiphysics.BOUNDARY), expected_boundary)
            i += 1
            
        i = 0
        for cond in current_model["IgaModelPart.surrogate_outer"].Conditions:
            expected_boundary = (i%2==0) #should be alternate between true and false
            self.assertAlmostEqual(cond.GetNodes()[1].X, expected_nodes_outer[i][0])
            self.assertAlmostEqual(cond.GetNodes()[1].Y, expected_nodes_outer[i][1])
            self.assertAlmostEqual(cond.GetNodes()[1].Z, expected_nodes_outer[i][2])
            self.assertEqual(cond.Is(KratosMultiphysics.BOUNDARY), expected_boundary)
            i += 1

        # test the creation of the breps
        self.assertEqual(current_model["IgaModelPart"].NumberOfGeometries(), 73)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().X, 2.0)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().Y, 3.0)
        self.assertAlmostEqual(current_model["IgaModelPart"].GetGeometry(1).Center().Z, 0.0)


if __name__ == '__main__':
    KratosUnittest.main()
