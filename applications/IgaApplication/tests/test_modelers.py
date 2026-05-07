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

class TestModelers(KratosUnittest.TestCase):
    def test_StrongSupportSurface(self):
        current_model = KratosMultiphysics.Model()
        modelers_list = KratosMultiphysics.Parameters(
        """ [{
            "modeler_name": "CadIoModeler",
            "Parameters": {
                "echo_level": 0,
                "cad_model_part_name": "IgaModelPart",
                "geometry_file_name": "modeler_tests/surface_geometry.cad.json"
            } }, {
            "modeler_name": "IgaModeler",
            "Parameters": {
                "echo_level":  0,
                "cad_model_part_name": "IgaModelPart",
                "analysis_model_part_name": "IgaModelPart",
                "physics_file_name": "modeler_tests/strong_support_physics.iga.json"
            } }] """ )

        run_modelers(current_model, modelers_list)

        support_1_model_part = current_model.GetModelPart("IgaModelPart.Support_1")
        support_1_Variation_model_part = current_model.GetModelPart("IgaModelPart.Support_1_Variation")
        support_2_model_part = current_model.GetModelPart("IgaModelPart.Support_2")
        support_2_Variation_model_part = current_model.GetModelPart("IgaModelPart.Support_2_Variation")

        # Check if all needed node are within the model parts
        self.assertEqual(support_1_model_part.NumberOfNodes(), 2)
        self.assertEqual(support_1_model_part.GetNodes()[6].Id, 6)
        self.assertEqual(support_1_model_part.GetNodes()[12].Id, 12)

        self.assertEqual(support_1_Variation_model_part.NumberOfNodes(), 2)
        self.assertEqual(support_1_Variation_model_part.GetNodes()[5].Id, 5)
        self.assertEqual(support_1_Variation_model_part.GetNodes()[11].Id, 11)

        self.assertEqual(support_2_model_part.NumberOfNodes(), 6)
        self.assertEqual(support_2_model_part.GetNodes()[7].Id, 7)
        self.assertEqual(support_2_model_part.GetNodes()[8].Id, 8)
        self.assertEqual(support_2_model_part.GetNodes()[9].Id, 9)
        self.assertEqual(support_2_model_part.GetNodes()[10].Id, 10)
        self.assertEqual(support_2_model_part.GetNodes()[11].Id, 11)
        self.assertEqual(support_2_model_part.GetNodes()[12].Id, 12)

        self.assertEqual(support_2_Variation_model_part.NumberOfNodes(), 6)
        self.assertEqual(support_2_Variation_model_part.GetNodes()[1].Id, 1)
        self.assertEqual(support_2_Variation_model_part.GetNodes()[2].Id, 2)
        self.assertEqual(support_2_Variation_model_part.GetNodes()[3].Id, 3)
        self.assertEqual(support_2_Variation_model_part.GetNodes()[4].Id, 4)
        self.assertEqual(support_2_Variation_model_part.GetNodes()[5].Id, 5)
        self.assertEqual(support_2_Variation_model_part.GetNodes()[6].Id, 6)

    def test_RefinementModeler_mesh(self):
        current_model = KratosMultiphysics.Model()
        modelers_list = KratosMultiphysics.Parameters(
        """ [{
            "modeler_name": "CadIoModeler",
            "Parameters": {
                "echo_level": 0,
                "cad_model_part_name": "IgaModelPart",
                "geometry_file_name": "modeler_tests/surface_geometry.cad.json"
            } }, {
            "modeler_name": "RefinementModeler",
            "Parameters": {
                "echo_level":  0,
                "refinements_file_name": "modeler_tests/h_refinements.iga.json"
            } }] """ )

        run_modelers(current_model, modelers_list)

        model_part = current_model.GetModelPart("IgaModelPart")

        # Check if all needed node are within the model parts
        self.assertEqual(model_part.NumberOfNodes(), 42)
        self.assertEqual(model_part.GetGeometry(1).GetGeometryPart(KratosMultiphysics.Geometry.BACKGROUND_GEOMETRY_INDEX).PointsNumber(), 42)
    
    def test_RefinementModeler_polynomial(self):
        current_model = KratosMultiphysics.Model()
        modelers_list = KratosMultiphysics.Parameters(
        """ [{
            "modeler_name": "CadIoModeler",
            "Parameters": {
                "echo_level": 0,
                "cad_model_part_name": "IgaModelPart",
                "geometry_file_name": "modeler_tests/surface_geometry.cad.json"
            } }, {
            "modeler_name": "RefinementModeler",
            "Parameters": {
                "echo_level":  0,
                "refinements_file_name": "modeler_tests/p_refinements.iga.json"
            } }] """ )

        run_modelers(current_model, modelers_list)

        model_part = current_model.GetModelPart("IgaModelPart")

        # Check if all needed node are within the model parts
        self.assertEqual(model_part.NumberOfNodes(), 30)
        self.assertEqual(model_part.GetGeometry(1).GetGeometryPart(KratosMultiphysics.Geometry.BACKGROUND_GEOMETRY_INDEX).PointsNumber(), 30)
    
    def test_RefinementModeler_both(self):
        current_model = KratosMultiphysics.Model()
        modelers_list = KratosMultiphysics.Parameters(
        """ [{
            "modeler_name": "CadIoModeler",
            "Parameters": {
                "echo_level": 0,
                "cad_model_part_name": "IgaModelPart",
                "geometry_file_name": "modeler_tests/surface_geometry.cad.json"
            } }, {
            "modeler_name": "RefinementModeler",
            "Parameters": {
                "echo_level":  0,
                "refinements_file_name": "modeler_tests/k_refinements.iga.json"
            } }] """ )

        run_modelers(current_model, modelers_list)

        model_part = current_model.GetModelPart("IgaModelPart")

        # Check if all needed node are within the model parts
        self.assertEqual(model_part.NumberOfNodes(), 72)
        self.assertEqual(model_part.GetGeometry(1).GetGeometryPart(KratosMultiphysics.Geometry.BACKGROUND_GEOMETRY_INDEX).PointsNumber(), 72)

    def test_RefinementModeler_multipatch(self):
        current_model = KratosMultiphysics.Model()
        modelers_list = KratosMultiphysics.Parameters(
        """ [{
            "modeler_name": "CadIoModeler",
            "Parameters": {
                "echo_level": 0,
                "cad_model_part_name": "IgaModelPart",
                "geometry_file_name": "modeler_tests/multipatch_geometry.cad.json"
            } }, {
            "modeler_name": "RefinementModeler",
            "Parameters": {
                "echo_level":  0,
                "refinements_file_name": "modeler_tests/multipatch_refinements.iga.json"
            } }] """ )

        run_modelers(current_model, modelers_list)

        model_part = current_model.GetModelPart("IgaModelPart")

        # Check if all needed node are within the model parts
        self.assertEqual(model_part.NumberOfNodes(), 50)
        self.assertEqual(model_part.GetGeometry(2).GetGeometryPart(KratosMultiphysics.Geometry.BACKGROUND_GEOMETRY_INDEX).PointsNumber(), 30)
        self.assertEqual(model_part.GetGeometry(3).GetGeometryPart(KratosMultiphysics.Geometry.BACKGROUND_GEOMETRY_INDEX).PointsNumber(), 20)

    def test_nurbs_geometry_2d_modeler_control_points(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsGeometryModeler",
            "Parameters": {
                "model_part_name" : "Mesh",
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

        model_part = current_model.GetModelPart("Mesh")
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

        # Check knots
        knots_u = geometry.KnotsU()
        self.assertVectorAlmostEqual(knots_u, [0.0, 0.0, 0.0, 0.0, 1.0/3.0, 2.0/3.0, 1.0, 1.0, 1.0, 1.0])
        knots_v = geometry.KnotsV()
        self.assertVectorAlmostEqual(knots_v, [0.0, 0.5,  1.0 ])

        # # Check polynomial degree
        degree_u = geometry.PolynomialDegreeU()
        self.assertEqual(degree_u, 4)
        degree_v = geometry.PolynomialDegreeV()
        self.assertEqual(degree_v, 1)

    def test_nurbs_geometry_2d_modeler_consistency(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsGeometryModeler",
            "Parameters": {
                "model_part_name" : "Mesh_01",
                "lower_point_xyz": [-1.0,-0.5,0.0],
                "upper_point_xyz": [1.0, 0.5,0.0],
                "lower_point_uvw": [0.0, 0.0, 0.0],
                "upper_point_uvw": [1.0, 1.0, 0.0],
                "polynomial_order" : [4, 5],
                "number_of_knot_spans" : [3,6]
            } }, {
            "modeler_name": "NurbsGeometryModeler",
            "Parameters": {
                "model_part_name" : "Mesh_02",
                "lower_point_xyz": [-1.0,-0.5,0.0],
                "upper_point_xyz": [1.0,0.5,0.0],
                "lower_point_uvw": [0.0, 0.0, 0.0],
                "upper_point_uvw": [1.0, 1.0, 0.0],
                "polynomial_order" : [1, 1],
                "number_of_knot_spans" : [1,1]
            }
        }]
        """)

        run_modelers(current_model, modeler_settings)
        model_part_01 = current_model.GetModelPart("Mesh_01")
        model_part_02 = current_model.GetModelPart("Mesh_02")

        geometry_01 = model_part_01.GetGeometry(1)
        geometry_02 = model_part_02.GetGeometry(1)

        # Check if geometry holds same nodes as model part
        for node_geo, node_mp in zip(geometry_01, model_part_01.Nodes):
            self.assertEqual(node_geo.Id, node_mp.Id )
            self.assertAlmostEqual(node_geo.X, node_mp.X)
            self.assertAlmostEqual(node_geo.Y, node_mp.Y)
            self.assertAlmostEqual(node_geo.Z, node_mp.Z)
        for node_geo, node_mp in zip(geometry_01, model_part_01.Nodes):
            self.assertEqual(node_geo.Id, node_mp.Id )
            self.assertAlmostEqual(node_geo.X, node_mp.X)
            self.assertAlmostEqual(node_geo.Y, node_mp.Y)
            self.assertAlmostEqual(node_geo.Z, node_mp.Z)

        ## Check knots
        # geometry 01
        knots_u_01 = geometry_01.KnotsU()
        self.assertVectorAlmostEqual(knots_u_01,
            [0.0, 0.0, 0.0, 0.0, 1.0/3.0, 2.0/3.0, 1.0, 1.0, 1.0, 1.0])
        knots_v_01 = geometry_01.KnotsV()
        self.assertVectorAlmostEqual(knots_v_01,
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0/6.0, 1.0/3.0, 1.0/2.0, 2.0/3.0, 5.0/6.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        # geometry 02
        knots_u_02 = geometry_02.KnotsU()
        self.assertVectorAlmostEqual(knots_u_02, [0.0, 1.0])
        knots_v_02 = geometry_02.KnotsV()
        self.assertVectorAlmostEqual(knots_v_02, [0.0, 1.0])

        ## Check Polynomial degree
        # Geometry 01
        self.assertEqual(geometry_01.PolynomialDegreeU(), 4)
        self.assertEqual(geometry_01.PolynomialDegreeV(), 5)
        # Geometry 02
        self.assertEqual(geometry_02.PolynomialDegreeU(), 1)
        self.assertEqual(geometry_02.PolynomialDegreeV(), 1)

        ## Check number of control points/nodes
        # Geomtry 01
        self.assertEqual(len(geometry_01),77)
        self.assertEqual(model_part_01.NumberOfNodes(), 77)
        # Geomtry 02
        self.assertEqual(len(geometry_02),4)
        self.assertEqual(model_part_02.NumberOfNodes(), 4)

        # Check if geoemtry is similar.
        param = KratosMultiphysics.Vector(3)
        param[0] = 0.0
        param[1] = 1.0
        param[2] = 0.0
        for _ in range(11):
            for _ in range(11):
                coord_01 = geometry_01.GlobalCoordinates(param)
                coord_02 = geometry_02.GlobalCoordinates(param)
                param[2] += 0.1
                self.assertVectorAlmostEqual(coord_01, coord_02)
            param[0] += 0.1
            param[1] -= 0.1
            param[2] = 0.0

    def test_nurbs_geometry_2d_modeler_parameter_space(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsGeometryModeler",
            "Parameters": {
                "model_part_name" : "Mesh_01",
                "lower_point_xyz": [-1.0,-0.5,0.0],
                "upper_point_xyz": [1.0, 0.5,0.0],
                "lower_point_uvw": [0.0, 0.0, 0.0],
                "upper_point_uvw": [1.0, 1.0, 0.0],
                "polynomial_order" : [4, 5],
                "number_of_knot_spans" : [3 ,6]
            } }, {
            "modeler_name": "NurbsGeometryModeler",
            "Parameters": {
                "model_part_name" : "Mesh_02",
                "lower_point_xyz": [-1.0,-0.5,0.0],
                "upper_point_xyz": [1.0,0.5,0.0],
                "lower_point_uvw": [-2.22, -5.5, 0.0],
                "upper_point_uvw": [5.0, 5.5, 0.0],
                "polynomial_order" : [4, 5],
                "number_of_knot_spans" : [3, 6]
            }
        }]
        """)

        run_modelers(current_model, modeler_settings)
        model_part_01 = current_model.GetModelPart("Mesh_01")
        model_part_02 = current_model.GetModelPart("Mesh_02")

        geometry_01 = model_part_01.GetGeometry(1)
        geometry_02 = model_part_02.GetGeometry(1)

        ## Check knots
        knots_u_01 = geometry_01.KnotsU()
        knots_u_02 = geometry_02.KnotsU()
        self.assertVectorAlmostEqual(knots_u_01,
            [0.0, 0.0, 0.0, 0.0, 1.0/3.0, 2.0/3.0, 1.0, 1.0, 1.0, 1.0])
        for u_01, u_02 in zip(knots_u_01, knots_u_02):
            u_02_mapped = (u_02 + 2.22) / (5.0 + 2.22)
            self.assertAlmostEqual(u_01, u_02_mapped)

        knots_v_01 = geometry_01.KnotsV()
        knots_v_02 = geometry_02.KnotsV()
        self.assertVectorAlmostEqual(knots_v_01,
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0/6.0, 1.0/3.0, 1.0/2.0, 2.0/3.0, 5.0/6.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        for v_01, v_02 in zip(knots_v_01, knots_v_02):
            v_02_mapped = (v_02 + 5.5) / (5.5 + 5.5)
            self.assertAlmostEqual(v_01, v_02_mapped)

        ## Check control points.
        for node_1, node_2 in zip(geometry_01, geometry_02):
            self.assertAlmostEqual(node_1.X, node_2.X)
            self.assertAlmostEqual(node_1.Y, node_2.Y)
            self.assertAlmostEqual(node_1.Z, node_2.Z)

        # Check if geoemtry is similar.
        param = KratosMultiphysics.Vector(3)
        param[0] = 0.0
        param[1] = 1.0
        param[2] = 0.0
        param2 = KratosMultiphysics.Vector(3)
        param2[2] = 0.0
        for _ in range(11):
            for _ in range(11):
                coord_01 = geometry_01.GlobalCoordinates(param)
                param2[0] = param[0]*(5.0 + 2.22) - 2.22
                param2[1] = param[1]*(5.5 + 5.5) - 5.5
                coord_02 = geometry_02.GlobalCoordinates(param2)
                param[2] += 0.1
                self.assertVectorAlmostEqual(coord_01, coord_02)
            param[0] += 0.1
            param[1] -= 0.1
            param[2] = 0.0

    def test_nurbs_geometry_3d_modeler_control_points(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsGeometryModeler",
            "Parameters": {
                "model_part_name" : "Mesh",
                "lower_point_xyz": [0.0,0.0,0.0],
                "upper_point_xyz": [1.0,1.0,1.0],
                "lower_point_uvw": [0.0,0.0,0.0],
                "upper_point_uvw": [1.0,1.0,1.0],
                "polynomial_order" : [4, 1, 1],
                "number_of_knot_spans" : [3,1,1]
            }
        }]
        """)
        run_modelers(current_model, modeler_settings)
        model_part = current_model.GetModelPart("Mesh")
        # Check control points
        nodes_mp = model_part.NodesArray(0)
        for i in range(4):
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
            self.assertAlmostEqual(nodes_mp[i].Y,1.0)
            self.assertAlmostEqual(nodes_mp[i].Z,0.0)
        for i in range(14,21):
            self.assertAlmostEqual(nodes_mp[i].Y,0.0)
            self.assertAlmostEqual(nodes_mp[i].Z,1.0)
        for i in range(21,28):
            self.assertAlmostEqual(nodes_mp[i].Y,1.0)
            self.assertAlmostEqual(nodes_mp[i].Z,1.0)

        geometry = model_part.GetGeometry(1)

        # Check if geometry holds same nodes as model part
        for i, node_geo in enumerate(geometry):
            self.assertEqual(node_geo.Id, nodes_mp[i].Id )
            self.assertEqual(node_geo.X, nodes_mp[i].X )
            self.assertEqual(node_geo.Y, nodes_mp[i].Y )
            self.assertEqual(node_geo.Z, nodes_mp[i].Z )

        # Check knots
        knots_u = geometry.KnotsU()
        self.assertVectorAlmostEqual(knots_u, [0.0, 0.0, 0.0, 0.0, 1.0/3.0, 2.0/3.0, 1.0, 1.0, 1.0, 1.0])
        knots_v = geometry.KnotsV()
        self.assertVectorAlmostEqual(knots_v, [0.0,  1.0 ])
        knots_w = geometry.KnotsW()
        self.assertVectorAlmostEqual(knots_w, [0.0, 1.0])
        # Check polynomial degree
        degree_u = geometry.PolynomialDegreeU()
        self.assertEqual(degree_u, 4)
        degree_v = geometry.PolynomialDegreeV()
        self.assertEqual(degree_v, 1)
        degree_w = geometry.PolynomialDegreeW()
        self.assertEqual(degree_w, 1)

    def test_nurbs_geometry_3d_modeler_consistency(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsGeometryModeler",
            "Parameters": {
                "model_part_name" : "Mesh_01",
                "lower_point_xyz": [-1.0,-0.5,1.0],
                "upper_point_xyz": [1.0,0.5,3.0],
                "lower_point_uvw": [0.0,0.0,0.0],
                "upper_point_uvw": [1.0,1.0,1.0],
                "polynomial_order" : [4, 5, 2],
                "number_of_knot_spans" : [3,6,7]
            } }, {
            "modeler_name": "NurbsGeometryModeler",
            "Parameters": {
                "model_part_name" : "Mesh_02",
                "lower_point_xyz": [-1.0,-0.5,1.0],
                "upper_point_xyz": [1.0,0.5,3.0],
                "lower_point_uvw": [0.0,0.0,0.0],
                "upper_point_uvw": [1.0,1.0,1.0],
                "polynomial_order" : [1, 1, 1],
                "number_of_knot_spans" : [1,1,1]
            }
        }]
        """)

        run_modelers(current_model, modeler_settings)
        model_part_01 = current_model.GetModelPart("Mesh_01")
        model_part_02 = current_model.GetModelPart("Mesh_02")

        geometry_01 = model_part_01.GetGeometry(1)
        geometry_02 = model_part_02.GetGeometry(1)

        # Check if geometry holds same nodes as model part
        for node_geo, node_mp in zip(geometry_01, model_part_01.Nodes):
            self.assertEqual(node_geo.Id, node_mp.Id )
            self.assertAlmostEqual(node_geo.X, node_mp.X)
            self.assertAlmostEqual(node_geo.Y, node_mp.Y)
            self.assertAlmostEqual(node_geo.Z, node_mp.Z)
        for node_geo, node_mp in zip(geometry_01, model_part_01.Nodes):
            self.assertEqual(node_geo.Id, node_mp.Id )
            self.assertAlmostEqual(node_geo.X, node_mp.X)
            self.assertAlmostEqual(node_geo.Y, node_mp.Y)
            self.assertAlmostEqual(node_geo.Z, node_mp.Z)

        ## Check knots
        # geometry 01
        knots_u_01 = geometry_01.KnotsU()
        self.assertVectorAlmostEqual(knots_u_01,
            [0.0, 0.0, 0.0, 0.0, 1.0/3.0, 2.0/3.0, 1.0, 1.0, 1.0, 1.0])
        knots_v_01 = geometry_01.KnotsV()
        self.assertVectorAlmostEqual(knots_v_01,
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0/6.0, 1.0/3.0, 1.0/2.0, 2.0/3.0, 5.0/6.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        knots_w_01 = geometry_01.KnotsW()
        self.assertVectorAlmostEqual(knots_w_01,
            [0.0, 0.0, 1.0/7.0, 2.0/7.0, 3.0/7.0, 4.0/7.0, 5.0/7.0, 6.0/7.0, 1.0, 1.0])
        # geometry 02
        knots_u_02 = geometry_02.KnotsU()
        self.assertVectorAlmostEqual(knots_u_02, [0.0, 1.0])
        knots_v_02 = geometry_02.KnotsV()
        self.assertVectorAlmostEqual(knots_v_02, [0.0, 1.0])
        knots_w_02 = geometry_02.KnotsW()
        self.assertVectorAlmostEqual(knots_w_02, [0.0, 1.0])

        ## Check Polynomial degree
        # Geometry 01
        self.assertEqual(geometry_01.PolynomialDegreeU(), 4)
        self.assertEqual(geometry_01.PolynomialDegreeV(), 5)
        self.assertEqual(geometry_01.PolynomialDegreeW(), 2)
        # Geometry 02
        self.assertEqual(geometry_02.PolynomialDegreeU(), 1)
        self.assertEqual(geometry_02.PolynomialDegreeV(), 1)
        self.assertEqual(geometry_02.PolynomialDegreeW(), 1)

        ## Check number of control points/nodes
        # Geomtry 01
        self.assertEqual(len(geometry_01),693)
        self.assertEqual(model_part_01.NumberOfNodes(), 693)
        # Geomtry 02
        self.assertEqual(len(geometry_02),8)
        self.assertEqual(model_part_02.NumberOfNodes(), 8)

        # Check if geoemtry is similar.
        param = KratosMultiphysics.Vector(3)
        param[0] = 0.0
        param[1] = 0.0
        param[2] = 0.0
        for _ in range(11):
            for _ in range(11):
                for _ in range(11):
                    coord_01 = geometry_01.GlobalCoordinates(param)
                    coord_02 = geometry_02.GlobalCoordinates(param)
                    param[2] += 0.1
                    self.assertVectorAlmostEqual(coord_01, coord_02)
                param[1] += 0.1
                param[2] = 0.0
            param[0] += 0.1
            param[1] = 0.0

    def test_nurbs_geometry_3d_modeler_parameter_space(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsGeometryModeler",
            "Parameters": {
                "model_part_name" : "Mesh_01",
                "lower_point_xyz": [-1.0,-0.5,1.0],
                "upper_point_xyz": [1.0,0.5,3.0],
                "lower_point_uvw": [0.0, 0.0, 0.0],
                "upper_point_uvw": [1.0, 1.0, 1.0],
                "polynomial_order" : [4, 5, 2],
                "number_of_knot_spans" : [3,6,7]
            } }, {
            "modeler_name": "NurbsGeometryModeler",
            "Parameters": {
                "model_part_name" : "Mesh_02",
                "lower_point_xyz": [-1.0, -0.5, 1.0],
                "upper_point_xyz": [1.0, 0.5, 3.0],
                "lower_point_uvw": [-2.22, -5.5, 1.0],
                "upper_point_uvw": [5.0, 5.5, 3.0],
                "polynomial_order" : [4, 5, 2],
                "number_of_knot_spans" : [3,6,7]
            }
        }]
        """)

        run_modelers(current_model, modeler_settings)
        model_part_01 = current_model.GetModelPart("Mesh_01")
        model_part_02 = current_model.GetModelPart("Mesh_02")

        geometry_01 = model_part_01.GetGeometry(1)
        geometry_02 = model_part_02.GetGeometry(1)

        ## Check knots
        knots_u_01 = geometry_01.KnotsU()
        knots_u_02 = geometry_02.KnotsU()
        self.assertVectorAlmostEqual(knots_u_01,
            [0.0, 0.0, 0.0, 0.0, 1.0/3.0, 2.0/3.0, 1.0, 1.0, 1.0, 1.0])
        for u_01, u_02 in zip(knots_u_01, knots_u_02):
            u_02_mapped = (u_02 + 2.22) / (5.0 + 2.22)
            self.assertAlmostEqual(u_01, u_02_mapped)

        knots_v_01 = geometry_01.KnotsV()
        knots_v_02 = geometry_02.KnotsV()
        self.assertVectorAlmostEqual(knots_v_01,
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0/6.0, 1.0/3.0, 1.0/2.0, 2.0/3.0, 5.0/6.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        for v_01, v_02 in zip(knots_v_01, knots_v_02):
            v_02_mapped = (v_02 + 5.5) / (5.5 + 5.5)
            self.assertAlmostEqual(v_01, v_02_mapped)

        knots_w_01 = geometry_01.KnotsW()
        knots_w_02 = geometry_02.KnotsW()
        self.assertVectorAlmostEqual(knots_w_01,
            [0.0, 0.0, 1.0/7.0, 2.0/7.0, 3.0/7.0, 4.0/7.0, 5.0/7.0, 6.0/7.0, 1.0, 1.0])
        for w_01, w_02 in zip(knots_w_01, knots_w_02):
            w_02_mapped = (w_02 - 1.0) / (3.0 - 1.0)
            self.assertAlmostEqual(w_01, w_02_mapped)

        ## Check control points.
        for node_1, node_2 in zip(geometry_01, geometry_02):
            self.assertAlmostEqual(node_1.X, node_2.X)
            self.assertAlmostEqual(node_1.Y, node_2.Y)
            self.assertAlmostEqual(node_1.Z, node_2.Z)

        # Check if geoemtry is similar.
        param = KratosMultiphysics.Vector(3)
        param[0] = 0.0
        param[1] = 0.0
        param[2] = 0.0
        param2 = KratosMultiphysics.Vector(3)
        for _ in range(11):
            for _ in range(11):
                for _ in range(11):
                    coord_01 = geometry_01.GlobalCoordinates(param)
                    param2[0] = param[0]*(5.0 + 2.22) - 2.22
                    param2[1] = param[1]*(5.5 + 5.5) - 5.5
                    param2[2] = param[2]*(3.0 - 1.0) + 1.0
                    coord_02 = geometry_02.GlobalCoordinates(param2)
                    param[2] += 0.1
                    self.assertVectorAlmostEqual(coord_01, coord_02)
                param[1] += 0.1
                param[2] = 0.0
            param[0] += 0.1
            param[1] = 0.0




if __name__ == '__main__':
    KratosUnittest.main()
