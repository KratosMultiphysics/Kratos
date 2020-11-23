# Importing the Kratos Library
import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.modeler_factory import KratosModelerFactory

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def run_modeler(current_model, modelers_list):
    factory = KratosModelerFactory()
    list_of_modelers = factory.ConstructListOfModelers(current_model, modelers_list)
    for modeler in list_of_modelers:
        modeler.SetupGeometryModel()

    for modeler in list_of_modelers:
        modeler.SetupModelPart()

class TestNurbsVolumeGridModeler(KratosUnittest.TestCase):

    def test_grid_modeler(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsVolumeGridModeler",
            "Parameters": {
                "model_part_name" : "Mesh",
                "lower_point": [0.0,0.0,0.0],
                "upper_point": [1.0,1.0,1.0],
                "polynomial_order" : [4, 1, 1],
                "number_of_knot_spans" : [3,1,1]
            }
        }]
        """)
        run_modeler(current_model, modeler_settings)
        model_part = current_model.GetModelPart("Mesh")
        # Check control points
        nodes = model_part.NodesArray(0)
        for i in range(4):
            self.assertAlmostEqual(nodes[i*7+0].X,0.0)
            self.assertAlmostEqual(nodes[i*7+1].X,0.083333333334)
            self.assertAlmostEqual(nodes[i*7+2].X,0.25)
            self.assertAlmostEqual(nodes[i*7+3].X,0.5)
            self.assertAlmostEqual(nodes[i*7+4].X,0.75)
            self.assertAlmostEqual(nodes[i*7+5].X,0.916666666667)
            self.assertAlmostEqual(nodes[i*7+6].X,1.0)
        for i in range(7):
            self.assertAlmostEqual(nodes[i].Y,0.0)
            self.assertAlmostEqual(nodes[i].Z,0.0)
        for i in range(7,14):
            self.assertAlmostEqual(nodes[i].Y,1.0)
            self.assertAlmostEqual(nodes[i].Z,0.0)
        for i in range(14,21):
            self.assertAlmostEqual(nodes[i].Y,0.0)
            self.assertAlmostEqual(nodes[i].Z,1.0)
        for i in range(21,28):
            self.assertAlmostEqual(nodes[i].Y,1.0)
            self.assertAlmostEqual(nodes[i].Z,1.0)

        geometry = model_part.GetGeometry(1)
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
    
    def test_grid_modeler_02(self):
        current_model = KratosMultiphysics.Model()
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsVolumeGridModeler",
            "Parameters": {
                "model_part_name" : "Mesh_01",
                "lower_point": [-1.0,-0.5,1.0],
                "upper_point": [1.0,0.5,3.0],
                "polynomial_order" : [4, 5, 2],
                "number_of_knot_spans" : [3,6,7]
            } }, {
            "modeler_name": "NurbsVolumeGridModeler",
            "Parameters": {
                "model_part_name" : "Mesh_02",
                "lower_point": [-1.0,-0.5,1.0],
                "upper_point": [1.0,0.5,3.0],
                "polynomial_order" : [1, 1, 1],
                "number_of_knot_spans" : [1,1,1]
            }
        }]
        """)
        
        run_modeler(current_model, modeler_settings)
        model_part_01 = current_model.GetModelPart("Mesh_01")
        model_part_02 = current_model.GetModelPart("Mesh_02")

        geometry_01 = model_part_01.GetGeometry(1)
        geometry_02 = model_part_02.GetGeometry(1)

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
        # geometry 01
        self.assertEqual(geometry_01.PolynomialDegreeU(), 4)        
        self.assertEqual(geometry_01.PolynomialDegreeV(), 5)        
        self.assertEqual(geometry_01.PolynomialDegreeW(), 2)
        # geometry 02
        self.assertEqual(geometry_02.PolynomialDegreeU(), 1)        
        self.assertEqual(geometry_02.PolynomialDegreeV(), 1)        
        self.assertEqual(geometry_02.PolynomialDegreeW(), 1)

        ## Check number of nodes
        # geomtry 01
        self.assertEqual(len(geometry_01),693)
        self.assertEqual(model_part_01.NumberOfNodes(), 693)
        # geomtry 02
        self.assertEqual(len(geometry_02),8)
        self.assertEqual(model_part_02.NumberOfNodes(), 8)  

        #Check if geoemtry is similar.
        param = KratosMultiphysics.Vector(3)
        param[0] = 0.0
        param[1] = 1.0
        param[2] = 0.0
        for i in range(11):
            for j in range(11):
                print(param)
                coord_01 = geometry_01.GlobalCoordinates(param)
                coord_02 = geometry_02.GlobalCoordinates(param)
                param[2] += 0.1
                print(coord_01)
                print(coord_02)
                self.assertVectorAlmostEqual(coord_01, coord_02)
            param[0] += 0.1
            param[1] -= 0.1
            param[2] = 0.0

if __name__ == "__main__":
    KratosUnittest.main()