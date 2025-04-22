import KratosMultiphysics as kratos
from KratosMultiphysics import KratosUnittest

# Import os module
import os

# Import math module
import math

def CalculateAnalyticalDistance(node):
    """Distance from a node to the surface of a  sphere of radius 0.5 centered at the origin."""
    distance = 0.5 - math.sqrt(node.X**2 + node.Y**2 + node.Z**2)
    return distance

class TestSurrogateBoundaryModelerCoarseSphere(KratosUnittest.TestCase):

    @classmethod
    def setUpClass(self): 
        """Prepare the test environment with model parts and the parameters"""
        self.current_model = kratos.Model()
        self.model_part = self.current_model.CreateModelPart("main_model_part")
        self.skin_model_part = self.current_model.CreateModelPart("skin_model_part")

        self.settings = kratos.Parameters(""" 
        {
            "output_model_part_name" : "main_model_part",
            "input_model_part_name" : "skin_model_part",
            "mdpa_file_name" : "auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_skin",
            "key_plane_generator": {
                "Parameters" : {
                    "voxel_sizes" : [0.15, 0.15, 0.15],
                    "min_point" : [-0.7, -0.7, -0.7],
                    "max_point" : [0.7, 0.7, 0.7]
                }
            },
            "coloring_settings_list": [
            {
                "type" : "cells_in_touch",
                "model_part_name": "skin_model_part",
                "color": -1,
                "input_entities": "conditions"
            },
            {
                "type" : "cells_with_inside_center",
                "model_part_name": "skin_model_part",
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
        """)
    
    def setUp(self):
        """Setup tasks before each test method."""
        pass

    def test_ComputeDistanceToSkin(self):
        surrogate_boundary = kratos.SurrogateBoundaryModeler(self.current_model, self.settings)

        surrogate_boundary.SetupGeometryModel()
        surrogate_boundary.PrepareGeometryModel()
        surrogate_boundary.SetupModelPart()

        main_model_part = self.current_model["main_model_part"]

        surrogate_boundary.ComputeSurrogateBoundary()
        self.assertEqual(len(surrogate_boundary.GetSurrogateBoundaryNodes()), 1000)

        for node in surrogate_boundary.GetSurrogateBoundaryNodes():
            if node.IsActive():
                distance = CalculateAnalyticalDistance(node.GetNodePtr())
                signed_distance = node.GetSignedDistance()
                
                self.assertAlmostEqual(distance, signed_distance, delta=4e-2)



if __name__ == '__main__':
    # Configure logging level and start the test runner
    KratosUnittest.main()