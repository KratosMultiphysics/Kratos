from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

class TestAxisProjection(KratosUnittest.TestCase):

    def _set_up_test_case(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        #create cube
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,2.0,0.0,0.0)
        mp.CreateNewNode(3,2.0,2.0,0.0)
        mp.CreateNewNode(4,0.0,2.0,0.0)

        mp.CreateNewNode(5,0.0,0.0,3.0)
        mp.CreateNewNode(6,2.0,0.0,3.0)
        mp.CreateNewNode(7,2.0,2.0,3.0)
        mp.CreateNewNode(8,0.0,2.0,3.0)

        element_name = "ShellThinElementCorotational3D4N"
        mp.CreateNewElement(element_name, 1, [1,5,6,2], mp.GetProperties()[0])
        mp.CreateNewElement(element_name, 2, [2,3,7,6], mp.GetProperties()[0])
        mp.CreateNewElement(element_name, 3, [4,3,7,8], mp.GetProperties()[0])
        mp.CreateNewElement(element_name, 4, [4,1,5,8], mp.GetProperties()[0])

        return mp



    def test_PlanarProjection(self):

        model_part = self._set_up_test_case()
        projection_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"  : "Structure",
            "echo_level"       : 1,
            "projection_type"  : "planar",
            "global_direction" : [0,0,1],
            "variable_name"    : "LOCAL_PRESTRESS_AXIS_1",
            "method_specific_settings" : { }
        }
        """)
        StructuralMechanicsApplication.ProjectVectorOnSurfaceUtility.Execute(model_part, projection_settings)

        for element_i in model_part.Elements:
            print(element_i.GetValue(StructuralMechanicsApplication.LOCAL_PRESTRESS_AXIS_1))



if __name__ == '__main__':
    KratosUnittest.main()
