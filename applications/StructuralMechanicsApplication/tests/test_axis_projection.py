from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

class TestAxisProjection(KratosUnittest.TestCase):

    def _set_up_test_case(self,mp,element_name):
        #create cube
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,2.0,0.0,0.0)
        mp.CreateNewNode(3,2.0,2.0,0.0)
        mp.CreateNewNode(4,0.0,2.0,0.0)

        mp.CreateNewNode(5,0.0,0.0,3.0)
        mp.CreateNewNode(6,2.0,0.0,3.0)
        mp.CreateNewNode(7,2.0,2.0,3.0)
        mp.CreateNewNode(8,0.0,2.0,3.0)

        mp.CreateNewElement(element_name, 1, [1,2,6,5], mp.GetProperties()[0])
        mp.CreateNewElement(element_name, 2, [2,3,7,6], mp.GetProperties()[0])
        mp.CreateNewElement(element_name, 3, [4,8,7,3], mp.GetProperties()[0])
        mp.CreateNewElement(element_name, 4, [4,1,5,8], mp.GetProperties()[0])



    def test_RadialProjection_Shell(self):
        self.execute_test_radial_projection("ShellThinElementCorotational3D4N")

    def test_PlanarProjection_Shell(self):
        self.execute_test_planar_projection("ShellThinElementCorotational3D4N")

    def test_RadialProjection_Membrane(self):
        self.execute_test_radial_projection("MembraneElement3D4N")

    def test_PlanarProjection_Membrane(self):
        self.execute_test_planar_projection("MembraneElement3D4N")

    def execute_test_radial_projection(self, element_name):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("solid_part")

        self._set_up_test_case(model_part,element_name)
        projection_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"  : "Structure",
            "echo_level"       : 1,
            "projection_type"  : "radial",
            "global_direction" : [0,0,1],
            "variable_name"    : "LOCAL_PRESTRESS_AXIS_1",
            "method_specific_settings" : { }
        }
        """)
        StructuralMechanicsApplication.ProjectVectorOnSurfaceUtility.Execute(model_part, projection_settings)

        expected_results = [[1.0,0.0,0.0],[0.0,1.0,0.0],[-1.0,0.0,0.0],[0.0,-1.0,0.0]]


        for i,element_i in enumerate(model_part.Elements):
            result_i = element_i.GetValue(StructuralMechanicsApplication.LOCAL_PRESTRESS_AXIS_1)
            for j in range(3): self.assertAlmostEqual(result_i[j],expected_results[i][j])

    def execute_test_planar_projection(self, element_name):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("solid_part")
        self._set_up_test_case(model_part, element_name)
        projection_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"  : "Structure",
            "echo_level"       : 1,
            "projection_type"  : "planar",
            "global_direction" : [1,0,1],
            "variable_name"    : "LOCAL_PRESTRESS_AXIS_1",
            "method_specific_settings" : { }
        }
        """)
        StructuralMechanicsApplication.ProjectVectorOnSurfaceUtility.Execute(model_part, projection_settings)

        res_scalar = 1.0/math.sqrt(2.0)
        expected_results = [[res_scalar,0.0,res_scalar],[0.0,0.0,1.0],[res_scalar,0.0,res_scalar],[0.0,0.0,1.0]]

        for i,element_i in enumerate(model_part.Elements):
            result_i = element_i.GetValue(StructuralMechanicsApplication.LOCAL_PRESTRESS_AXIS_1)
            for j in range(3): self.assertAlmostEqual(result_i[j],expected_results[i][j])


if __name__ == '__main__':
    KratosUnittest.main()
