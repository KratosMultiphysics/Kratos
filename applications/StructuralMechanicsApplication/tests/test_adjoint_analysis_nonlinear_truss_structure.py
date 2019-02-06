from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import os
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_mechanics_analysis
import structural_response_function_factory
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import KratosMultiphysics.kratos_utilities as kratos_utils
import time as timer

class AdjointSensitivityNonlinearTruss(KratosUnittest.TestCase):

    def _removeH5Files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                kratos_utils.DeleteFileIfExisting(name)

    def testThreeDimensionalStructure(self):
        with open("response_function_tests/adjoint_strain_energy_response_parameters_3D_truss.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters( parameter_file.read())

        model = KratosMultiphysics.Model()
        response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["kratos_response_settings"], model)
        model_part_adjoint = response_function.adjoint_model_part
        response_function.RunCalculation(calculate_gradient=True)
        
        node_sensitivity = []
        for i in range(1,6):
            node_sensitivity.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY) )
            print("node_sensitivity[i]", node_sensitivity[i - 1])

        ## Testing Adjoint sensitivity values with Finite difference results
        self.assertAlmostEqual(node_sensitivity[0][0], 0.003973418927216699, 2)
        self.assertAlmostEqual(node_sensitivity[0][1], 0.00012844054708693875, 2)
        self.assertAlmostEqual(node_sensitivity[0][2], 12.575784792545617, 1)
        self.assertAlmostEqual(node_sensitivity[1][0], -6.742208012155969, 1)
        self.assertAlmostEqual(node_sensitivity[1][1], -3.371206929614345, 1)
        self.assertAlmostEqual(node_sensitivity[1][2], -3.143893783885687, 1)
        self.assertAlmostEqual(node_sensitivity[2][0], 6.742758621669508, 1)
        self.assertAlmostEqual(node_sensitivity[2][1], -3.3712069225089176, 1)
        self.assertAlmostEqual(node_sensitivity[2][2], -3.143893783885687, 1)

        self._removeH5Files("primal_output_truss")

    def testMisesTrussStructure(self):
        with open("response_function_tests/adjoint_strain_energy_response_parameters_truss.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters( parameter_file.read())

        model = KratosMultiphysics.Model()
        response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["kratos_response_settings"], model)
        model_part_adjoint = response_function.adjoint_model_part
        response_function.RunCalculation(calculate_gradient=True)
        
        node_sensitivity = []
        for i in range(1,4):
            node_sensitivity.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY) )
            print("node_sensitivity[i]", node_sensitivity[i - 1])

        # ## Testing Adjoint sensitivity values with Finite difference results
        self.assertAlmostEqual(node_sensitivity[0][0],0.196686219489095, 0)
        self.assertAlmostEqual(node_sensitivity[0][1], -5.9577956088574515, 1)
        self.assertAlmostEqual(node_sensitivity[0][2], 0.0002899730233707487, 1)
        self.assertAlmostEqual(node_sensitivity[1][0], -14.449074176425823, 1)
        self.assertAlmostEqual(node_sensitivity[1][1], 3.0951781273103047, 0)
        self.assertAlmostEqual(node_sensitivity[1][2], 0.00014498642286753238, 1)
        self._removeH5Files("primal_output_truss")

if __name__ == '__main__':
    KratosUnittest.main()