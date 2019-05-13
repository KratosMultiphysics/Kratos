from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import os
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import structural_response_function_factory
import KratosMultiphysics.kratos_utilities as kratos_utils

class AdjointSensitivityNonlinearTruss(KratosUnittest.TestCase):
    def _removeH5AndRestFiles(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                kratos_utils.DeleteFileIfExisting(name)

    def testStructureWithFreeNodes(self):
        with open("response_function_tests/adjoint_strain_energy_response_parameters_truss_free_nodes.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters( parameter_file.read())

        model = KratosMultiphysics.Model()
        response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["kratos_response_settings"], model)
        model_part_adjoint = response_function.adjoint_model_part
        response_function.RunCalculation(calculate_gradient=True)

        node_primal_displacement = []
        node_adjoint_displacement = []
        node_sensitivity = []
        for i in range(1,4):
            node_sensitivity.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY) )
            node_adjoint_displacement.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT) )
            node_primal_displacement.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT) )


        ## Testing Primal Displacement
        self.assertAlmostEqual(node_primal_displacement[0][1], -0.3800988288518232, 5)
        self.assertAlmostEqual(node_primal_displacement[1][0], -0.01835460117547531, 5)
        self.assertAlmostEqual(node_primal_displacement[1][1], -0.152655188858641, 5)

        ## Testing Adjoint Displacement
        self.assertAlmostEqual(node_adjoint_displacement[0][1], -0.5919920188465935, 5)
        self.assertAlmostEqual(node_adjoint_displacement[1][0], -0.03786923343778928, 5)
        self.assertAlmostEqual(node_adjoint_displacement[1][1], -0.18677491827383022, 5)

        # ## Testing Adjoint sensitivity values with Finite difference results
        self.assertAlmostEqual(node_sensitivity[0][0], 0.01611482787211571, 0)
        self.assertAlmostEqual(node_sensitivity[0][1], -6.537811436491125, 1)
        self.assertAlmostEqual(node_sensitivity[1][0], -2.915941523706777, 1)
        self.assertAlmostEqual(node_sensitivity[1][1], 2.667856723981288, 1)
        self.assertAlmostEqual(node_sensitivity[2][0], -1.5577322715287776, 1)
        self.assertAlmostEqual(node_sensitivity[2][1], 0.6202639069030624, 1)

        ## Testing the Implementation
        self.assertAlmostEqual(node_sensitivity[0][0], 0.007062394882514367, 3)
        self.assertAlmostEqual(node_sensitivity[0][1], -6.529229670502979, 3)
        self.assertAlmostEqual(node_sensitivity[1][0], -2.8968539820576353, 3)
        self.assertAlmostEqual(node_sensitivity[1][1], 2.644373654128033, 3)
        self.assertAlmostEqual(node_sensitivity[2][0], -1.5669897149417131, 3)
        self.assertAlmostEqual(node_sensitivity[2][1], 0.626750518458256, 3)

        self._removeH5AndRestFiles("primal_output_truss")
        self._removeH5AndRestFiles("test_restart_file")

    def testStructureWithMultipleLoads(self):
        with open("response_function_tests/adjoint_strain_energy_response_parameters_truss_multiple_loads.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters( parameter_file.read())

        model = KratosMultiphysics.Model()
        response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["kratos_response_settings"], model)
        model_part_adjoint = response_function.adjoint_model_part
        response_function.RunCalculation(calculate_gradient=True)

        node_primal_displacement = []
        node_adjoint_displacement = []
        node_sensitivity = []
        for i in range(1,4):
            node_sensitivity.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY) )
            node_adjoint_displacement.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT) )
            node_primal_displacement.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT) )

        ## Testing Primal Displacement
        self.assertAlmostEqual(node_primal_displacement[0][1], -0.25026254024475564, 5)
        self.assertAlmostEqual(node_primal_displacement[1][0], 0.010468441473508287, 5)
        self.assertAlmostEqual(node_primal_displacement[1][1], -0.18251911698105952, 5)

        ## Testing Adjoint Displacement
        self.assertAlmostEqual(node_adjoint_displacement[0][1], -0.333990330226329, 5)
        self.assertAlmostEqual(node_adjoint_displacement[1][0], 0.019849352705151938, 5)
        self.assertAlmostEqual(node_adjoint_displacement[1][1], -0.2700670998184397, 4)

        # ## Testing Adjoint sensitivity values with Finite difference results
        self.assertAlmostEqual(node_sensitivity[0][0], 0.008305695109456224, 0)
        self.assertAlmostEqual(node_sensitivity[0][1], -0.056400936410128104, 3)
        self.assertAlmostEqual(node_sensitivity[1][0], 1.9167212722770441, 1)
        self.assertAlmostEqual(node_sensitivity[1][1], -1.9154809716148689, 1)
        self.assertAlmostEqual(node_sensitivity[2][0], -2.5470450155662405, 1)
        self.assertAlmostEqual(node_sensitivity[2][1], 1.956868369346054, 1)

        # Testing the Implementation
        self.assertAlmostEqual(node_sensitivity[0][0], 0.0019208582102663886, 2)
        self.assertAlmostEqual(node_sensitivity[0][1], -0.056662177222219556, 3)
        self.assertAlmostEqual(node_sensitivity[1][0],  1.924937643575508, 3)
        self.assertAlmostEqual(node_sensitivity[1][1],  -1.9293686239981187, 3)
        self.assertAlmostEqual(node_sensitivity[2][0],  -2.5522806761345747, 3)
        self.assertAlmostEqual(node_sensitivity[2][1], 1.9604805174396722, 3)

        self._removeH5AndRestFiles("primal_output_truss")
        self._removeH5AndRestFiles("test_restart_file")

    def testThreeDimensionalStructure(self):
        with open("response_function_tests/adjoint_strain_energy_response_parameters_3D_truss.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters( parameter_file.read())

        model = KratosMultiphysics.Model()
        response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["kratos_response_settings"], model)
        model_part_adjoint = response_function.adjoint_model_part
        response_function.RunCalculation(calculate_gradient=True)

        node_primal_displacement = []
        node_adjoint_displacement = []
        node_sensitivity = []
        for i in range(1,4):
            node_sensitivity.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY) )
            node_adjoint_displacement.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT) )
            node_primal_displacement.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT) )

        ## Testing Primal Displacement
        self.assertAlmostEqual(node_primal_displacement[0][2], -1.3113583781634623, 5)

        ## Testing Adjoint Displacement
        self.assertAlmostEqual(node_adjoint_displacement[0][2], -2.603827109923794, 5)

        ## Testing Adjoint sensitivity values with Finite difference results
        self.assertAlmostEqual(node_sensitivity[0][0], 0.004232732919717819, 2)
        self.assertAlmostEqual(node_sensitivity[0][1], 0.00012836380847147666, 2)
        self.assertAlmostEqual(node_sensitivity[0][2], 12.578750211389432, 2)
        self.assertAlmostEqual(node_sensitivity[1][0], -6.738579212139938, 1)
        self.assertAlmostEqual(node_sensitivity[1][1], -3.369399338737366, 1)
        self.assertAlmostEqual(node_sensitivity[1][2], -3.1446327909634415, 2)
        self.assertAlmostEqual(node_sensitivity[2][0], 6.739157221602453, 1)
        self.assertAlmostEqual(node_sensitivity[2][1], -3.3693993330530243, 1)
        self.assertAlmostEqual(node_sensitivity[2][2], -3.1446327938056124, 2)

        # Testing the Implementation
        self.assertAlmostEqual(node_sensitivity[0][0], 0.0002460614756534252, 3)
        self.assertAlmostEqual(node_sensitivity[0][1], 0.00047000270732744376, 3)
        self.assertAlmostEqual(node_sensitivity[0][2], 12.574351908036515, 3)
        self.assertAlmostEqual(node_sensitivity[1][0], -6.725441867648498, 3)
        self.assertAlmostEqual(node_sensitivity[1][1], -3.3626750259447986, 3)
        self.assertAlmostEqual(node_sensitivity[1][2], -3.1435332945081744, 3)
        self.assertAlmostEqual(node_sensitivity[2][0], 6.725718336704316, 3)
        self.assertAlmostEqual(node_sensitivity[2][1], -3.3626700765001027, 3)
        self.assertAlmostEqual(node_sensitivity[2][2], -3.1435037362864993, 3)

        self._removeH5AndRestFiles("primal_output_truss")
        self._removeH5AndRestFiles("test_restart_file")

    def testMisesTrussStructure(self):
        with open("response_function_tests/adjoint_strain_energy_response_parameters_truss.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters( parameter_file.read())

        model = KratosMultiphysics.Model()
        response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["kratos_response_settings"], model)
        model_part_adjoint = response_function.adjoint_model_part
        response_function.RunCalculation(calculate_gradient=True)

        node_primal_displacement = []
        node_adjoint_displacement = []
        node_sensitivity = []
        for i in range(1,4):
            node_sensitivity.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY) )
            node_adjoint_displacement.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT) )
            node_primal_displacement.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT) )

        ## Testing Primal Displacement
        self.assertAlmostEqual(node_primal_displacement[0][1], -0.28480876916354353, 5)

        ## Testing Adjoint Displacement
        self.assertAlmostEqual(node_adjoint_displacement[0][1], -0.6471369201631015, 4)

        # ## Testing Adjoint sensitivity values with Finite difference results
        self.assertAlmostEqual(node_sensitivity[0][1], -5.925022561292791, 1)
        self.assertAlmostEqual(node_sensitivity[1][0], -14.444793340473437, 1)
        self.assertAlmostEqual(node_sensitivity[1][1], 3.00866391906851, 0)

        ## Testing the Implementation
        self.assertAlmostEqual(node_sensitivity[0][1], -5.94525099835189, 3)
        self.assertAlmostEqual(node_sensitivity[1][0], -14.4823766392479, 3)
        self.assertAlmostEqual(node_sensitivity[1][1], 2.9728211538028146, 3)

        self._removeH5AndRestFiles("primal_output_truss")
        self._removeH5AndRestFiles("test_restart_file")

    def testThickShellStructure(self):
        with open("response_function_tests/adjoint_strain_energy_response_parameters_nonlinear_structure.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters( parameter_file.read())

        model = KratosMultiphysics.Model()
        response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["kratos_response_settings"], model)
        model_part_adjoint = response_function.adjoint_model_part
        response_function.RunCalculation(calculate_gradient=True)

        node_primal_displacement = []
        node_adjoint_displacement = []
        node_sensitivity = []
        for i in range(4,10):
            node_sensitivity.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY) )
            node_adjoint_displacement.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT) )
            node_primal_displacement.append( model_part_adjoint.GetNode(i).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT) )

        ## Testing Primal Displacement
        self.assertAlmostEqual(node_primal_displacement[5][1], 0.3594939823791553, 5)
        self.assertAlmostEqual(node_primal_displacement[5][2], 0.7873034436481292, 5)

        ## Testing Adjoint Displacement
        self.assertAlmostEqual(node_adjoint_displacement[5][1], 0.32234729813716384, 5)
        self.assertAlmostEqual(node_adjoint_displacement[5][2], 0.28211849486597773, 5)

        ## Testing the Implementation
        self.assertAlmostEqual(node_sensitivity[0][0], 0.11103597470787194, 5)
        self.assertAlmostEqual(node_sensitivity[0][1], 0.43866198490768604, 4)
        self.assertAlmostEqual(node_sensitivity[0][2], 0.27514421409649337, 4)
        self.assertAlmostEqual(node_sensitivity[1][2],  0.22171521858225285, 4)

        self.assertAlmostEqual(node_sensitivity[2][0], -0.09996784337216753, 5)
        self.assertAlmostEqual(node_sensitivity[2][1], 0.40379398660077026, 4)
        self.assertAlmostEqual(node_sensitivity[2][2], 0.2306230957482736, 4)

        self.assertAlmostEqual(node_sensitivity[3][1], -0.2956722449688707, 4)
        self.assertAlmostEqual(node_sensitivity[3][2], -0.3135344440155797, 4)

        self.assertAlmostEqual(node_sensitivity[4][1], -0.296582114991923, 4)
        self.assertAlmostEqual(node_sensitivity[4][2], -0.21369097199038822, 4)

        self.assertAlmostEqual(node_sensitivity[5][1], -0.31173182957715273, 4)
        self.assertAlmostEqual(node_sensitivity[5][2], -0.3433681880518615, 4)



        self._removeH5AndRestFiles("rectangular_plate_structure_primal")
        self._removeH5AndRestFiles("test_restart_file")

## Most of the time required for the test is used in reading and wrting out file
if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    nightSuite = suites['nightly']
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([AdjointSensitivityNonlinearTruss]))
    KratosUnittest.runTests(suites)
    allSuite = suites['all']
    allSuite.addTests(nightSuite)
    KratosUnittest.runTests(suites)