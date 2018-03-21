from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
os.environ['OMP_NUM_THREADS'] = "1"
#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_mechanics_analysis

class TestCase(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def _remove_file(self, file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)

    def _solve_primal_problem(self):
        with open("./adjoint_sensitivity_test_martin/adjoint_shell_structure_3d4n/linear_shell_test_parameters.json",'r') as parameter_file:
            ProjectParametersPrimal = Parameters( parameter_file.read())
        parameter_file.close()
        primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal)

        primal_analysis.Run()

    def test_local_stress_response(self): 
        # Solve primal problem (only here necessary. The other tests corresponding to the same primal problem.)
        self._solve_primal_problem() 

        # Create the adjoint solver
        self._solve_primal_problem()
        with open("./adjoint_sensitivity_test_martin/adjoint_shell_structure_3d4n/linear_shell_test_local_stress_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())
        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)
        adjoint_analysis.Run()

        # Check sensitivities for the parameter THICKNESS
        reference_values = [0.41386635771251334, -2.601255910505967, 0.6628211584248311]
        sensitivities_to_check = []
        element_list = [1,5,8]
        for element_id in element_list:
            sensitivities_to_check.append(adjoint_analysis.GetModelPart().Elements[element_id].GetValue(KratosMultiphysics.StructuralMechanicsApplication.THICKNESS_SENSITIVITY))
    
        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

    def test_nodal_displacement_response(self):           
        # Create the adjoint solver
        with open("./adjoint_sensitivity_test_martin/adjoint_shell_structure_3d4n/linear_shell_test_nodal_disp_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())
        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)
        adjoint_analysis.Run()

        # Check sensitivities for the parameter THICKNESS
        reference_values = [-0.002801100145774953, -0.0027781711818140294, -0.00046158741678306966]
        sensitivities_to_check = []
        element_list = [1,5,8]
        for element_id in element_list:
            sensitivities_to_check.append(adjoint_analysis.GetModelPart().Elements[element_id].GetValue(KratosMultiphysics.StructuralMechanicsApplication.THICKNESS_SENSITIVITY))
    
        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

    def test_strain_energy_response(self):           
        # Create the adjoint solver
        with open("./adjoint_sensitivity_test_martin/adjoint_shell_structure_3d4n/linear_shell_test_strain_energy_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())
        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)
        adjoint_analysis.Run()

        # Check sensitivities for the parameter THICKNESS
        reference_values = [-0.003271375240263905, -0.005556342364194191, -0.001972358301203816]
        sensitivities_to_check = []
        element_list = [1,5,8]
        for element_id in element_list:
            sensitivities_to_check.append(adjoint_analysis.GetModelPart().Elements[element_id].GetValue(KratosMultiphysics.StructuralMechanicsApplication.THICKNESS_SENSITIVITY))
    
        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

        # Delete *.h5 only after last test case because primal solution is used in each test case
        self._remove_file("./adjoint_sensitivity_test_martin/adjoint_shell_structure_3d4n/rectangular_plate_0.h5")
        self._remove_file("./adjoint_sensitivity_test_martin/adjoint_shell_structure_3d4n/rectangular_plate.time")

    def tearDown(self):
        pass
    #TODO: add this test in test_StructualMechanicsApllication.py

if __name__ == '__main__':
    KratosUnittest.main()