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
        with open("./adjoint_sensitivity_test_martin/adjoint_beam_structure/beam_test_parameters.json",'r') as parameter_file:
            ProjectParametersPrimal = Parameters( parameter_file.read())
        parameter_file.close()
        primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal)

        primal_analysis.Run()

    def test_local_stress_response(self):
        # Solve primal problem (only in one test case necessary)
        self._solve_primal_problem()

        # Create the adjoint solver
        with open("./adjoint_sensitivity_test_martin/adjoint_beam_structure/beam_test_local_stress_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())
        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)
        adjoint_analysis.Run()

        #out1 = self.adjoint_analysis.GetModelPart().Elements[1].CalculateOnIntegrationPoints(StructuralMechanicsApplication.I22_SENSITIVITY,self.adjoint_analysis.GetModelPart().ProcessInfo)

    def test_nodal_displacement_response(self):
        # Create the adjoint solver
        with open("./adjoint_sensitivity_test_martin/adjoint_beam_structure/beam_test_nodal_disp_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())
        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)

        adjoint_analysis.Run()

    def test_train_energy_response(self):
        # Create the adjoint solver
        with open("./adjoint_sensitivity_test_martin/adjoint_beam_structure/beam_test_strain_energy_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())
        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)

        adjoint_analysis.Run()

        # Delete *.h5 only after last test case because primal solution is used in each test case
        self._remove_file("./adjoint_sensitivity_test_martin/adjoint_beam_structure/beam_test_0.h5")
        self._remove_file("./adjoint_sensitivity_test_martin/adjoint_beam_structure/Beam_structure.time")

    def tearDown(self):
        pass
      

if __name__ == '__main__':
    KratosUnittest.main()
