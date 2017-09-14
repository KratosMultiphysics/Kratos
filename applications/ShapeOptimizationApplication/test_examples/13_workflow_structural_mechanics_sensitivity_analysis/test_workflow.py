import os
os.environ['OMP_NUM_THREADS'] = "6"
from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_MainKratos

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class TestCase(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def _remove_file(self, file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)

    def _create_test(self, parameter_file_name):
        with open(parameter_file_name + '_parameters.json', 'r') as parameter_file:
            project_parameters = Parameters(parameter_file.read())
            parameter_file.close()
        test = test_MainKratos.MainKratos(project_parameters)
        return test

    def solve(self, parameter_file_name):
        test = self._create_test(parameter_file_name)
        test.Solve()

    def test_Beam(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve structure
            self.solve('beam_test')
            #self.solve('test_beam_sensitivities/beam_test')
            # solve adjoint
            print("Primal solution process finished")
            test = self._create_test('beam_test_adjoint')
            test.Solve()
            print("Adjoint solution process finished and sensitivities are computed")    
           # Sensitivity = [[]]
           # Sensitivity[0].append(test.main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_X))
           # Sensitivity[0].append(test.main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))
    def circle_plate(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve structure
            self.solve('circle_plate_test')
            #self.solve('test_beam_sensitivities/beam_test')
            # solve adjoint
            print("Primal solution process finished!")
            test = self._create_test('circle_plate_test_adjoint')
            test.Solve()
            print("Adjoint solution process finished and sensitivities are computed")    
           # Sensitivity = [[]]
           # Sensitivity[0].append(test.main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_X))
           # Sensitivity[0].append(test.main_model_part.GetNode(1968).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))       

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
