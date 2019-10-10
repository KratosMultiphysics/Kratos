import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from os.path import join

class TestPyKratosParameters(KratosUnittest.TestCase):
    def test_pykratos_parameters(self):
        parameter_file_name = "test_parameters.json"
        cs_data_structure = ImportDataStructure(parameter_file_name)

        file_name = join('pykratos', 'test_parameters_class.json')
        with open(file_name, 'r') as file:
            par = cs_data_structure.Parameters(file.read())



        # test getters and setters
        par.SetInt('test_int', 2)
        par.SetDouble('test_float', 2.2)
        par.SetBool('test_bool', False)
        par.SetString('test_str', 'two')

        self.assertEqual(par['test_int'].GetInt(), 2)
        self.assertEqual(par['test_float'].GetDouble(), 2.2)
        self.assertEqual(par['test_bool'].GetBool(), False)
        self.assertEqual(par['test_str'].GetString(), 'two')

        # *** what to test??




        # *** how do lists/arrays work wrt Parameters class?


        print('\n\nREACHED END OF TEST')
        self.assertTrue(False)
        # *** stdout is only shown when test encounters an error


if __name__ == '__main__':
    KratosUnittest.main()