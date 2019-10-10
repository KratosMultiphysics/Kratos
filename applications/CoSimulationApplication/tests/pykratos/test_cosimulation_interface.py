import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from os.path import join

class TestCoSimulationInterface(KratosUnittest.TestCase):
    def test_cosimulation_interface(self):
        parameter_file_name = "test_parameters.json"
        cs_data_structure = ImportDataStructure(parameter_file_name)

        # *** TODO

        # print output of test
        if False:
            print('\nTestCoSimulationInterface successful.\n')
            self.assertTrue(False)


if __name__ == '__main__':
    KratosUnittest.main()