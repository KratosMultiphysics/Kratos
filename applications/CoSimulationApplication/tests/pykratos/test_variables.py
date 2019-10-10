import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

from os.path import join
import copy

class TestPyKratosVariables(KratosUnittest.TestCase):
    def test_pykratos_variables(self):
        parameter_file_name = "test_parameters.json"
        cs_data_structure = ImportDataStructure(parameter_file_name)

        pres = "PRESSURE"
        disp = "DISPLACEMENT"

        # calling global variables (different ways)
        var_pres = vars(KM)[pres]
        var_pres_1 = KM.__dict__[pres]
        var_pres_2 = KM.KratosGlobals.GetVariable(pres)
        var_pres_3 = copy.deepcopy(var_pres)
        print(id(var_pres), id(var_pres_1), id(var_pres_2), id(var_pres_3))

        self.assertEqual(var_pres, var_pres_1)
        self.assertEqual(var_pres, var_pres_2)
        self.assertNotEqual(var_pres, var_pres_3)

        var_disp = vars(KM)[disp]
        var_disp_1 = KM.__dict__[disp]
        var_disp_2 = KM.KratosGlobals.GetVariable(disp)
        var_disp_3 = copy.deepcopy(var_disp)
        print(id(var_disp), id(var_disp_1), id(var_disp_2), id(var_disp_3))

        self.assertEqual(var_disp, var_disp_1)
        self.assertEqual(var_disp, var_disp_2)
        self.assertNotEqual(var_disp, var_disp_3)

        # dynamically add global variables
        """
        I don't think this is possible.
        Variables are defined as globals in the module.
        But it is not possible to add stuff to the module, 
        from outside the module, so no globals can be added.
        """

        # different types of variables? *** --> check in test_cosimulation_interface
        print(var_pres)
        print(var_disp)

        print(type(var_pres))
        print(type(var_disp))


        # print output of test
        if False:
            print('\nTestPyKratosVariables successful.\n')
            self.assertTrue(False)


if __name__ == '__main__':
    KratosUnittest.main()