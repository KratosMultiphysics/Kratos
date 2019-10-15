import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from os.path import join

import numpy as np

class TestCoSimulationInterface(KratosUnittest.TestCase):
    def test_cosimulation_interface(self):
        parameter_file_name = "test_parameters.json"
        cs_data_structure = ImportDataStructure(parameter_file_name)

        # get reference to global Variables
        pressure = vars(KM)['PRESSURE']
        temperature = vars(KM)['TEMPERATURE']
        density = vars(KM)['DENSITY']
        force = vars(KM)['FORCE']
        displacement = vars(KM)['DISPLACEMENT']

        # read in Parameters
        parameter_file_name = "pykratos/test_cosimulation_interface.json"
        with open(parameter_file_name, 'r') as parameter_file:
            par = cs_data_structure.Parameters(parameter_file.read())

        # create Model, ModelParts and Nodes from Parameters
        model = cs_data_structure.Model()
        model_parts = {}

        mp_keys = []
        for par_interface in [par['interface_a'], par['interface_b']]:
            mp_keys += par_interface.keys()
        mp_keys = list(set(mp_keys))

        counter = 0
        node_id = 0
        for key in mp_keys:
            model_parts[key] = model.CreateModelPart(key)
            for par_interface in [par['interface_a'], par['interface_b']]:
                if key in par_interface.keys():
                    for var in par_interface[key].list():
                        model_parts[key].AddNodalSolutionStepVariable(vars(KM)[var.GetString()])
            for i in range(2 + counter):
                model_parts[key].CreateNewNode(str(node_id), float(i), float(counter), 0.)
                node_id += 1
            counter += 1

        # create CoSimulationInterfaces from Parameters
        interface_a = CoSimulationInterface(model, par['interface_a'])
        interface_b = CoSimulationInterface(model, par['interface_b'])

        # initialize some values in Nodes ***
        for node in model_parts['mp_1'].Nodes:
            node.SetSolutionStepValue(temperature, 0, 1.)
            node.SetSolutionStepValue(force, 0, [2., 2., 2.])

        for node in model_parts['mp_2'].Nodes:
            node.SetSolutionStepValue(density, 0, 3.)

        for node in model_parts['mp_3'].Nodes:
            node.SetSolutionStepValue(displacement, 0, [4., 4., 4.])

        # test GetNumpyArray and SetNumpyArray
        for interface in [interface_a, interface_b]:
            data = interface.GetNumpyArray()
            input_data = np.arange(data.size).reshape(data.shape) * 1.
            interface.SetNumpyArray(input_data)
            output_data = interface.GetNumpyArray()
            self.assertListEqual(input_data.tolist(), output_data.tolist())

        # test __add__ *** TO DO
        data = interface_a.GetNumpyArray()
        data_1 = np.random.randint(0, 10, data.shape)
        data_2 = np.random.randint(0, 10, data.shape)
        interface_1 = interface_a.deepcopy()
        interface_2 = interface_a.deepcopy()
        interface_1.SetNumpyArray(data_1)
        interface_2.SetNumpyArray(data_2)

        data_3 = data_1 + data_2
        interface_3 = interface_1 + interface_2
        self.assertListEqual(data_3.tolist(), interface_3.GetPythonList())

        data_3 = data_1 + 5
        interface_3 = interface_1 + 5
        self.assertListEqual(data_3.tolist(), interface_3.GetPythonList())
        interface_3 = 5 + interface_1
        self.assertListEqual(data_3.tolist(), interface_3.GetPythonList())

        # test __sub__
        data_3 = data_1 - data_2
        interface_3 = interface_1 - interface_2
        self.assertListEqual(data_3.tolist(), interface_3.GetPythonList())

        data_3 = data_1 - 5
        interface_3 = interface_1 - 5
        self.assertListEqual(data_3.tolist(), interface_3.GetPythonList())

        data_3 = 5 - data_1
        interface_3 = 5 - interface_1
        self.assertListEqual(data_3.tolist(), interface_3.GetPythonList())

        # test __mul__
        data_3 = data_1 * data_2
        interface_3 = interface_1 * interface_2
        self.assertListEqual(data_3.tolist(), interface_3.GetPythonList())

        data_3 = data_1 * 5
        interface_3 = interface_1 * 5
        self.assertListEqual(data_3.tolist(), interface_3.GetPythonList())
        interface_3 = 5 * interface_1
        self.assertListEqual(data_3.tolist(), interface_3.GetPythonList())

        # test __truediv__
        data_3 = data_1 / 5
        interface_3 = interface_1 / 5
        self.assertListEqual(data_3.tolist(), interface_3.GetPythonList())

        # test __pow__
        data_3 = data_1 ** 5
        interface_3 = interface_1 ** 5
        self.assertListEqual(data_3.tolist(), interface_3.GetPythonList())

        # print output of test
        if False:
            print('\nTestCoSimulationInterface successful.\n')
            self.assertTrue(False)


if __name__ == '__main__':
    KratosUnittest.main()