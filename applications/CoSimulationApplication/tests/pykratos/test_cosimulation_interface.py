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

        # *** TODO

        # create variables
        pressure = vars(KM)['PRESSURE']  # DoubleVariable
        displacement = vars(KM)['DISPLACEMENT']  # Array3Variable

        # read in Parameters
        parameter_file_name = "pykratos/test_cosimulation_interface.json"
        with open(parameter_file_name, 'r') as parameter_file:
            par = cs_data_structure.Parameters(parameter_file.read())

        # create Model, ModelParts and Nodes
        model = cs_data_structure.Model()
        model_parts = {}
        counter = 0
        node_id = 0
        for key in [_.GetString() for _ in par['model_parts']]:
            model_parts[key] = model.CreateModelPart(key)
            model_parts[key].AddNodalSolutionStepVariable(pressure)
            model_parts[key].AddNodalSolutionStepVariable(displacement)

            for i in range(2 + counter):
                model_parts[key].CreateNewNode(str(node_id), float(i), float(counter), 0.)
                node_id += 1
            counter += 1

        # create CoSimulationInterfaces
        json_p = '{"mp_1": "PRESSURE", "mp_2": "PRESSURE", "mp_3": "PRESSURE"}'
        par_interface_p = cs_data_structure.Parameters(json_p)
        interface_p = CoSimulationInterface(model, par_interface_p)

        json_d ='{"mp_1": "DISPLACEMENT", "mp_2": "DISPLACEMENT", "mp_3": "DISPLACEMENT"}'
        par_interface_d = cs_data_structure.Parameters(json_d)
        interface_d = CoSimulationInterface(model, par_interface_d)

        # change data in interfaces
        p = interface_p.GetNumpyArray()
        interface_p.SetNumpyArray(np.arange(p.size) ** 2.)

        d = interface_d.GetNumpyArray()
        interface_d.SetNumpyArray(np.arange(d.size).reshape(d.shape) * 1.)

        # print information about nodes
        print('\nNew nodes:')
        for key in [_.GetString() for _ in par['model_parts']]:
            for node in model_parts[key].Nodes:
                print(f'Id = {node.Id}, id = {id(node)}, ' +
                      f'p = {node.GetSolutionStepValue(pressure)}, ' +
                      f'delta = {node.GetSolutionStepValue(displacement)}')

        # check ndarray operations *** TODO


        """
        why specify Variables in dict? 
        because we know in advance which Variable is used for
        input interface and which one for output interface
        and the modelparts should be the same on both interfaces
        (can't imagine when this is not the case)
        so you're just creating a possibililty to make mistakes...
        
        so: we should only specify the modelparts that are used
        on the interfaces, the rest is determined by the solver
        
        also: change CoSimulationInterface.__init__ so that it accepts
        a param object with a list of modelpart-names and
        only one variable specified??
        """

        # print output of test
        if False:
            print('\nTestCoSimulationInterface successful.\n')
            self.assertTrue(False)


if __name__ == '__main__':
    KratosUnittest.main()