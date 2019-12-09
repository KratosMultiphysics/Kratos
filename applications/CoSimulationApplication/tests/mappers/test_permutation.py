import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
import os

import time
from contextlib import contextmanager
@contextmanager
def timer(name=None, t=0, n=0, ms=False):
    startTime = time.time()
    yield
    elapsedTime = time.time() - startTime
    if ms:
        s = '\n' * n + '\t' * t + f'{elapsedTime * 1000:.2f}ms'
        s.replace(',', ' ')
    else:
        s = '\n' * n + '\t' * t + f'{elapsedTime:.1f}s'
    if name is not None:
        s += f' - {name}'
    s += '\n' * n
    print(s)

class TestMapperPermutation(KratosUnittest.TestCase):
    def test_mapper_permutation(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_permutation.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())


        """
        What to test??
        --------------
        
        - initialize from both sides, check output of both methods
        - 
        
        """

        var_from = vars(KM)["TEMPERATURE"]
        model_from = cs_data_structure.Model()
        model_part_from = model_from.CreateModelPart('wall_from')
        model_part_from.AddNodalSolutionStepVariable(var_from)

        for i in range(11):
            model_part_from.CreateNewNode(i, i, 2 * i, 3 * i)

        mapper = cs_tools.CreateInstance(parameters['mapper'])

        model_part_to = mapper.Initialize(model_part_from, None)


        print(model_part_to)


        # mapper((model_part_from, var_from), (model_part_to, var_to))

        # var_to = vars(KM)["PRESSURE"]
        # model_to = cs_data_structure.Model()
        # model_part_to = model_to.CreateModelPart('wall_to')
        # model_part_to.AddNodalSolutionStepVariable(var_to)
        #
        # for i in range(4):
        #     model_part_to.CreateNewNode(i, 0.0, 0.0, i / 3)
        #
        #
        #





        # mapper.Initialize(None, model_part_to)



        #
        # values = []
        # for i_to, node in enumerate(model_part_to.Nodes):
        #     values.append(node.GetSolutionStepValue(var_to))
        #
        # self.assertListEqual(values, [0., 9., 49., 100.])


if __name__ == '__main__':
    KratosUnittest.main()
