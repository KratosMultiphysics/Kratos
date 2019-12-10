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

class TestMapperCombined(KratosUnittest.TestCase):
    def test_mapper_combined(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_combined.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())

        # compare 3 mappers: nearest, combined_a, combined_b
        if True:
            # create model_part_from
            var_from = vars(KM)["TEMPERATURE"]
            model_from = cs_data_structure.Model()
            model_part_from = model_from.CreateModelPart('wall_from')
            model_part_from.AddNodalSolutionStepVariable(var_from)

            for i in range(100):
                node = model_part_from.CreateNewNode(i, np.random.rand(), np.random.rand(), np.random.rand())
                node.SetSolutionStepValue(var_from, 0, np.random.rand())

            # create model_part_to
            var_to = vars(KM)["PRESSURE"]
            model_to = cs_data_structure.Model()
            model_part_to = model_to.CreateModelPart('wall_to')
            model_part_to.AddNodalSolutionStepVariable(var_to)

            for i in range(100):
                model_part_to.CreateNewNode(i, np.random.rand(), np.random.rand(), np.random.rand())

            # create mappers, get output data
            data = []
            for mapper_name in ['mapper_nearest', 'mapper_combined_a', 'mapper_combined_b']:
                mapper = cs_tools.CreateInstance(parameters[mapper_name])
                mapper.Initialize(model_part_from, model_part_to)
                mapper((model_part_from, var_from), (model_part_to, var_to))

                values = []
                for i_to, node in enumerate(model_part_to.Nodes):
                    values.append(node.GetSolutionStepValue(var_to))
                data.append(values)

            # check output data
            self.assertListEqual(data[0], data[1])
            self.assertListEqual(data[0], data[2])


if __name__ == '__main__':
    KratosUnittest.main()
