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

        # check if method Initialize works
        if True:
            var = vars(KM)["TEMPERATURE"]
            model = cs_data_structure.Model()
            model_part_from = model.CreateModelPart('wall_from')
            model_part_from.AddNodalSolutionStepVariable(var)

            model_part_from.CreateNewNode(0, 0., 1., 2.)

            # model_part_from given
            mapper = cs_tools.CreateInstance(parameters['mapper'])
            model_part_to = mapper.Initialize(model_part_from, forward=True)
            node = model_part_to.Nodes[0]
            self.assertListEqual([node.X, node.Y, node.Z], [2., 0., 1.])

            # model_part_to given
            mapper = cs_tools.CreateInstance(parameters['mapper'])
            model_part_from = mapper.Initialize(model_part_to, forward=False)
            node = model_part_from.Nodes[0]
            self.assertListEqual([node.X, node.Y, node.Z], [0., 1., 2.])

        # check if method __call__ works
        if True:
            var = vars(KM)["TEMPERATURE"]
            model = cs_data_structure.Model()
            model_part_from = model.CreateModelPart('wall_from')
            model_part_from.AddNodalSolutionStepVariable(var)

            for i in range(10):
                node = model_part_from.CreateNewNode(i, i * 1., i * 2., i * 3.)
                node.SetSolutionStepValue(var, 0, np.random.rand())

            mapper = cs_tools.CreateInstance(parameters['mapper'])
            model_part_to = mapper.Initialize(model_part_from, forward=True)
            mapper((model_part_from, var), (model_part_to, var))

            for node_from, node_to in zip(model_part_from.Nodes, model_part_to.Nodes):
                val_from = node_from.GetSolutionStepValue(var)
                val_to = node_to.GetSolutionStepValue(var)
                self.assertEqual(val_from, val_to)


if __name__ == '__main__':
    KratosUnittest.main()
