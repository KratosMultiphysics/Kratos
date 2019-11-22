import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
import time

class TestMapperNearest(KratosUnittest.TestCase):
    def test_mapper_nearest(self):
        self.assertAlmostEqual(1, 1)

        parameter_file_name = 'test_nearest.json'
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())

        # test 1D problem
        var_from = vars(KM)["TEMPERATURE"]
        model_from = cs_data_structure.Model()
        model_part_from = model_from.CreateModelPart('wall_from')
        model_part_from.AddNodalSolutionStepVariable(var_from)

        for i in range(11):
            node = model_part_from.CreateNewNode(i, 0.0, 0.0, 0.1 * i)
            node.SetSolutionStepValue(var_from, 0, i ** 2)

        var_to = vars(KM)["PRESSURE"]
        model_to = cs_data_structure.Model()
        model_part_to = model_to.CreateModelPart('wall_to')
        model_part_to.AddNodalSolutionStepVariable(var_to)

        for i in range(4):
            model_part_to.CreateNewNode(i, 0.0, 0.0, i / 3)

        mapper = cs_tools.CreateInstance(parameters['mapper'])
        mapper.Initialize(model_part_from, model_part_to)
        mapper((model_part_from, var_from), (model_part_to, var_to))

        values = []
        for i_to, node in enumerate(model_part_to.Nodes):
            values.append(node.GetSolutionStepValue(var_to))

        self.assertListEqual(values, [0., 9., 49., 100.])

        # test 1D problem with vectors?



        # test 3D problem with perturbed grid
        var_from = vars(KM)["TEMPERATURE"]
        model_from = cs_data_structure.Model()
        model_part_from = model_from.CreateModelPart('wall_from')
        model_part_from.AddNodalSolutionStepVariable(var_from)

        x = np.linspace(0, 10, 11)
        y = np.linspace(0, 20, 21)
        z = np.linspace(0, 30, 31)
        X, Y, Z = np.meshgrid(x, y, z)
        X, Y, Z = X.flatten(), Y.flatten(), Z.flatten()

        for i in range(X.size):
            node = model_part_from.CreateNewNode(i, X[i], Y[i], Z[i])
            node.SetSolutionStepValue(var_from, 0, np.random.rand())

        dX = (np.random.rand(X.size) - .5) * .49
        dY = (np.random.rand(X.size) - .5) * .49
        dZ = (np.random.rand(X.size) - .5) * .49

        Xb, Yb, Zb = X + dX, Y + dY, Z + dZ

        var_to = vars(KM)["PRESSURE"]
        model_to = cs_data_structure.Model()
        model_part_to = model_to.CreateModelPart('wall_to')
        model_part_to.AddNodalSolutionStepVariable(var_to)

        for i in range(Xb.size):
            model_part_to.CreateNewNode(i, Xb[i], Yb[i], Zb[i])

        mapper = cs_tools.CreateInstance(parameters['mapper'])
        mapper.Initialize(model_part_from, model_part_to)
        mapper((model_part_from, var_from), (model_part_to, var_to))

        for node_from, node_to in zip(model_part_from.Nodes, model_part_to.Nodes):
            self.assertEqual(node_from.GetSolutionStepValue(var_from),
                             node_to.GetSolutionStepValue(var_to))


if __name__ == '__main__':
    KratosUnittest.main()
