import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


class TestMapperPermutation(KratosUnittest.TestCase):
    def test_mapper_permutation(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_axisymmetric_2d_to_3d.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())

        # check if method Initialize works
        if True:
            var = vars(KM)["TEMPERATURE"]
            model = cs_data_structure.Model()
            model_part_from = model.CreateModelPart('wall_from')
            model_part_from.AddNodalSolutionStepVariable(var)

            n_from = 10
            x_from = np.linspace(0, 2 * np.pi, n_from)
            y_from = 1. + 0.2 * np.sin(x_from)
            z_from = np.zeros(n_from)
            v_from = x_from * y_from
            for i in range(n_from):
                node = model_part_from.CreateNewNode(i, x_from[i], y_from[i], z_from[i])
                node.SetSolutionStepValue(var, 0, v_from[i])

            # model_part_from given (forward initialization)
            mapper = cs_tools.CreateInstance(parameters['mapper'])
            model_part_to = mapper.Initialize(model_part_from, forward=True)
            mapper((model_part_from, var), (model_part_to, var))

            n_to = model_part_to.NumberOfNodes()
            x_to = np.zeros(n_to)
            y_to = np.zeros(n_to)
            z_to = np.zeros(n_to)
            v_to = np.zeros(n_to)

            for i, node in enumerate(model_part_to.Nodes):
                x_to[i], y_to[i], z_to[i] = node.X0, node.Y0, node.Z0
                v_to[i] = node.GetSolutionStepValue(var)

            # plot results
            c_from = cm.jet((v_from - v_from.min()) / (v_from.max() - v_from.min()))
            c_to = cm.jet((v_to - v_from.min()) / (v_from.max() - v_from.min()))

            fig = plt.figure()

            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(x_from, y_from, z_from, s=50, c=c_from, depthshade=True, marker='s')
            ax.scatter(x_to, y_to, z_to, s=20, c=c_to, depthshade=True)

            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            plt.show()
            plt.close()

            # node = model_part_to.Nodes[0]
            # self.assertListEqual([node.X, node.Y, node.Z], [2., 0., 1.])

            # model_part_to given
            if False:
                mapper = cs_tools.CreateInstance(parameters['mapper'])
                model_part_from = mapper.Initialize(model_part_to, forward=False)
                node = model_part_from.Nodes[0]
                self.assertListEqual([node.X, node.Y, node.Z], [0., 1., 2.])

        # check if method __call__ works for Double Variable
        if False:
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

            # for node_from, node_to in zip(model_part_from.Nodes, model_part_to.Nodes):
            #     val_from = node_from.GetSolutionStepValue(var)
            #     val_to = node_to.GetSolutionStepValue(var)
            #     self.assertEqual(val_from, val_to)

        # check if method __call__ works for Array Variable
        if False:
            var = vars(KM)["DISPLACEMENT"]
            model = cs_data_structure.Model()
            model_part_from = model.CreateModelPart('wall_from')
            model_part_from.AddNodalSolutionStepVariable(var)

            for i in range(10):
                node = model_part_from.CreateNewNode(i, i * 1., i * 2., i * 3.)
                node.SetSolutionStepValue(var, 0, list(np.random.rand(3)))

            mapper = cs_tools.CreateInstance(parameters['mapper'])
            model_part_to = mapper.Initialize(model_part_from, forward=True)
            mapper((model_part_from, var), (model_part_to, var))

            for node_from, node_to in zip(model_part_from.Nodes, model_part_to.Nodes):
                val_from = node_from.GetSolutionStepValue(var)
                val_from = list(np.array(val_from)[mapper.permutation])
                val_to = node_to.GetSolutionStepValue(var)
                self.assertListEqual(val_from, val_to)


if __name__ == '__main__':
    KratosUnittest.main()
