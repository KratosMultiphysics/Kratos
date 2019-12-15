import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

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


class TestMapperRadialBasis(KratosUnittest.TestCase):

    def test_mapper_radial_basis(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_radial_basis.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())

        # 1D: values on straight line, irregular grid spacing
        var_from = vars(KM)["TEMPERATURE"]
        model_from = cs_data_structure.Model()
        model_part_from = model_from.CreateModelPart('wall_from')
        model_part_from.AddNodalSolutionStepVariable(var_from)

        n_from = 66
        z_from = np.sqrt(np.linspace(.1, .9, n_from))
        v_from = 1.2 - 1.7 * z_from
        for i in range(n_from):
            node = model_part_from.CreateNewNode(i, 0., 0., z_from[i])
            node.SetSolutionStepValue(var_from, 0, v_from[i])

        var_to = vars(KM)["PRESSURE"]
        model_to = cs_data_structure.Model()
        model_part_to = model_to.CreateModelPart('wall_to')
        model_part_to.AddNodalSolutionStepVariable(var_to)

        n_to = 99
        z_to = np.sqrt(np.linspace(0, 1, n_to))
        for i in range(n_to):
            model_part_to.CreateNewNode(i, 0., 0., z_to[i])

        mapper = cs_tools.CreateInstance(parameters['mapper_1d'])
        mapper.Initialize(model_part_from, model_part_to)
        mapper((model_part_from, var_from), (model_part_to, var_to))

        v_to = np.zeros(n_to)
        for i, node in enumerate(model_part_to.Nodes):
            v_to[i] = node.GetSolutionStepValue(var_to)

        # for i in range(n_to):
        #     self.assertAlmostEqual(v_to[i] / (1.2 - 1.7 * z_to[i]), 1., delta=1e-8)

        if True:
            # visualization of results
            plt.plot(z_from, v_from, 'b', label='from', linewidth=2.5)
            plt.scatter(z_from, v_from, s=20, color='b')
            plt.plot(z_to, v_to, 'r', label='to', linewidth=1.5)
            plt.scatter(z_to, v_to, s=15, color='r')
            plt.legend()
            plt.show()
            plt.close()


        # 2D: values on straight line, irregular grid spacing
        var_from = vars(KM)["TEMPERATURE"]
        model_from = cs_data_structure.Model()
        model_part_from = model_from.CreateModelPart('wall_from')
        model_part_from.AddNodalSolutionStepVariable(var_from)

        n_from = 15
        x_from = np.sqrt(np.linspace(.1, .9, n_from))
        y_from = 7.2 + 3.3 * x_from
        v_from = 1.2 - 1.7 * x_from + 2.5 * y_from
        for i in range(n_from):
            node = model_part_from.CreateNewNode(i, x_from[i], y_from[i], 0.)
            node.SetSolutionStepValue(var_from, 0, v_from[i])

        var_to = vars(KM)["PRESSURE"]
        model_to = cs_data_structure.Model()
        model_part_to = model_to.CreateModelPart('wall_to')
        model_part_to.AddNodalSolutionStepVariable(var_to)

        n_to = 30
        x_to = np.sqrt(np.linspace(.2, .8, n_to))
        y_to = 7.2 + 3.3 * x_to
        for i in range(n_to):
            model_part_to.CreateNewNode(i, x_to[i], y_to[i], 0.)

        mapper = cs_tools.CreateInstance(parameters['mapper_2d'])
        mapper.Initialize(model_part_from, model_part_to)
        mapper((model_part_from, var_from), (model_part_to, var_to))

        v_to = np.zeros(n_to)
        for i, node in enumerate(model_part_to.Nodes):
            v_to[i] = node.GetSolutionStepValue(var_to)

        # for i in range(n_to):
        #     self.assertAlmostEqual(v_to[i] / (1.2 - 1.7 * x_to[i] + 2.5 * y_to[i]), 1., delta=1e-8)

        if True:
            # visualization of results
            _, ax = plt.subplots(ncols=2)

            ax[0].plot(x_from, v_from, 'b', label='from', linewidth=2.5)
            ax[0].scatter(x_from, v_from, s=20, color='b')
            ax[0].plot(x_to, v_to, 'r', label='to', linewidth=1.5)
            ax[0].scatter(x_to, v_to, s=15, color='r')

            ax[1].plot(y_from, v_from, 'b', label='from', linewidth=2.5)
            ax[1].scatter(y_from, v_from, s=20, color='b')
            ax[1].plot(y_to, v_to, 'r', label='to', linewidth=1.5)
            ax[1].scatter(y_to, v_to, s=15, color='r')

            for a in ax:
                a.legend()
            plt.show()
            plt.close()


if __name__ == '__main__':
    KratosUnittest.main()
