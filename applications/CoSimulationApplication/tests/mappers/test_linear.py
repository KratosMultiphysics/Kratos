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


class TestMapperLinear(KratosUnittest.TestCase):
    def test_mapper_linear(self):
        self.test_mapper_linear_1d()
        self.test_mapper_linear_2d()
        self.test_mapper_linear_3d()

    def test_mapper_linear_1d(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_linear.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())

        # values on straight line, irregular grid spacing
        var_from = vars(KM)["TEMPERATURE"]
        model_from = cs_data_structure.Model()
        model_part_from = model_from.CreateModelPart('wall_from')
        model_part_from.AddNodalSolutionStepVariable(var_from)

        n_from = 6
        z_from = np.sqrt(np.linspace(.1, .9, n_from))
        v_from = 1.2 - 1.7 * z_from
        for i in range(n_from):
            node = model_part_from.CreateNewNode(i, 0., 0., z_from[i])
            node.SetSolutionStepValue(var_from, 0, v_from[i])

        var_to = vars(KM)["PRESSURE"]
        model_to = cs_data_structure.Model()
        model_part_to = model_to.CreateModelPart('wall_to')
        model_part_to.AddNodalSolutionStepVariable(var_to)

        n_to = 9
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

        if 0:
            # visualization of results
            plt.plot(z_from, v_from, 'b', label='from', linewidth=2.5)
            plt.scatter(z_from, v_from, s=20, color='b')
            plt.plot(z_to, v_to, 'r', label='to', linewidth=1.5)
            plt.scatter(z_to, v_to, s=15, color='r')
            plt.legend()
            plt.show()
            plt.close()

    def test_mapper_linear_2d(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_linear.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())

        # values on straight line, irregular grid spacing
        var_from = vars(KM)["TEMPERATURE"]
        model_from = cs_data_structure.Model()
        model_part_from = model_from.CreateModelPart('wall_from')
        model_part_from.AddNodalSolutionStepVariable(var_from)

        n_from = 6
        x_from = np.sqrt(np.linspace(.1, .9, n_from))
        # x_from[4:6] = x_from[1:3]
        y_from = 7.2 + 3.3 * x_from
        v_from = 1.2 - 1.7 * x_from + 2.5 * y_from
        for i in range(n_from):
            node = model_part_from.CreateNewNode(i, x_from[i], y_from[i], 0.)
            node.SetSolutionStepValue(var_from, 0, v_from[i])

        var_to = vars(KM)["PRESSURE"]
        model_to = cs_data_structure.Model()
        model_part_to = model_to.CreateModelPart('wall_to')
        model_part_to.AddNodalSolutionStepVariable(var_to)

        n_to = 9
        x_to = np.sqrt(np.linspace(.11, .95, n_to))
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

        self.assertRaises()

        if 0:
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

    # def test_mapper_linear_3d(self):
    #     parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_linear.json')
    #     cs_data_structure = ImportDataStructure(parameter_file_name)
    #     with open(parameter_file_name, 'r') as parameter_file:
    #         parameters = cs_data_structure.Parameters(parameter_file.read())
    #
    #     # values on flat plane, irregular grid spacing
    #     var_from = vars(KM)["TEMPERATURE"]
    #     model_from = cs_data_structure.Model()
    #     model_part_from = model_from.CreateModelPart('wall_from')
    #     model_part_from.AddNodalSolutionStepVariable(var_from)
    #
    #     nx_from = 105
    #     ny_from = 95
    #     n_from = nx_from * ny_from
    #     x_from = (np.ones((nx_from, ny_from)) * np.sqrt(np.linspace(.5, 1, nx_from)).reshape(-1, 1)).flatten()
    #     y_from = (np.ones((nx_from, ny_from)) * np.sqrt(np.linspace(.5, 1, ny_from)).reshape(1, -1)).flatten()
    #     z_from = 1.1 * x_from - 0.9 * y_from
    #     v_from = 1.2 - 1.7 * x_from + 2.5 * y_from + 0.5 * z_from
    #     for i in range(n_from):
    #         node = model_part_from.CreateNewNode(i, x_from[i], y_from[i], z_from[i])
    #         node.SetSolutionStepValue(var_from, 0, v_from[i])
    #
    #     var_to = vars(KM)["PRESSURE"]
    #     model_to = cs_data_structure.Model()
    #     model_part_to = model_to.CreateModelPart('wall_to')
    #     model_part_to.AddNodalSolutionStepVariable(var_to)
    #
    #     nx_to = 95
    #     ny_to = 105
    #     n_to = nx_to * ny_to
    #     x_to = (np.ones((nx_to, ny_to)) * np.sqrt(np.linspace(.5, 1, nx_to)).reshape(-1, 1)).flatten()
    #     y_to = (np.ones((nx_to, ny_to)) * np.sqrt(np.linspace(.5, 1, ny_to)).reshape(1, -1)).flatten()
    #     z_to = 1.1 * x_to - 0.9 * y_to
    #     for i in range(n_to):
    #         model_part_to.CreateNewNode(i, x_to[i], y_to[i], z_to[i])
    #
    #     mapper = cs_tools.CreateInstance(parameters['mapper_3d'])
    #     with timer('init'):
    #         mapper.Initialize(model_part_from, model_part_to)
    #     with timer('map'):
    #         mapper((model_part_from, var_from), (model_part_to, var_to))
    #
    #     v_to = np.zeros(n_to)
    #     for i, node in enumerate(model_part_to.Nodes):
    #         v_to[i] = node.GetSolutionStepValue(var_to)
    #
    #     if 0:
    #         vmin, vmax = v_to.min(), v_to.max()
    #         fig = plt.figure()
    #         ax = fig.add_subplot(111, projection='3d')
    #         ax.scatter(x_to, y_to, z_to, c=v_to, vmin=vmin, vmax=vmax, cmap=cm.coolwarm)
    #         _ = ax.scatter(x_from, y_from, z_from, s=50, c=v_from, vmin=vmin, vmax=vmax, cmap=cm.coolwarm)
    #         fig.colorbar(_, ax=ax)
    #         ax.set_xlabel('x')
    #         ax.set_ylabel('y')
    #         ax.set_zlabel('z')
    #         plt.show()
    #         plt.close('all')
    #
    #     # for i in range(n_to):
    #     #     self.assertAlmostEqual(v_to[i] / (1.2 - 1.7 * x_to[i] + 2.5 * y_to[i] + 0.5 * z_to[i]), 1., delta=1e-8)
    #
    #
    #     # interpolation on 3D sinc function
    #     var_from = vars(KM)["TEMPERATURE"]
    #     model_from = cs_data_structure.Model()
    #     model_part_from = model_from.CreateModelPart('wall_from')
    #     model_part_from.AddNodalSolutionStepVariable(var_from)
    #
    #     nx_from = 95
    #     ny_from = 105
    #     n_from = nx_from * ny_from
    #     x_from = (np.ones((nx_from, ny_from)) * np.linspace(-10, 10, nx_from).reshape(-1, 1)).flatten()
    #     y_from = (np.ones((nx_from, ny_from)) * np.linspace(-10, 10, ny_from).reshape(1, -1)).flatten()
    #     z_from = np.sinc(np.sqrt(x_from ** 2 + y_from ** 2) / np.pi)
    #     v_from = z_from + x_from + y_from
    #     for i in range(n_from):
    #         node = model_part_from.CreateNewNode(i, x_from[i], y_from[i], z_from[i])
    #         node.SetSolutionStepValue(var_from, 0, v_from[i])
    #
    #     var_to = vars(KM)["PRESSURE"]
    #     model_to = cs_data_structure.Model()
    #     model_part_to = model_to.CreateModelPart('wall_to')
    #     model_part_to.AddNodalSolutionStepVariable(var_to)
    #
    #     nx_to = 105
    #     ny_to = 95
    #     n_to = nx_to * ny_to
    #     x_to = (np.ones((nx_to, ny_to)) * np.linspace(-10, 10, nx_to).reshape(-1, 1)).flatten()
    #     y_to = (np.ones((nx_to, ny_to)) * np.linspace(-10, 10, ny_to).reshape(1, -1)).flatten()
    #     z_to = np.sinc(np.sqrt(x_to ** 2 + y_to ** 2) / np.pi)
    #     for i in range(n_to):
    #         model_part_to.CreateNewNode(i, x_to[i], y_to[i], z_to[i])
    #
    #     mapper = cs_tools.CreateInstance(parameters['mapper_3d'])
    #     with timer('init'):
    #         mapper.Initialize(model_part_from, model_part_to)
    #     with timer('map'):
    #         mapper((model_part_from, var_from), (model_part_to, var_to))
    #
    #     v_to = np.zeros(n_to)
    #     for i, node in enumerate(model_part_to.Nodes):
    #         v_to[i] = node.GetSolutionStepValue(var_to)
    #
    #     if 0:
    #         vmin, vmax = v_from.min(), v_from.max()
    #
    #         fig = plt.figure()
    #         ax = fig.add_subplot(111, projection='3d')
    #         _ = ax.scatter(x_from, y_from, z_from, c=v_from, vmin=vmin, vmax=vmax, cmap=cm.coolwarm)
    #         # _ = ax.plot_trisurf(x_from, y_from, z_from, color=v_from, vmin=vmin, vmax=vmax,
    #         #                     cmap=cm.coolwarm)
    #         fig.colorbar(_, ax=ax)
    #
    #         fig = plt.figure()
    #         ax = fig.add_subplot(111, projection='3d')
    #         _ = ax.scatter(x_to, y_to, z_to, c=v_to, vmin=vmin, vmax=vmax, cmap=cm.coolwarm)
    #         # _ = ax.plot_trisurf(x_to, y_to, z_to, color=v_to, vmin=vmin, vmax=vmax,
    #         #                     cmap=cm.coolwarm)
    #         fig.colorbar(_, ax=ax)
    #         # ax.set_xlabel('x')
    #         # ax.set_ylabel('y')
    #         # ax.set_zlabel('z')
    #         plt.show()
    #         plt.close('all')
    #
    #
    #     # interpolation of 3D sine function on sphere
    #     def f(x, y, z):
    #         # out = 1.5 * x + 2 * y + 2.5 * z
    #         out = np.sin(x) * np.sin(y) * np.sin(z)
    #         return out
    #
    #     r = np.pi
    #
    #     # ModelPart to
    #     var_from = vars(KM)["TEMPERATURE"]
    #     model_from = cs_data_structure.Model()
    #     model_part_from = model_from.CreateModelPart('wall_from')
    #     model_part_from.AddNodalSolutionStepVariable(var_from)
    #
    #     n_theta = 95
    #     n_phi = 105
    #     dtheta = np.pi / n_theta
    #     dphi = np.pi / (n_phi - 1)
    #     theta = np.ones((n_theta, n_phi)) * np.linspace(0, 2 * np.pi - dtheta, n_theta).reshape(-1, 1)
    #     phi = np.ones((n_theta, n_phi)) * np.linspace(dphi, np.pi - dphi, n_phi).reshape(1, -1)
    #
    #     n_from = n_theta * n_phi
    #     x_from = r * np.cos(theta) * np.sin(phi)
    #     y_from = r * np.sin(theta) * np.sin(phi)
    #     z_from = r * np.cos(phi)
    #     v_from = f(x_from, y_from, z_from)
    #     for i in range(n_from):
    #         node = model_part_from.CreateNewNode(i, x_from.flatten()[i], y_from.flatten()[i], z_from.flatten()[i])
    #         node.SetSolutionStepValue(var_from, 0, v_from.flatten()[i])
    #
    #     # ModelPart from
    #     var_to = vars(KM)["PRESSURE"]
    #     model_to = cs_data_structure.Model()
    #     model_part_to = model_to.CreateModelPart('wall_to')
    #     model_part_to.AddNodalSolutionStepVariable(var_to)
    #
    #     n_theta = 105
    #     n_phi = 95
    #     dtheta = np.pi / n_theta
    #     dphi = np.pi / (n_phi - 1)
    #     theta = np.ones((n_theta, n_phi)) * np.linspace(0, 2 * np.pi - dtheta, n_theta).reshape(-1, 1)
    #     phi = np.ones((n_theta, n_phi)) * np.linspace(dphi, np.pi - dphi, n_phi).reshape(1, -1)
    #
    #     n_to = n_theta * n_phi
    #     x_to = r * np.cos(theta) * np.sin(phi)
    #     y_to = r * np.sin(theta) * np.sin(phi)
    #     z_to = r * np.cos(phi)
    #     v_to_ref = f(x_to, y_to, z_to)
    #     for i in range(n_to):
    #         model_part_to.CreateNewNode(i, x_to.flatten()[i], y_to.flatten()[i], z_to.flatten()[i])
    #
    #     # map values
    #     mapper = cs_tools.CreateInstance(parameters['mapper_3d'])
    #     with timer('init', ms=True):
    #         mapper.Initialize(model_part_from, model_part_to)
    #     with timer('map', ms=True):
    #         mapper((model_part_from, var_from), (model_part_to, var_to))
    #
    #     v_to = np.zeros(n_to)
    #     for i, node in enumerate(model_part_to.Nodes):
    #         v_to[i] = node.GetSolutionStepValue(var_to)
    #     v_to = v_to.reshape(n_theta, n_phi)
    #
    #     # plot results
    #     if 0:
    #         c_from = cm.jet((v_from - v_from.min()) / (v_from.max() - v_from.min()))
    #         c_to = cm.jet((v_to - v_from.min()) / (v_from.max() - v_from.min()))
    #         error = np.abs(v_to - v_to_ref)
    #         c_error = cm.jet(error / error.max())
    #
    #         fig = plt.figure(figsize=(18, 6))
    #         plt.suptitle(f'max error = {error.max():.2g}     ({v_from.min():.1f} < v_from < {v_from.max():.1f})')
    #
    #         ax_from = fig.add_subplot(131, projection='3d')
    #         ax_from.set_title('from')
    #         ax_from.plot_surface(x_from, y_from, z_from, facecolors=c_from,
    #                              rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)
    #
    #         ax_to = fig.add_subplot(132, projection='3d')
    #         ax_to.set_title('to')
    #         ax_to.plot_surface(x_to, y_to, z_to, facecolors=c_to,
    #                            rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)
    #
    #         ax_error = fig.add_subplot(133, projection='3d')
    #         ax_error.set_title('to (error)')
    #         ax_error.plot_surface(x_to, y_to, z_to, facecolors=c_error,
    #                               rstride=1, cstride=1, antialiased=False, shade=False)
    #
    #         for ax in [ax_from, ax_to, ax_error]:
    #             ax.set_xlabel('x')
    #             ax.set_ylabel('y')
    #             ax.set_zlabel('z')
    #
    #         plt.tight_layout()
    #         plt.show()
    #         plt.close('all')


if __name__ == '__main__':
    KratosUnittest.main()
