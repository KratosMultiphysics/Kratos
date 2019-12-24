import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


class TestMapperNearest(KratosUnittest.TestCase):
    def test_mapper_nearest(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_nearest.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())
        par_mapper = parameters['mapper']

        gui = 0  # *** gui gives problems when running all tests?

        # 1D case: square-root grid + linear function
        """
        n_from = 14, n_to = 5 
            => max error = 0.0085
        """
        n_from, n_to = 14, 5
        par_mapper['settings'].SetArray('directions', ['Z'])

        case = Case1D(cs_data_structure, n_from, n_to)
        case.map(cs_tools, par_mapper)
        self.assertTrue(case.check(tolerance=0.01))
        if gui:
            case.plot()

        # 2D case: circle + linear function
        """
        n_from = 33, n_to = 22 
            => max error = 0.55
        """
        n_from, n_to = 33, 22
        par_mapper['settings'].SetArray('directions', ['X', 'Y'])

        case = Case2D(cs_data_structure, n_from, n_to)
        case.map(cs_tools, par_mapper)
        self.assertTrue(case.check(tolerance=1.))
        if gui:
            case.plot()

        # 3D case : sphere + sine function
        """
        n_theta_from, n_phi_from = 40, 20
        n_theta_to, n_phi_to = 50, 21
            => max error = 0.16

        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
            => max error = 0.13
        """
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])

        case = Case3DSphere(cs_data_structure, n_theta_from, n_phi_from, n_theta_to, n_phi_to)
        case.map(cs_tools, par_mapper)
        self.assertTrue(case.check(tolerance=0.2))
        if gui:
            case.plot()

        # 3D case: sinc + linear vector function
        """
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
            => max error = 0.13
        """
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])

        case = Case3DSinc(cs_data_structure, n_x_from, n_y_from, n_x_to, n_y_to)
        case.map(cs_tools, par_mapper)
        for tmp in case.check(tolerance=0.2):
            self.assertTrue(tmp)
        if gui:
            case.plot()


class Case1D:
    # 1D case: square-root grid + linear function
    def __init__(self, cs_data_structure, n_from, n_to):
        self.n_from = n_from
        self.n_to = n_to

        model = cs_data_structure.Model()

        # ModelPart from
        self.var_from = vars(KM)["TEMPERATURE"]
        self.model_part_from = model.CreateModelPart('wall_from')
        self.model_part_from.AddNodalSolutionStepVariable(self.var_from)
        self.z_from = np.linspace(0, 10, self.n_from) ** .5
        self.v_from = self.fun(self.z_from)
        for i in range(self.n_from):
            node = self.model_part_from.CreateNewNode(i, 0., 0., self.z_from[i])
            node.SetSolutionStepValue(self.var_from, 0, self.v_from[i])

        # ModelPart to
        self.var_to = vars(KM)["PRESSURE"]
        self.model_part_to = model.CreateModelPart('wall_to')
        self.model_part_to.AddNodalSolutionStepVariable(self.var_to)
        self.z_to = np.linspace(0, 10, self.n_to) ** .5
        for i in range(self.n_to):
            self.model_part_to.CreateNewNode(i, 0., 0., self.z_to[i])

    def map(self, cs_tools, parameters):
        mapper = cs_tools.CreateInstance(parameters)
        mapper.Initialize(self.model_part_from, self.model_part_to)
        mapper((self.model_part_from, self.var_from),
               (self.model_part_to, self.var_to))

        self.v_to_fun = self.fun(self.z_to)
        self.v_to = np.zeros(self.n_to)
        for i, node in enumerate(self.model_part_to.Nodes):
            self.v_to[i] = node.GetSolutionStepValue(self.var_to)
        self.v_error = np.abs(self.v_to - self.v_to_fun)

    def check(self, tolerance):
        criterion = (self.v_error < tolerance)
        return criterion.all()

    def plot(self):
        _, ax = plt.subplots(ncols=2, sharex=True, figsize=(15, 6))
        plt.title(f'max error = {self.v_error.max():.2g}')

        ax[0].plot(self.z_from, self.v_from, label='from', marker='o')
        ax[0].plot(self.z_to, self.v_to, label='to', marker='o')
        ax[1].plot(self.z_to, self.v_error, label='error', marker='o')
        for a in ax:
            a.legend()
            a.set_xlabel('z')
            a.set_ylabel('f(z)')

        plt.tight_layout()
        plt.show()
        plt.close()

    def fun(self, z):
        return z / 10


class Case2D:
    # 2D case: circle + linear function
    def __init__(self, cs_data_structure, n_from, n_to):
        self.n_from = n_from
        self.n_to = n_to

        model = cs_data_structure.Model()

        # ModelPart from
        self.var_from = vars(KM)["TEMPERATURE"]
        self.model_part_from = model.CreateModelPart('wall_from')
        self.model_part_from.AddNodalSolutionStepVariable(self.var_from)

        dtheta = 2 * np.pi / self.n_from
        self.theta_from = np.linspace(0, 2 * np.pi - dtheta, self.n_from)

        self.x_from, self.y_from = self.get_cartesian(self.theta_from)
        self.v_from = self.fun(self.x_from, self.y_from)
        for i in range(self.n_from):
            node = self.model_part_from.CreateNewNode(i, self.x_from[i], self.y_from[i], 0.)
            node.SetSolutionStepValue(self.var_from, 0, self.v_from[i])

        # ModelPart to
        self.var_to = vars(KM)["PRESSURE"]
        self.model_part_to = model.CreateModelPart('wall_to')
        self.model_part_to.AddNodalSolutionStepVariable(self.var_to)

        dtheta = 2 * np.pi / self.n_to
        self.theta_to = np.linspace(0, 2 * np.pi - dtheta, self.n_to)

        self.x_to, self.y_to = self.get_cartesian(self.theta_to)
        for i in range(self.n_to):
            self.model_part_to.CreateNewNode(i, self.x_to[i], self.y_to[i], 0.)

    def map(self, cs_tools, parameters):
        mapper = cs_tools.CreateInstance(parameters)
        mapper.Initialize(self.model_part_from, self.model_part_to)
        mapper((self.model_part_from, self.var_from),
               (self.model_part_to, self.var_to))

        self.v_to_fun = self.fun(self.x_to, self.y_to)
        self.v_to = np.zeros(self.n_to)
        for i, node in enumerate(self.model_part_to.Nodes):
            self.v_to[i] = node.GetSolutionStepValue(self.var_to)
        self.v_error = np.abs(self.v_to - self.v_to_fun)

    def check(self, tolerance):
        criterion = (self.v_error < tolerance)
        return criterion.all()

    def plot(self):
        _, ax = plt.subplots(ncols=2, sharex=True, figsize=(15, 6))

        ax[0].plot(self.theta_from * 180 / np.pi, self.v_from, label='from', marker='o')
        ax[0].plot(self.theta_to * 180 / np.pi, self.v_to, label='to', marker='o')

        ax[1].plot(self.theta_to * 180 / np.pi, self.v_error, label='error', marker='o')

        for a in ax:
            a.legend()
            a.set_xlabel(r'$\theta$')
            a.set_ylabel(r'f($\theta$)')

        plt.tight_layout()
        plt.show()
        plt.close()

    def fun(self, x, y):
        return 2 * x + 3 * y

    def get_cartesian(self, theta):
        r = 2.
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        return x, y


class Case3DSphere:
    # 3D case: sphere + sine function
    def __init__(self, cs_data_structure, n_theta_from, n_phi_from, n_theta_to, n_phi_to):
        self.n_theta_from = n_theta_from
        self.n_phi_from = n_phi_from
        self.n_from = n_theta_from * n_phi_from
        self.n_theta_to = n_theta_to
        self.n_phi_to = n_phi_to  # for bounding box: not too far from n_phi_from!
        self.n_to = n_theta_to * n_phi_to

        model = cs_data_structure.Model()

        # ModelPart from
        self.var_from = vars(KM)["TEMPERATURE"]
        self.model_part_from = model.CreateModelPart('wall_from')
        self.model_part_from.AddNodalSolutionStepVariable(self.var_from)

        shape = (self.n_theta_from, self.n_phi_from)
        dtheta = 2 * np.pi / self.n_theta_from
        dphi = np.pi / (self.n_phi_from - 1)
        theta = np.ones(shape) * np.linspace(0, 2 * np.pi - dtheta, self.n_theta_from).reshape(-1, 1)
        phi = np.ones(shape) * np.linspace(dphi, np.pi - dphi, self.n_phi_from).reshape(1, -1)

        self.x_from, self.y_from, self.z_from = self.get_cartesian(theta, phi)
        self.v_from = self.fun(self.x_from, self.y_from, self.z_from)
        for i in range(self.n_from):
            node = self.model_part_from.CreateNewNode(i, self.x_from.flatten()[i],
                                self.y_from.flatten()[i], self.z_from.flatten()[i])
            node.SetSolutionStepValue(self.var_from, 0, self.v_from.flatten()[i])

        # ModelPart to
        self.var_to = vars(KM)["PRESSURE"]
        self.model_part_to = model.CreateModelPart('wall_to')
        self.model_part_to.AddNodalSolutionStepVariable(self.var_to)

        shape = (self.n_theta_to, self.n_phi_to)
        dtheta = 2 * np.pi / self.n_theta_to
        dphi = np.pi / (self.n_phi_to - 1)
        theta = np.ones(shape) * np.linspace(0, 2 * np.pi - dtheta, self.n_theta_to).reshape(-1, 1)
        phi = np.ones(shape) * np.linspace(dphi, np.pi - dphi, self.n_phi_to).reshape(1, -1)

        self.x_to, self.y_to, self.z_to = self.get_cartesian(theta, phi)
        for i in range(self.n_to):
            self.model_part_to.CreateNewNode(i, self.x_to.flatten()[i],
                            self.y_to.flatten()[i], self.z_to.flatten()[i])

    def map(self, cs_tools, parameters):
        mapper = cs_tools.CreateInstance(parameters)
        mapper.Initialize(self.model_part_from, self.model_part_to)
        mapper((self.model_part_from, self.var_from),
               (self.model_part_to, self.var_to))

        self.v_to_fun = self.fun(self.x_to, self.y_to, self.z_to)
        self.v_to = np.zeros(self.n_to)
        for i, node in enumerate(self.model_part_to.Nodes):
            self.v_to[i] = node.GetSolutionStepValue(self.var_to)
        self.v_to = self.v_to.reshape(self.x_to.shape)
        self.v_error = np.abs(self.v_to - self.v_to_fun)

    def check(self, tolerance):
        criterion = (self.v_error < tolerance)
        return criterion.all()

    def plot(self):
        v_min = min(self.v_from.min(), self.v_to.min())
        v_max = max(self.v_from.max(), self.v_to.max())
        c_from = cm.jet((self.v_from - v_min) / (v_max - v_min))
        c_to = cm.jet((self.v_to - v_min) / (v_max - v_min))
        c_error = cm.jet(self.v_error / self.v_error.max())

        fig = plt.figure(figsize=(18, 6))
        plt.suptitle(f'max error = {self.v_error.max():.2g}     ({v_min:.1f} < v < {v_max:.1g})')

        ax_from = fig.add_subplot(131, projection='3d')
        ax_from.set_title('from')
        ax_from.plot_surface(self.x_from, self.y_from, self.z_from, facecolors=c_from,
                             rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

        ax_to = fig.add_subplot(132, projection='3d')
        ax_to.set_title('to')
        ax_to.plot_surface(self.x_to, self.y_to, self.z_to, facecolors=c_to,
                           rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

        ax_error = fig.add_subplot(133, projection='3d')
        ax_error.set_title('to (error)')
        ax_error.plot_surface(self.x_to, self.y_to, self.z_to, facecolors=c_error,
                              rstride=1, cstride=1, antialiased=False, shade=False)

        for ax in [ax_from, ax_to, ax_error]:
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')

        plt.tight_layout()
        plt.show()
        plt.close()

    def fun(self, x, y, z):
        return np.sin(x) * np.sin(y) * np.sin(z)

    def get_cartesian(self, theta, phi):
        r = np.pi
        x = r * np.cos(theta) * np.sin(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(phi)
        return x, y, z


class Case3DSinc:
    # 3D case: sinc + linear vector function
    def __init__(self, cs_data_structure, n_x_from, n_y_from, n_x_to, n_y_to):
        self.n_x_from = n_x_from
        self.n_y_from = n_y_from
        self.n_from = n_x_from * n_y_from
        self.n_x_to = n_x_to
        self.n_y_to = n_y_to
        self.n_to = n_x_to * n_y_to

        model = cs_data_structure.Model()

        # ModelPart from
        self.var_from = vars(KM)["VELOCITY"]
        self.model_part_from = model.CreateModelPart('wall_from')
        self.model_part_from.AddNodalSolutionStepVariable(self.var_from)

        shape = (self.n_x_from, self.n_y_from)
        self.x_from = np.ones(shape) * np.linspace(-1, 1, self.n_x_from).reshape(-1, 1)
        self.y_from = np.ones(shape) * np.linspace(-1, 1, self.n_y_from).reshape(1, -1)
        self.z_from = np.sinc(np.sqrt((10 * self.x_from) ** 2 + (self.y_from * 10) ** 2) / np.pi)
        self.v_from = self.fun(self.x_from, self.y_from, self.z_from)
        for i in range(self.n_from):
            node = self.model_part_from.CreateNewNode(i, self.x_from.flatten()[i],
                                    self.y_from.flatten()[i], self.z_from.flatten()[i])
            hist = []
            for j in range(3):
                hist.append(self.v_from[j].flatten()[i])
            node.SetSolutionStepValue(self.var_from, 0, hist)

        # ModelPart to
        self.var_to = vars(KM)["FORCE"]
        self.model_part_to = model.CreateModelPart('wall_to')
        self.model_part_to.AddNodalSolutionStepVariable(self.var_to)

        shape = (self.n_x_to, self.n_y_to)
        self.x_to = np.ones(shape) * np.linspace(-1, 1, self.n_x_to).reshape(-1, 1)
        self.y_to = np.ones(shape) * np.linspace(-1, 1, self.n_y_to).reshape(1, -1)
        self.z_to = np.sinc(np.sqrt((10 * self.x_to) ** 2 + (self.y_to * 10) ** 2) / np.pi)
        self.v_to = self.fun(self.x_to, self.y_to, self.z_to)
        for i in range(self.n_to):
            self.model_part_to.CreateNewNode(i, self.x_to.flatten()[i],
                                 self.y_to.flatten()[i], self.z_to.flatten()[i])

    def map(self, cs_tools, parameters):
        mapper = cs_tools.CreateInstance(parameters)
        mapper.Initialize(self.model_part_from, self.model_part_to)
        mapper((self.model_part_from, self.var_from),
               (self.model_part_to, self.var_to))

        self.v_to_fun = self.fun(self.x_to, self.y_to, self.z_to)
        self.v_to = [np.zeros(self.n_to), np.zeros(self.n_to), np.zeros(self.n_to)]
        for i, node in enumerate(self.model_part_to.Nodes):
            hist = node.GetSolutionStepValue(self.var_to)
            for j in range(3):
                self.v_to[j][i] = hist[j]

        self.v_error = []
        for j in range(3):
            self.v_to[j] = self.v_to[j].reshape(self.x_to.shape)
            self.v_error.append(np.abs(self.v_to[j] - self.v_to_fun[j]))

    def check(self, tolerance):
        out = []
        for j in range(3):
            criterion = (self.v_error[j] < tolerance)
            out.append(criterion.all())
        return out

    def plot(self):
        fig = plt.figure(figsize=(18, 10))
        for j in range(3):
            v_min = min(self.v_from[j].min(), self.v_to[j].min())
            v_max = max(self.v_from[j].max(), self.v_to[j].max())
            c_from = cm.jet((self.v_from[j] - v_min) / (v_max - v_min))
            c_to = cm.jet((self.v_to[j] - v_min) / (v_max - v_min))
            c_error = cm.jet(self.v_error[j] / self.v_error[j].max())

            ax_from = fig.add_subplot(3, 3, 1 + j * 3, projection='3d')
            ax_from.set_title('from')
            ax_from.plot_surface(self.x_from, self.y_from, self.z_from, facecolors=c_from,
                                 rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

            ax_to = fig.add_subplot(3, 3, 2 + j * 3, projection='3d')
            ax_to.set_title('to')
            ax_to.plot_surface(self.x_to, self.y_to, self.z_to, facecolors=c_to,
                               rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

            ax_error = fig.add_subplot(3, 3, 3 + j * 3, projection='3d')
            ax_error.set_title(f'max error = {self.v_error[j].max():.2g} ({v_min:.1f} < v < {v_max:.1g})')
            ax_error.plot_surface(self.x_to, self.y_to, self.z_to, facecolors=c_error,
                                  rstride=1, cstride=1, antialiased=False, shade=False)

            for ax in [ax_from, ax_to, ax_error]:
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_zlabel('z')

        plt.tight_layout()
        plt.get_current_fig_manager().window.showMaximized()
        plt.show()
        plt.close()

    def fun(self, x, y, z):
        return [y + z, z + x, x + y]


if __name__ == '__main__':
    KratosUnittest.main()
