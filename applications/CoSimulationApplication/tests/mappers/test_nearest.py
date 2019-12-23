import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
import os
from copy import deepcopy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


class TestMapperNearest(KratosUnittest.TestCase):
    def test_mapper_nearest(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_nearest.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())
        par_mapper_0 = parameters['mapper']

        gui = 1

        # 1D: square-root grid + linear function
        """
        n_from = 14, n_to = 5 
            => max error = 0.0085
        """
        n_from, n_to = 14, 5
        case = Case1D(cs_data_structure, n_from, n_to)

        par_mapper = deepcopy(par_mapper_0)
        par_mapper['settings'].SetArray('directions', ['Z'])
        mapper = cs_tools.CreateInstance(par_mapper)
        mapper.Initialize(case.model_part_from, case.model_part_to)
        mapper((case.model_part_from, case.var_from),
               (case.model_part_to, case.var_to))

        self.assertTrue(case.check(tolerance=0.01))

        if gui:
            case.plot()

        # 3D: sphere + sine function
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
        case = Case3DSphere(cs_data_structure, n_theta_from, n_phi_from, n_theta_to, n_phi_to)

        par_mapper = deepcopy(par_mapper_0)
        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])
        mapper = cs_tools.CreateInstance(par_mapper)
        mapper.Initialize(case.model_part_from, case.model_part_to)
        mapper((case.model_part_from, case.var_from),
               (case.model_part_to, case.var_to))

        self.assertTrue(case.check(tolerance=0.2))

        if gui:
            case.plot()



        if False:
            args_from, args_to, fun, shape_to = case_3d_sphere(cs_data_structure)

            par_mapper = deepcopy(par_mapper_0)
            par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])
            mapper = cs_tools.CreateInstance(par_mapper)
            mapper.Initialize(args_from[0], args_to[0])
            mapper(args_from, args_to)

            self.assertTrue(case_3d_sphere_check(args_to, fun))

            if gui:
                case_3d_sphere_plot(args_from, args_to, fun, shape_to)


        # *** do a check with vector variables? perhaps other 3D case?
        # *** e.g. each vector variable is linear in 1 direction


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

    def check(self, tolerance):
        v_to = np.zeros(self.n_to)
        for i, node in enumerate(self.model_part_to.Nodes):
            v_to[i] = node.GetSolutionStepValue(self.var_to)

        v_error = np.abs(v_to - self.fun(self.z_to))
        criterion = (v_error < tolerance)
        return criterion.all()

    def plot(self):
        v_to_fun = self.fun(self.z_to)
        v_to = np.zeros(self.n_to)
        for i, node in enumerate(self.model_part_to.Nodes):
            v_to[i] = node.GetSolutionStepValue(self.var_to)


        _, ax = plt.subplots(ncols=2, sharex=True, figsize=(15, 6))

        ax[0].plot(self.z_from, self.v_from, label='from', marker='o')
        ax[0].plot(self.z_to, v_to, label='to', marker='o')

        ax[1].plot(self.z_to, np.abs(v_to_fun - v_to), label='error', marker='o')

        for a in ax:
            a.legend()
            a.set_xlabel('z')
            a.set_ylabel('f(z)')

        plt.tight_layout()
        plt.show()
        plt.close()

    def fun(self, z):
        return z / 10

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
        r = np.pi

        # ModelPart from
        self.var_from = vars(KM)["TEMPERATURE"]
        self.model_part_from = model.CreateModelPart('wall_from')
        self.model_part_from.AddNodalSolutionStepVariable(self.var_from)

        shape = (self.n_theta_from, self.n_phi_from)
        dtheta = np.pi / self.n_theta_from
        dphi = np.pi / (self.n_phi_from - 1)
        theta = np.ones(shape) * np.linspace(0, 2 * np.pi - dtheta, self.n_theta_from).reshape(-1, 1)
        phi = np.ones(shape) * np.linspace(dphi, np.pi - dphi, self.n_phi_from).reshape(1, -1)

        self.x_from = r * np.cos(theta) * np.sin(phi)
        self.y_from = r * np.sin(theta) * np.sin(phi)
        self.z_from = r * np.cos(phi)
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
        dtheta = np.pi / self.n_theta_to
        dphi = np.pi / (self.n_phi_to - 1)
        theta = np.ones(shape) * np.linspace(0, 2 * np.pi - dtheta, self.n_theta_to).reshape(-1, 1)
        phi = np.ones(shape) * np.linspace(dphi, np.pi - dphi, self.n_phi_to).reshape(1, -1)

        self.x_to = r * np.cos(theta) * np.sin(phi)
        self.y_to = r * np.sin(theta) * np.sin(phi)
        self.z_to = r * np.cos(phi)
        for i in range(self.n_to):
            self.model_part_to.CreateNewNode(i, self.x_to.flatten()[i],
                            self.y_to.flatten()[i], self.z_to.flatten()[i])

    def check(self, tolerance):
        v_to_fun = self.fun(self.x_to, self.y_to, self.z_to)
        v_to = np.zeros(self.n_to)
        for i, node in enumerate(self.model_part_to.Nodes):
            v_to[i] = node.GetSolutionStepValue(self.var_to)
        v_to = v_to.reshape(self.x_to.shape)

        v_error = np.abs(v_to - v_to_fun)
        criterion = (v_error < tolerance)
        return criterion.all()

    def plot(self):
        v_to_fun = self.fun(self.x_to, self.y_to, self.z_to)
        v_to = np.zeros(self.n_to)
        for i, node in enumerate(self.model_part_to.Nodes):
            v_to[i] = node.GetSolutionStepValue(self.var_to)
        v_to = v_to.reshape(self.x_to.shape)

        v_min = min(self.v_from.min(), v_to.min())
        v_max = max(self.v_from.max(), v_to.max())
        c_from = cm.jet((self.v_from - v_min) / (v_max - v_min))
        c_to = cm.jet((v_to - v_min) / (v_max - v_min))
        v_error = np.abs(v_to - v_to_fun)
        c_error = cm.jet(v_error / v_error.max())

        fig = plt.figure(figsize=(18, 6))
        plt.suptitle(f'max error = {v_error.max():.2g}     ({v_min:.1f} < v < {v_max:.1g})')

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




# *** make class object? easier to import! less double work

if __name__ == '__main__':
    KratosUnittest.main()
