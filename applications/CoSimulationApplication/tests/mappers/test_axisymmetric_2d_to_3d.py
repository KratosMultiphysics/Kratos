import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


class TestMapperAxisymmetric2DTo3D(KratosUnittest.TestCase):
    def test_mapper_axisymmetric_2d_to_3d(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_axisymmetric_2d_to_3d.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())

        # check if method Initialize works
        if True:
            n_t = parameters['mapper']['settings']['n_tangential'].GetInt()

            # create 2D model_part_in
            model = cs_data_structure.Model()
            model_part_in = model.CreateModelPart('wall_in')

            n_in = 10
            x_in = np.linspace(0, 2 * np.pi, n_in)
            y_in = 1. + 0.2 * np.sin(x_in)
            for i in range(n_in):
                model_part_in.CreateNewNode(i, x_in[i], y_in[i], 0.)

            # create reference geometry for 3D model_part_out
            n_out_ref = n_in * n_t
            x_out_ref = np.zeros(n_out_ref)
            y_out_ref = np.zeros(n_out_ref)
            z_out_ref = np.zeros(n_out_ref)

            i_to = 0
            for i_from in range(n_in):
                for j in range(n_t):
                    theta = j * 2 * np.pi / n_t
                    x_out_ref[i_to] = x_in[i_from]
                    y_out_ref[i_to] = np.cos(theta) * y_in[i_from]
                    z_out_ref[i_to] = np.sin(theta) * y_in[i_from]
                    i_to += 1
                i_from += 1

            # initialize mapper to get model_part_out
            mapper = cs_tools.CreateInstance(parameters['mapper'])
            model_part_out = mapper.Initialize(model_part_in, forward=True)

            # get mapped geometry from 3D model_part_out
            n_out = model_part_out.NumberOfNodes()
            x_out = np.zeros(n_out)
            y_out = np.zeros(n_out)
            z_out = np.zeros(n_out)
            for i, node in enumerate(model_part_out.Nodes):
                x_out[i], y_out[i], z_out[i] = node.X0, node.Y0, node.Z0

            # compare mapped and reference geometries
            self.assertEqual(n_out, n_out_ref)
            self.assertListEqual(list(x_out), list(x_out_ref))
            self.assertListEqual(list(y_out), list(y_out_ref))
            self.assertListEqual(list(z_out), list(z_out_ref))

        # check if method __call__ works
        if True:
            def fun_s(x):
                return 1. + 2.5 * x

            def fun_v(x, y, z):
                theta = np.arctan2(z, y)
                f_x = 1. + 2.5 * x
                f_y = f_x * 0.5 * np.cos(theta)
                f_z = f_x * 0.5 * np.sin(theta)
                return [f_x, f_y, f_z]

            # create model_part_from (2D)
            var_s = vars(KM)["TEMPERATURE"]
            var_v = vars(KM)["VELOCITY"]
            model = cs_data_structure.Model()
            model_part_from = model.CreateModelPart('wall_from')
            model_part_from.AddNodalSolutionStepVariable(var_s)
            model_part_from.AddNodalSolutionStepVariable(var_v)

            n = 10
            for i in range(n):
                x = .3 * i
                node = model_part_from.CreateNewNode(i, x, 1. + 0.2 * np.sin(x), 0.)
                node.SetSolutionStepValue(var_s, 0, fun_s(node.X0))
                node.SetSolutionStepValue(var_v, 0, fun_v(node.X0, node.Y0, node.Z0))

            # initialize mapper to get model_part_to (3D)
            mapper = cs_tools.CreateInstance(parameters['mapper'])
            model_part_to = mapper.Initialize(model_part_from, forward=True)

            # check mapped values for Double Variable
            mapper((model_part_from, var_s), (model_part_to, var_s))
            for node in model_part_to.Nodes:
                self.assertAlmostEqual(node.GetSolutionStepValue(var_s),
                                       fun_s(node.X0), delta=1e-8)

            # check mapped values for Array Variable
            mapper((model_part_from, var_v), (model_part_to, var_v))
            for node in model_part_to.Nodes:
                for v1, v2 in zip(list(node.GetSolutionStepValue(var_v)),
                                  fun_v(node.X0, node.Y0, node.Z0)):
                    self.assertAlmostEqual(v1, v2, delta=1e-8)

        # extra: visual check of whole method
        if True:
            if os.getcwd() == os.path.dirname(os.path.abspath(__file__)):
                print('\n\nrunning visual check for whole method')

                # create model_part_from (2D)
                var_s = vars(KM)["TEMPERATURE"]
                var_v = vars(KM)["VELOCITY"]
                model = cs_data_structure.Model()
                model_part_from = model.CreateModelPart('wall_from')
                model_part_from.AddNodalSolutionStepVariable(var_s)
                model_part_from.AddNodalSolutionStepVariable(var_v)

                n = 10
                for i in range(n):
                    model_part_from.CreateNewNode(i, i / n, 1. + 0.2 * np.sin(2 * np.pi * i / n), 0.)

                # get model_part_to (3D) from mapper
                mapper = cs_tools.CreateInstance(parameters['mapper'])
                model_part_to = mapper.Initialize(model_part_from, forward=True)

                # for model_part_from (2D): get geometry, set historical variables
                n_from = model_part_from.NumberOfNodes()
                x_from = np.zeros(n_from)
                y_from = np.zeros(n_from)
                z_from = np.zeros(n_from)
                v_s_from = np.zeros(n_from)
                v_v_from = np.zeros((n_from, 3))
                for i, node in enumerate(model_part_from.Nodes):
                    x_from[i], y_from[i], z_from[i] = node.X0, node.Y0, node.Z0

                    v_s_from[i] = x_from[i]
                    v_v_from[i] = np.array([x_from[i] * .5, x_from[i], 0.])

                    node.SetSolutionStepValue(var_s, 0, v_s_from[i])
                    node.SetSolutionStepValue(var_v, 0, tuple(v_v_from[i]))

                # map scalar and vector variables
                mapper((model_part_from, var_s), (model_part_to, var_s))
                mapper((model_part_from, var_v), (model_part_to, var_v))

                # for model_part_to (3D): get geometry, get historical variables
                n_to = model_part_to.NumberOfNodes()
                x_to = np.zeros(n_to)
                y_to = np.zeros(n_to)
                z_to = np.zeros(n_to)
                v_s_to = np.zeros(n_to)
                v_v_to = np.zeros((n_to, 3))
                for i, node in enumerate(model_part_to.Nodes):
                    x_to[i], y_to[i], z_to[i] = node.X0, node.Y0, node.Z0
                    v_s_to[i] = node.GetSolutionStepValue(var_s)
                    v_v_to[i, :] = np.array(node.GetSolutionStepValue(var_v))

                # create plot for visual check TODO doesn't work yet
                c_from = cm.jet((v_s_from - v_s_from.min()) / (v_s_from.max() - v_s_from.min()))
                c_to = cm.jet((v_s_to - v_s_from.min()) / (v_s_from.max() - v_s_from.min()))

                fig = plt.figure()

                ax_s = fig.add_subplot(121, projection='3d')
                ax_s.set_title('check geometry and scalar mapping')
                ax_s.scatter(x_from, y_from, z_from, s=50, c=c_from, depthshade=True, marker='s')
                ax_s.scatter(x_to, y_to, z_to, s=20, c=c_to, depthshade=True)

                ax_v = fig.add_subplot(122, projection='3d')
                ax_v.set_title('check vector mapping')
                ax_v.quiver(x_from, y_from, z_from, v_v_from[:, 0], v_v_from[:, 1], v_v_from[:, 2],
                            pivot='tail', arrow_length_ratio=0.1, normalize=False, length=0.2, colors='r', linewidth=2)
                ax_v.quiver(x_to, y_to, z_to, v_v_to[:, 0], v_v_to[:, 1], v_v_to[:, 2],
                            pivot='tail', arrow_length_ratio=0.1, normalize=False, length=0.2)

                for ax in [ax_s, ax_v]:
                    ax.set_xlabel('x')
                    ax.set_ylabel('y')
                    ax.set_zlabel('z')

                plt.get_current_fig_manager().window.showMaximized()
                plt.show()
                plt.close()


if __name__ == '__main__':
    KratosUnittest.main()
