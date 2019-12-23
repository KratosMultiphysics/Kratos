import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
import os
from copy import deepcopy
import matplotlib.pyplot as plt


class TestMapperNearest(KratosUnittest.TestCase):
    def test_mapper_nearest(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_nearest.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())
        par_mapper_0 = parameters['mapper']

        gui = True

        # 1D: square-root grid + linear function
        args_from, args_to, fun = case_1d(cs_data_structure)

        par_mapper = deepcopy(par_mapper_0)
        par_mapper['settings'].SetArray('directions', ['Z'])
        mapper = cs_tools.CreateInstance(par_mapper)
        mapper.Initialize(args_from[0], args_to[0])
        mapper(args_from, args_to)

        self.assertTrue(case_1d_check(args_to, fun))

        if gui:
            case_1d_plot(args_from, args_to, fun)

        # 3D:




        # *** old code
        print('\nHERE STARTS OLD CODE\n')

        if False:
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

            mapper = cs_tools.CreateInstance(parameters['mapper_1d'])
            mapper.Initialize(model_part_from, model_part_to)
            mapper((model_part_from, var_from), (model_part_to, var_to))

            values = []
            for i_to, node in enumerate(model_part_to.Nodes):
                values.append(node.GetSolutionStepValue(var_to))

            self.assertListEqual(values, [0., 9., 49., 100.])


            # test 3D problem with perturbed grid, for scalar and vector variables
            var1_from = vars(KM)["TEMPERATURE"]
            var3_from = vars(KM)["DISPLACEMENT"]
            model_from = cs_data_structure.Model()
            model_part_from = model_from.CreateModelPart('wall_from')
            model_part_from.AddNodalSolutionStepVariable(var1_from)
            model_part_from.AddNodalSolutionStepVariable(var3_from)

            n = 10  # n = 10 --> 7e3 points; n = 20 --> 5e4 points
            x = np.linspace(0, n, n + 1)
            y = np.linspace(0, 2 * n, 2 * n + 1)
            z = np.linspace(0, 3 * n, 3 * n + 1)
            X, Y, Z = np.meshgrid(x, y, z)
            X, Y, Z = X.flatten(), Y.flatten(), Z.flatten()
            print(f'#DoFs = {X.size}')

            for i in range(X.size):
                node = model_part_from.CreateNewNode(i, X[i], Y[i], Z[i])
                node.SetSolutionStepValue(var1_from, 0, np.random.rand())
                node.SetSolutionStepValue(var3_from, 0, np.random.rand(3).tolist())

            dX = (np.random.rand(X.size) - .5) * .49
            dY = (np.random.rand(X.size) - .5) * .49
            dZ = (np.random.rand(X.size) - .5) * .49

            Xb, Yb, Zb = X + dX, Y + dY, Z + dZ

            var1_to = vars(KM)["PRESSURE"]
            var3_to = vars(KM)["FORCE"]
            model_to = cs_data_structure.Model()
            model_part_to = model_to.CreateModelPart('wall_to')
            model_part_to.AddNodalSolutionStepVariable(var1_to)
            model_part_to.AddNodalSolutionStepVariable(var3_to)

            for i in range(Xb.size):
                model_part_to.CreateNewNode(i, Xb[i], Yb[i], Zb[i])

            mapper1 = cs_tools.CreateInstance(parameters['mapper_3d'])
            mapper1.Initialize(model_part_from, model_part_to)
            mapper1((model_part_from, var1_from), (model_part_to, var1_to))

            mapper3 = cs_tools.CreateInstance(parameters['mapper_3d'])
            mapper3.Initialize(model_part_from, model_part_to)
            mapper3((model_part_from, var3_from), (model_part_to, var3_to))

            for node_from, node_to in zip(model_part_from.Nodes, model_part_to.Nodes):
                self.assertEqual(node_from.GetSolutionStepValue(var1_from),
                                 node_to.GetSolutionStepValue(var1_to))
                self.assertListEqual(node_from.GetSolutionStepValue(var3_from),
                                 node_to.GetSolutionStepValue(var3_to))


def case_1d(cs_data_structure):
    # 1D: square-root grid + linear function
    def fun(z):
        return z / 10

    var_from = vars(KM)["TEMPERATURE"]
    model_from = cs_data_structure.Model()
    model_part_from = model_from.CreateModelPart('wall_from')
    model_part_from.AddNodalSolutionStepVariable(var_from)
    n_from = 14
    z_from = np.linspace(0, 10, n_from) ** .5
    v_from = fun(z_from)
    for i in range(n_from):
        node = model_part_from.CreateNewNode(i, 0., 0., z_from[i])
        node.SetSolutionStepValue(var_from, 0, v_from[i])

    var_to = vars(KM)["PRESSURE"]
    model_to = cs_data_structure.Model()
    model_part_to = model_to.CreateModelPart('wall_to')
    model_part_to.AddNodalSolutionStepVariable(var_to)
    n_to = 5
    z_to = np.linspace(0, 10, n_to) ** .5
    v_to = np.zeros_like(z_to)
    for i in range(n_to):
        model_part_to.CreateNewNode(i, 0., 0., z_to[i])

    return (model_part_from, var_from), (model_part_to, var_to), fun

def case_1d_check(args_to, fun):
    model_part_to, var_to = args_to

    n_to = model_part_to.NumberOfNodes()
    z_to = np.zeros(n_to)
    v_to = np.zeros(n_to)
    for i_to, node in enumerate(model_part_to.Nodes):
        z_to[i_to] = node.Z
        v_to[i_to] = node.GetSolutionStepValue(var_to)

    v_error = np.abs(v_to - fun(z_to))
    criterion = (v_error < 0.02)
    return criterion.all()

def case_1d_plot(args_from, args_to, fun):
    model_part_from, var_from = args_from
    model_part_to, var_to = args_to

    n_from = model_part_from.NumberOfNodes()
    z_from = np.zeros(n_from)
    v_from = np.zeros(n_from)
    for i, node in enumerate(model_part_from.Nodes):
        z_from[i] = node.Z
        v_from[i] = node.GetSolutionStepValue(var_from)

    n_to = model_part_to.NumberOfNodes()
    z_to = np.zeros(n_to)
    v_to = np.zeros(n_to)
    for i, node in enumerate(model_part_to.Nodes):
        z_to[i] = node.Z
        v_to[i] = node.GetSolutionStepValue(var_to)
    v_to_fun = fun(z_to)

    _, ax = plt.subplots(ncols=2, sharex=True, figsize=(15, 6))

    ax[0].plot(z_from, v_from, label='from', marker='o')
    ax[0].plot(z_to, v_to, label='to', marker='o')

    ax[1].plot(z_to, np.abs(v_to_fun - v_to), label='error', marker='o')

    for a in ax:
        a.legend()

    plt.tight_layout()
    plt.show()
    plt.close()

def case_3d_sphere(cs_data_structure):


    pass



if __name__ == '__main__':
    KratosUnittest.main()
