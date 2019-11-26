import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
import matplotlib.pyplot as plt

class TestMapperLinear(KratosUnittest.TestCase):
    def test_mapper_linear_1d(self):

        parameter_file_name = 'test_linear_1d.json'
        # parameter_file_name = 'test_nearest.json'

        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())

        # values on straight line, irregular grid spacing
        var_from = vars(KM)["TEMPERATURE"]
        model_from = cs_data_structure.Model()
        model_part_from = model_from.CreateModelPart('wall_from')
        model_part_from.AddNodalSolutionStepVariable(var_from)

        n_from = 30
        z_from, v_from = np.zeros(n_from), np.zeros(n_from)
        for i in range(n_from):
            z_from[i] = .1 + .8 * np.sqrt(i / (n_from - 1))
            v_from[i] = -3.3 * z_from[i]
            node = model_part_from.CreateNewNode(i, 0.0, 0.0, z_from[i])
            node.SetSolutionStepValue(var_from, 0, v_from[i])

        var_to = vars(KM)["PRESSURE"]
        model_to = cs_data_structure.Model()
        model_part_to = model_to.CreateModelPart('wall_to')
        model_part_to.AddNodalSolutionStepVariable(var_to)

        n_to = 25
        for i in range(n_to):
            model_part_to.CreateNewNode(i, 0.0, 0.0, .01 + i / (n_to - 1))

        mapper = cs_tools.CreateInstance(parameters['mapper'])
        mapper.Initialize(model_part_from, model_part_to)
        mapper((model_part_from, var_from), (model_part_to, var_to))

        z_to, v_to = np.zeros(n_to), np.zeros(n_to)
        for i, node in enumerate(model_part_to.Nodes):
            z_to[i] = node.Z
            v_to[i] = node.GetSolutionStepValue(var_to)

        for i in range(n_to):
            self.assertAlmostEqual(-3.3 * z_to[i], v_to[i], delta=1e-8)

        if False:
            # visualization of results
            plt.plot(z_from, v_from, 'b', label='from', linewidth=2.5)
            plt.scatter(z_from, v_from, s=20, color='b')
            plt.plot(z_to, v_to, 'r', label='to', linewidth=1.5)
            plt.scatter(z_to, v_to, s=15, color='r')
            plt.legend()
            plt.show()
            plt.close()


if __name__ == '__main__':
    KratosUnittest.main()
