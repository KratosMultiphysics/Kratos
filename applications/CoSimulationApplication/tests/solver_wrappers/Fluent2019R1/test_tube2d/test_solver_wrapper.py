import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
from copy import deepcopy
import multiprocessing
from matplotlib import pyplot
import os

class TestSolverWrapperFluent2019R1Tube2D(KratosUnittest.TestCase):
    def test_solver_wrapper(self):
        # to reset Fluent case, run file setup_case.sh

        print('Starting tests for Fluent Tube2D.')

        parameter_file_name = os.path.join(os.path.dirname(__file__),
                                           'test_solver_wrapper.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())
        par_solver_0 = parameters['solver_wrappers'][0]

        # if running from this folder
        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            par_solver_0['settings'].SetString('working_directory', 'CFD')

        displacement = vars(KM)['DISPLACEMENT']

        # test if nodes are moved to the correct position
        if True:
            # adapt Parameters, create solver
            par_solver = deepcopy(par_solver_0)
            par_solver['settings'].SetInt('flow_iterations', 5)
            solver = cs_tools.CreateInstance(par_solver)

            # give value to DISPLACEMENT variable
            mp = solver.model['beamoutside_nodes']
            for node in mp.Nodes:
                dy = 0.0005 * np.sin(2 * np.pi / 0.05 * node.X)
                node.SetSolutionStepValue(displacement, 0, [0., dy, 0.])

            # update position by iterating once in solver
            solver.Initialize()
            solver.InitializeSolutionStep()
            solver.SolveSolutionStep(solver.GetInterfaceInput())
            solver.FinalizeSolutionStep()
            solver.Finalize()

            # create solver to check new coordinates
            par_solver['settings'].SetInt('timestep_start', 1)
            solver = cs_tools.CreateInstance(par_solver)
            solver.Initialize()
            solver.Finalize()

            # check if correct displacement was given
            mp = solver.model['beamoutside_nodes']
            for node in mp.Nodes:
                y_goal = 0.005 + 0.0005 * np.sin(2 * np.pi / 0.05 * node.X)
                self.assertAlmostEqual(node.Y, y_goal, delta=1e-16)

        # test if different partitioning gives the same ModelParts
        if True:
            # create two solvers with different flow solver partitioning
            par_solver = deepcopy(par_solver_0)
            model_parts = []
            for cores in [1, multiprocessing.cpu_count()]:
                par_solver['settings'].SetInt('cores', cores)
                solver = cs_tools.CreateInstance(par_solver)
                solver.Initialize()
                solver.Finalize()
                model_parts.append(deepcopy(solver.model['beamoutside_nodes']))

            # compare Nodes in ModelParts between both solvers
            mp1, mp2 = model_parts
            for node1 in mp1.Nodes:
                node2 = mp2.GetNode(node1.Id)
                self.assertEqual(node1.X, node2.X)
                self.assertEqual(node1.Y, node2.Y)
                self.assertEqual(node1.Z, node2.Z)

        # test if same coordinates always gives same pressure & traction
        if True:
            # adapt Parameters, create solver
            par_solver = deepcopy(par_solver_0)
            par_solver['settings'].SetInt('cores', multiprocessing.cpu_count())
            par_solver['settings'].SetInt('flow_iterations', 500)
            solver = cs_tools.CreateInstance(par_solver)
            solver.Initialize()
            solver.InitializeSolutionStep()

            # change grid to position 1
            mp = solver.model['beamoutside_nodes']
            for node in mp.Nodes:
                node.Y = 0.005 + 0.0005 * np.sin(2 * np.pi / 0.05 * node.X)
            output1 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()

            # change grid to position 2
            for node in mp.Nodes:
                node.Y = 0.005 - 0.0005 * np.sin(2 * np.pi / 0.05 * node.X)
            output2 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()

            # change grid back to position 1
            for node in mp.Nodes:
                node.Y = 0.005 + 0.0005 * np.sin(2 * np.pi / 0.05 * node.X)
            output3 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()

            solver.FinalizeSolutionStep()
            solver.Finalize()

            # normalize data and compare
            a1 = output1.GetNumpyArray()
            a2 = output2.GetNumpyArray()
            a3 = output3.GetNumpyArray()

            mean = np.mean(a1)
            ref = np.abs(a1 - mean).max()

            a1n = (a1 - mean) / ref
            a2n = (a2 - mean) / ref
            a3n = (a3 - mean) / ref

            self.assertNotAlmostEqual(np.sum(np.abs(a1n - a2n)) / a1n.size, 0., delta=1e-12)
            for i in range(a1.size):
                self.assertAlmostEqual(a1n[i] - a3n[i], 0., delta=1e-12)

        # test if correct number of displacements is applied
        if True:
            # adapt Parameters, create solver
            par_solver = deepcopy(par_solver_0)
            par_solver['settings'].SetInt('flow_iterations', 5)
            solver = cs_tools.CreateInstance(par_solver)

            # give value to DISPLACEMENT variable
            mp = solver.model['beamoutside_nodes']
            for node in mp.Nodes:
                dy = 0.0001 * np.sin(2 * np.pi / 0.05 * node.X)
                node.SetSolutionStepValue(displacement, 0, [0., dy, 0.])

            # run solver for some timesteps and iterations
            solver.Initialize()
            timesteps = 3
            iterations = 4
            for i in range(timesteps):
                solver.InitializeSolutionStep()
                for j in range(iterations):
                    solver.SolveSolutionStep(solver.GetInterfaceInput())
                solver.FinalizeSolutionStep()
            solver.Finalize()

            # create solver to check coordinates at last timestep
            par_solver['settings'].SetInt('timestep_start', timesteps)
            solver = cs_tools.CreateInstance(par_solver)
            solver.Initialize()
            solver.Finalize()

            # check if displacement was applied correct number of times
            mp = solver.model['beamoutside_nodes']
            for node in mp.Nodes:
                n = timesteps * iterations
                y_goal = 0.005 + n * 0.0001 * np.sin(2 * np.pi / 0.05 * node.X)
                self.assertAlmostEqual(node.Y, y_goal, delta=1e-16)

        # test if restart option works correctly
        if True:
            # adapt Parameters, create solver
            par_solver = deepcopy(par_solver_0)
            par_solver['settings'].SetInt('cores', multiprocessing.cpu_count())
            par_solver['settings'].SetInt('flow_iterations', 500)
            solver = cs_tools.CreateInstance(par_solver)

            # give value to DISPLACEMENT variable
            mp = solver.model['beamoutside_nodes']
            for node in mp.Nodes:
                dy = 0.0002 * np.sin(2 * np.pi / 0.05 * node.X)
                node.SetSolutionStepValue(displacement, 0, [0., dy, 0.])

            # run solver for 2 timesteps
            solver.Initialize()
            for i in range(4):
                solver.InitializeSolutionStep()
                for j in range(2):
                    solver.SolveSolutionStep(solver.GetInterfaceInput())
                solver.FinalizeSolutionStep()
            solver.Finalize()

            # get data for solver without restart
            interface1 = solver.GetInterfaceOutput().deepcopy()
            data1 = interface1.GetNumpyArray().copy()

            # create solver which restarts at timestep 2
            par_solver['settings'].SetInt('timestep_start', 2)
            solver = cs_tools.CreateInstance(par_solver)

            # give value to DISPLACEMENT variable
            mp = solver.model['beamoutside_nodes']
            for node in mp.Nodes:
                dy = 0.0002 * np.sin(2 * np.pi / 0.05 * node.X)
                node.SetSolutionStepValue(displacement, 0, [0., dy, 0.])

            # run solver for 2 more timesteps
            solver.Initialize()
            for i in range(2):
                solver.InitializeSolutionStep()
                for j in range(2):
                    solver.SolveSolutionStep(solver.GetInterfaceInput())
                solver.FinalizeSolutionStep()
            solver.Finalize()

            # get data for solver with restart
            interface2 = solver.GetInterfaceOutput().deepcopy()
            data2 = interface2.GetNumpyArray().copy()

            # compare coordinates of Nodes
            mp1 = interface1.model['beamoutside_nodes']
            mp2 = interface2.model['beamoutside_nodes']
            for node1 in mp1.Nodes:
                node2 = mp2.GetNode(node1.Id)
                self.assertAlmostEqual(node1.X, node2.X, delta=1e-16)
                self.assertAlmostEqual(node1.Y, node2.Y, delta=1e-16)
                self.assertAlmostEqual(node1.Z, node2.Z, delta=1e-16)

            # normalize pressure and traction data and compare
            mean = np.mean(data1)
            ref = np.abs(data1 - mean).max()

            data1n = (data1 - mean) / ref
            data2n = (data2 - mean) / ref

            for i in range(data1n.size):
                # print(f'data1n = {data1n[i]:.1e}, data2n = {data2n[i]:.1e}, diff = {data1n[i] - data2n[i]:.1e}')
                self.assertAlmostEqual(data1n[i] - data2n[i], 0., delta=1e-14)

        print('Finishing tests for Fluent Tube2D.')
        # self.assertTrue(False)


if __name__ == '__main__':
    KratosUnittest.main()
