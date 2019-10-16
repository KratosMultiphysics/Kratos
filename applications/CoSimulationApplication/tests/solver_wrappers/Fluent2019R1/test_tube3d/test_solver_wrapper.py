import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
from copy import deepcopy
import multiprocessing
import os

class TestSolverWrapperFluent2019R1Tube3D(KratosUnittest.TestCase):
    def test_solver_wrapper(self):
        # to reset Fluent case, run file setup_case.sh

        print('Starting tests for Fluent Tube3D.')

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

        def convert(dr, X, Y):
            theta = np.arctan2(Y, X)
            dx = dr * np.cos(theta)
            dy = dr * np.sin(theta)
            return dx, dy

        # test if nodes are moved to the correct position
        if True:
            # adapt Parameters, create solver
            par_solver = deepcopy(par_solver_0)
            par_solver['settings'].SetInt('flow_iterations', 5)
            solver = cs_tools.CreateInstance(par_solver)

            # give value to DISPLACEMENT variable
            mp = solver.model['wall_nodes']
            for node in mp.Nodes:
                dr = 0.0005 * np.sin(2 * np.pi / 0.05 * node.X)
                dy, dz = convert(dr, node.Y, node.Z)
                node.SetSolutionStepValue(displacement, 0, [0., dy, dz])
            mp0 = solver.GetInterfaceInput().deepcopy().model['wall_nodes']
            # mp0 = deepcopy(mp)

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
            mp = solver.model['wall_nodes']
            for node in mp.Nodes:
                node0 = mp0.GetNode(node.Id)
                disp = node0.GetSolutionStepValue(displacement)
                y_goal = node0.Y0 + disp[1]
                z_goal = node0.Z0 + disp[2]
                self.assertAlmostEqual(node.Y, y_goal, delta=1e-16)
                self.assertAlmostEqual(node.Z, z_goal, delta=1e-16)

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
                model_parts.append(deepcopy(solver.model['wall_nodes']))

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
            mp = solver.model['wall_nodes']
            for node in mp.Nodes:
                dr = 0.0005 * np.sin(2 * np.pi / 0.05 * node.X0)
                dy, dz = convert(dr, node.Y0, node.Z0)
                node.Y = node.Y0 + dy
                node.Z = node.Z0 + dz
            output1 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()

            # change grid to position 2
            for node in mp.Nodes:
                dr = 0.0005 * np.sin(2 * np.pi / 0.05 * node.X0)
                dy, dz = convert(dr, node.Y0, node.Z0)
                node.Y = node.Y0 - dy
                node.Z = node.Z0 - dz
            output2 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()

            # change grid back to position 1
            for node in mp.Nodes:
                dr = 0.0005 * np.sin(2 * np.pi / 0.05 * node.X0)
                dy, dz = convert(dr, node.Y0, node.Z0)
                node.Y = node.Y0 + dy
                node.Z = node.Z0 + dz
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
            mp = solver.model['wall_nodes']
            for node in mp.Nodes:
                dr = 0.0001 * np.sin(2 * np.pi / 0.05 * node.X0)
                dy, dz = convert(dr, node.Y0, node.Z0)
                node.SetSolutionStepValue(displacement, 0, [0., dy, dz])
            mp0 = solver.GetInterfaceInput().deepcopy().model['wall_nodes']

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
            mp = solver.model['wall_nodes']
            for node in mp.Nodes:
                n = timesteps * iterations
                node0 = mp0.GetNode(node.Id)
                disp = node0.GetSolutionStepValue(displacement)
                y_goal = node0.Y0 + n * disp[1]
                z_goal = node0.Z0 + n * disp[2]
                self.assertAlmostEqual(node.Y, y_goal, delta=1e-16)
                self.assertAlmostEqual(node.Z, z_goal, delta=1e-16)

        # test if restart option works correctly
        if True:
            # adapt Parameters, create solver
            par_solver = deepcopy(par_solver_0)
            par_solver['settings'].SetInt('cores', multiprocessing.cpu_count())
            par_solver['settings'].SetInt('flow_iterations', 500)
            solver = cs_tools.CreateInstance(par_solver)

            # give value to DISPLACEMENT variable
            mp = solver.model['wall_nodes']
            for node in mp.Nodes:
                dr = 0.0002 * np.sin(2 * np.pi / 0.05 * node.X0)
                dy, dz = convert(dr, node.Y0, node.Z0)
                node.SetSolutionStepValue(displacement, 0, [0., dy, dz])
            mp0 = solver.GetInterfaceInput().deepcopy().model['wall_nodes']

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
            mp = solver.model['wall_nodes']
            for node in mp.Nodes:
                node0 = mp0.GetNode(node.Id)
                node.SetSolutionStepValue(displacement, 0,
                                          node0.GetSolutionStepValue(displacement))

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
            mp1 = interface1.model['wall_nodes']
            mp2 = interface2.model['wall_nodes']
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
                self.assertAlmostEqual(data1n[i] - data2n[i], 0., delta=1e-14)

        print('Finishing tests for Fluent Tube3D.')


if __name__ == '__main__':
    KratosUnittest.main()
