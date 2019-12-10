import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
from copy import deepcopy
import multiprocessing
import os
import subprocess


class TestSolverWrapperAbaqus614(KratosUnittest.TestCase):
    def test_solver_wrapper_abaqus_614(self):
        self.test_solver_wrapper_abaqus_614_tube2d()
        self.test_solver_wrapper_abaqus_614_tube3d()

    def test_solver_wrapper_abaqus_614_tube2d(self):
        print('Starting tests for Abaqus Tube2D.')

        parameter_file_name = os.path.join(os.path.dirname(__file__),
                                           'test_614_tube2D', 'test_solver_wrapper.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())
        par_solver_0 = parameters['solver_wrappers'][0]

        # if running from this folder
        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            par_solver_0['settings'].SetString('working_directory', 'test_614_tube2D/CSM')
            par_solver_0['settings'].SetString('input_file', 'test_614_tube2D/Base.inp')

        # Create hostfile for Abaqus
        os.system("cd test_614_tube2D/CSM; ./makeHostFile.sh")

        par_solver = deepcopy(par_solver_0)

        pressure = vars(KM)['PRESSURE']
        traction = vars(KM)['TRACTION']

        p = 10000
        shear_x = 0
        shear_y = 0
        shear_z = 0

        # Create solver0
        if True:
            # Create the solver (__init__)
            print("Creating an AbaqusSolver")
            AbaqusSolver0 = cs_tools.CreateInstance(par_solver_0)
            print("AbaqusSolver0 created")

        # Test start and restart
        if True:
            mp = AbaqusSolver0.model['BEAMINSIDEMOVING_load_points']
            for node in mp.Nodes:
                # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                # print(node.Y)
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            AbaqusSolver0.Initialize()

            # Step 1, Coupling 1
            AbaqusSolver0.InitializeSolutionStep()
            output1_1 = AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput())
            os.system("cp -r test_614_tube2D/CSM/CSM_Time1.odb test_614_tube2D/CSM/CSM_Time1_Iter1.odb")
            # Step 1, Coupling 2
            output1_2 = AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput()).deepcopy()
            AbaqusSolver0.FinalizeSolutionStep()

            # Compare output, as input hasn't changed these should be the same
            # normalize data and compare
            a1 = output1_1.GetNumpyArray()
            a2 = output1_2.GetNumpyArray()

            mean = np.mean(a1)
            ref = np.abs(a1 - mean).max()

            a1n = (a1 - mean) / ref
            a2n = (a2 - mean) / ref

            for i in range(a1.size):
                self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)

            # Step 2 and 3
            for i in range(2):
                AbaqusSolver0.InitializeSolutionStep()
                AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput())
                AbaqusSolver0.FinalizeSolutionStep()
            # Step 4
            AbaqusSolver0.InitializeSolutionStep()
            output_single_run = AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput()).deepcopy()
            AbaqusSolver0.FinalizeSolutionStep()
            AbaqusSolver0.Finalize()

            os.system("cp test_614_tube2D/CSM/CSM_Time4Surface0Output.dat test_614_tube2D/CSM/CSM_Time4Surface0Output_Single.dat")

            # With restart
            # create solver which restarts at timestep 2
            par_solver['settings'].SetInt('timestep_start', 2)
            AbaqusSolver1 = cs_tools.CreateInstance(par_solver)
            mp = AbaqusSolver1.model['BEAMINSIDEMOVING_load_points']
            for node in mp.Nodes:
                # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                # print(node.Y)
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            AbaqusSolver1.Initialize()

            for i in range(2):
                AbaqusSolver1.InitializeSolutionStep()
                output_restart = AbaqusSolver1.SolveSolutionStep(AbaqusSolver1.GetInterfaceInput()).deepcopy()
                AbaqusSolver1.FinalizeSolutionStep()
            AbaqusSolver1.Finalize()

            # Compare output, as input hasn't changed these should be the same
            # normalize data and compare
            a1 = output_single_run.GetNumpyArray()
            a2 = output_restart.GetNumpyArray()

            mean = np.mean(a1)
            ref = np.abs(a1 - mean).max()

            a1n = (a1 - mean) / ref
            a2n = (a2 - mean) / ref

            for i in range(a1.size):
                # print(f"{a1[i]} ?= {a2[i]}")
                self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)

        if True:
            # Test whether using 4 cpus gives the same results as using a single one.
            par_solver = deepcopy(par_solver_0)
            par_solver["settings"].SetInt("cores", 4)
            AbaqusSolver2 = cs_tools.CreateInstance(par_solver)
            mp = AbaqusSolver2.model['BEAMINSIDEMOVING_load_points']
            for node in mp.Nodes:
                # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                # print(node.Y)
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

            AbaqusSolver2.Initialize()
            for i in range(4):
                AbaqusSolver2.InitializeSolutionStep()
                output_4cores = AbaqusSolver2.SolveSolutionStep(AbaqusSolver2.GetInterfaceInput()).deepcopy()
                AbaqusSolver2.FinalizeSolutionStep()
            AbaqusSolver2.Finalize()

            # Compare output, as input hasn't changed these should be the same
            # normalize data and compare
            a4 = output_4cores.GetNumpyArray()
            print(a4.shape)
            a4n = (a4 - mean)/ref

            for i in range(a1.size):
                self.assertAlmostEqual(a2n[i] - a4n[i], 0., delta=1e-12)
                self.assertAlmostEqual(a1n[i] - a4n[i], 0., delta=1e-12)

        if True:
            # Test whether shear is also applied
            shear_y = p
            mp = AbaqusSolver2.model['BEAMINSIDEMOVING_load_points']
            for node in mp.Nodes:
                # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                # print(node.Y)
                node.SetSolutionStepValue(pressure, 0, p)
                node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])
            AbaqusSolver2.Initialize()
            for i in range(4):
                AbaqusSolver2.InitializeSolutionStep()
                output_shear = AbaqusSolver2.SolveSolutionStep(AbaqusSolver2.GetInterfaceInput()).deepcopy()
                AbaqusSolver2.FinalizeSolutionStep()
            AbaqusSolver2.Finalize()

            a5 = output_shear.GetNumpyArray()

            mean_disp_y_no_shear = 0
            mean_disp_y_shear = 0

            for i in range(1,a1.size,3):
                mean_disp_y_no_shear += a2[i]
                mean_disp_y_shear += a5[i]

            mean_disp_y_no_shear /= (a1.size/3)
            mean_disp_y_shear /= (a1.size/3)

            print(f"Mean y-displacement without shear = {mean_disp_y_no_shear} m")
            print(f"Mean y-displacement with shear = {mean_disp_y_shear} m")

            self.assertNotAlmostEqual(mean_disp_y_no_shear - mean_disp_y_shear, 0., delta=1e-12)

    def test_solver_wrapper_abaqus_614_tube3d(self):
        print('Starting tests for Abaqus Tube3D.')
        #Axial direction is the x-direction

        parameter_file_name = os.path.join(os.path.dirname(__file__),
                                           'test_614_tube3D', 'test_solver_wrapper.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())
        par_solver_0 = parameters['solver_wrappers'][0]

        # if running from this folder
        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            par_solver_0['settings'].SetString('working_directory', 'test_614_tube3D/CSM')
            par_solver_0['settings'].SetString('input_file', 'test_614_tube3D/Base.inp')

        # Create hostfile for Abaqus
        os.system("cd test_614_tube3D/CSM; ./makeHostFile.sh")

        par_solver = deepcopy(par_solver_0)

        pressure = vars(KM)['PRESSURE']
        traction = vars(KM)['TRACTION']

        p = 10000
        shear_x = 0
        shear_y = 0
        shear_z = 0

        # Create solver0
        if True:
            # Create the solver (__init__)
            print("Creating an AbaqusSolver")
            AbaqusSolver0 = cs_tools.CreateInstance(par_solver_0)
            print("AbaqusSolver0 created")

            # Test start and restart
            if True:
                mp = AbaqusSolver0.model['WALLOUTSIDE_load_points']
                for node in mp.Nodes:
                    # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                    # print(node.Y)
                    node.SetSolutionStepValue(pressure, 0, p)
                    node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

                AbaqusSolver0.Initialize()

                # Step 1, Coupling 1
                AbaqusSolver0.InitializeSolutionStep()
                output1_1 = AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput())
                os.system("cp -r test_614_tube3D/CSM/CSM_Time1.odb test_614_tube3D/CSM/CSM_Time1_Iter1.odb")
                # Step 1, Coupling 2
                output1_2 = AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput()).deepcopy()
                AbaqusSolver0.FinalizeSolutionStep()

                # Compare output, as input hasn't changed these should be the same
                # normalize data and compare
                a1 = output1_1.GetNumpyArray()
                a2 = output1_2.GetNumpyArray()

                mean = np.mean(a1)
                ref = np.abs(a1 - mean).max()

                a1n = (a1 - mean) / ref
                a2n = (a2 - mean) / ref

                for i in range(a1.size):
                    self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)

                # Step 2 and 3
                for i in range(2):
                    AbaqusSolver0.InitializeSolutionStep()
                    AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput())
                    AbaqusSolver0.FinalizeSolutionStep()
                # Step 4
                AbaqusSolver0.InitializeSolutionStep()
                output_single_run = AbaqusSolver0.SolveSolutionStep(AbaqusSolver0.GetInterfaceInput()).deepcopy()
                AbaqusSolver0.FinalizeSolutionStep()
                AbaqusSolver0.Finalize()

                os.system(
                    "cp test_614_tube3D/CSM/CSM_Time4Surface0Output.dat "
                    "test_614_tube3D/CSM/CSM_Time4Surface0Output_Single.dat")

                # With restart
                # create solver which restarts at timestep 2
                par_solver['settings'].SetInt('timestep_start', 2)
                AbaqusSolver1 = cs_tools.CreateInstance(par_solver)
                mp = AbaqusSolver1.model['WALLOUTSIDE_load_points']
                for node in mp.Nodes:
                    # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                    # print(node.Y)
                    node.SetSolutionStepValue(pressure, 0, p)
                    node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

                AbaqusSolver1.Initialize()

                for i in range(2):
                    AbaqusSolver1.InitializeSolutionStep()
                    output_restart = AbaqusSolver1.SolveSolutionStep(AbaqusSolver1.GetInterfaceInput()).deepcopy()
                    AbaqusSolver1.FinalizeSolutionStep()
                AbaqusSolver1.Finalize()

                # Compare output, as input hasn't changed these should be the same
                # normalize data and compare
                a1 = output_single_run.GetNumpyArray()
                a2 = output_restart.GetNumpyArray()

                mean = np.mean(a1)
                ref = np.abs(a1 - mean).max()

                a1n = (a1 - mean) / ref
                a2n = (a2 - mean) / ref

                for i in range(a1.size):
                    # print(f"{a1[i]} ?= {a2[i]}")
                    self.assertAlmostEqual(a1n[i] - a2n[i], 0., delta=1e-12)

            if True:
                # Test whether using 4 cpus gives the same results as using a single one.
                par_solver = deepcopy(par_solver_0)
                par_solver["settings"].SetInt("cores", 4)
                AbaqusSolver2 = cs_tools.CreateInstance(par_solver)
                mp = AbaqusSolver2.model['WALLOUTSIDE_load_points']
                for node in mp.Nodes:
                    # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                    # print(node.Y)
                    node.SetSolutionStepValue(pressure, 0, p)
                    node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])

                AbaqusSolver2.Initialize()
                for i in range(4):
                    AbaqusSolver2.InitializeSolutionStep()
                    output_4cores = AbaqusSolver2.SolveSolutionStep(AbaqusSolver2.GetInterfaceInput()).deepcopy()
                    AbaqusSolver2.FinalizeSolutionStep()
                AbaqusSolver2.Finalize()

                # Compare output, as input hasn't changed these should be the same
                # normalize data and compare
                a4 = output_4cores.GetNumpyArray()
                print(a4.shape)
                a4n = (a4 - mean) / ref

                for i in range(a1.size):
                    self.assertAlmostEqual(a2n[i] - a4n[i], 0., delta=1e-12)
                    self.assertAlmostEqual(a1n[i] - a4n[i], 0., delta=1e-12)

            if True:
                # Test whether shear is also applied (x is the axial direction)
                shear_x = p
                mp = AbaqusSolver2.model['WALLOUTSIDE_load_points']
                for node in mp.Nodes:
                    # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
                    # print(node.Y)
                    node.SetSolutionStepValue(pressure, 0, p)
                    node.SetSolutionStepValue(traction, 0, [shear_x, shear_y, shear_z])
                AbaqusSolver2.Initialize()
                for i in range(4):
                    AbaqusSolver2.InitializeSolutionStep()
                    output_shear = AbaqusSolver2.SolveSolutionStep(AbaqusSolver2.GetInterfaceInput()).deepcopy()
                    AbaqusSolver2.FinalizeSolutionStep()
                AbaqusSolver2.Finalize()

                a5 = output_shear.GetNumpyArray()

                mean_disp_x_no_shear = 0
                mean_disp_x_shear = 0

                for i in range(0, a1.size, 3):
                    mean_disp_x_no_shear += a2[i]
                    mean_disp_x_shear += a5[i]

                mean_disp_x_no_shear /= (a1.size / 3)
                mean_disp_x_shear /= (a1.size / 3)

                print(f"Mean x-displacement without shear = {mean_disp_x_no_shear} m")
                print(f"Mean x-displacement with shear = {mean_disp_x_shear} m")

                self.assertNotAlmostEqual(mean_disp_x_no_shear - mean_disp_x_shear, 0., delta=1e-12)



if __name__ == '__main__':
    KratosUnittest.main()
