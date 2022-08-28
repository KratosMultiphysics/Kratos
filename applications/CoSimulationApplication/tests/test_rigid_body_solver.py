import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC

from KratosMultiphysics import kratos_utilities
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.sdof.sdof_solver import SDoFSolver
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.sdof.sdof_solver_wrapper import Create as CreateSDofSolverWrapper

from KratosMultiphysics.CoSimulationApplication.solver_wrappers.rigid_body.rigid_body_solver import RigidBodySolver
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.rigid_body import rigid_body_input_check as input_check

import os
import numpy as np

class TestRigidBodySolver(KratosUnittest.TestCase):
    def setUp(self):
        self.model = KM.Model()
        self.default_parameters = KM.Parameters('''{
            "problem_data": {
                "problem_name": "RigidBodyStandalone",
                "start_time": 0.0,
                "end_time": 1.0,
                "echo_level" : 0
            },
            "solver_settings": {
                "domain_size": 3,
                "echo_level": 0,
                "buffer_size": 3,
                "model_import_settings": {
                    "input_type": "none",
                    "input_filename": "Main",
                    "restart_load_file_label": "11.0",
                    "input_output_path": "restart"
                },
                "time_integration_parameters": {
                    "rho_inf": 1,
                    "time_step": 0.01
                },
                "active_dofs": [
                    {
                        "dof": "displacement_x"
                    }
                ]
            },
            "output_processes": [],
            "processes": {
                "gravity": [],
                "initial_conditions_process_list": [],
                "boundary_conditions_process_list": [],
                "auxiliar_process_list": []
            }
        }''')


    def tearDown(self):
        kratos_utilities.DeleteFileIfExisting("result.dat")

    def __CompareResults(self, reference, result):
        pass



    def test_InitializeSolutionVariables(self):
        # Read from solver_settings
        solver_settings = self.default_parameters["solver_settings"]
        domain_size = solver_settings["domain_size"].GetInt()
        echo_level = solver_settings["echo_level"].GetInt()
        rho_inf = solver_settings["time_integration_parameters"]["rho_inf"].GetDouble()
        time_step = solver_settings["time_integration_parameters"]["time_step"].GetDouble()
        # Read from model
        simulation = RigidBodySolver(self.model, self.default_parameters)
        
        self.assertEqual(simulation.domain_size, domain_size)
        self.assertEqual(simulation.echo_level, echo_level)
        self.assertEqual(simulation.rho_inf, rho_inf)
        self.assertEqual(simulation.delta_t, time_step)

    def test_InitializeDofsVariables1(self):
        Parameters = KM.Parameters('''{
            "solver_settings":{
                "active_dofs": [
                    {
                        "dof": "displacement_x",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 40.0,
                            "stiffness": 400.0,
                            "damping": 4.0
                        }
                    },
                    {
                        "dof": "displacement_y",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 50.0,
                            "stiffness": 500.0,
                            "damping": 5.0
                        }
                    },
                    {
                        "dof": "displacement_z",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 100.0,
                            "stiffness": 1000.0,
                            "damping": 10.0
                        }
                    },
                    {
                        "dof": "rotation_x",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 2.0,
                            "stiffness": 200.0,
                            "damping": 2.0
                        }
                    },
                    {
                        "dof": "rotation_y",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 10.0,
                            "stiffness": 1000.0,
                            "damping": 10.0
                        }
                    },
                    {
                        "dof": "rotation_z",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 5.0,
                            "stiffness": 500.0,
                            "damping": 5.0
                        }
                    }
                ]
            }
        }''')
        Parameters.RecursivelyAddMissingParameters(self.default_parameters)

        simulation = RigidBodySolver(self.model, Parameters)
        
        Constrains = {'displacement_x': False, 'displacement_y': False, 'displacement_z': False, 'rotation_x': False, 'rotation_y': False, 'rotation_z': False}
        M = [[ 40.,  0. ,   0. ,  0.,   0.,   0.],
            [  0. , 50. ,   0. ,  0. ,  0. ,  0.],
            [  0. ,  0. , 100. ,  0. ,  0. ,  0.],
            [  0. ,  0. ,   0. ,  2. ,  0. ,  0.],
            [  0. ,  0. ,   0. ,  0. , 10. ,  0.],
            [  0. ,  0. ,   0. ,  0. ,  0. ,  5.]]
        C = [[ 4.,  0.,  0.,  0.,  0.,  0.],
             [ 0.,  5.,  0.,  0.,  0.,  0.],
             [ 0.,  0., 10.,  0.,  0.,  0.],
             [ 0.,  0.,  0.,  2.,  0.,  0.],
             [ 0.,  0.,  0.,  0., 10.,  0.],
             [ 0.,  0.,  0.,  0.,  0.,  5.]]
        K = [[ 400.,    0.,    0.,    0.,    0.,    0.],
             [   0.,  500.,    0.,    0.,    0.,    0.],
             [   0.,    0., 1000.,    0.,    0.,    0.],
             [   0.,    0.,    0.,  200.,    0.,    0.],
             [   0.,    0.,    0.,    0., 1000.,    0.],
             [   0.,    0.,    0.,    0.,    0.,  500.]]

        for i, dof in enumerate(['displacement_x', 'displacement_y', 'displacement_z', 'rotation_x', 'rotation_y', 'rotation_z']):
            self.assertEqual(simulation.is_constrained[dof], Constrains[dof])
            for j in range(6):
                self.assertEqual(simulation.M[i][j], M[i][j])
                self.assertEqual(simulation.C[i][j], C[i][j])
                self.assertEqual(simulation.K[i][j], K[i][j])


    def test_InitializeDofsVariables2(self):
        Parameters = KM.Parameters('''{
            "solver_settings":{
                "active_dofs": [
                    {
                        "dof": "displacement_x",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 40.0,
                            "stiffness": 400.0,
                            "damping": 4.0
                        }
                    },
                    {
                        "dof": "displacement_y",
                        "constrained": true,
                        "system_parameters": {
                            "mass": 50.0,
                            "stiffness": 500.0,
                            "damping": 5.0
                        }
                    },
                    {
                        "dof": "displacement_z",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 100.0,
                            "stiffness": 1000.0,
                            "damping": 10.0
                        }
                    },
                    {
                        "dof": "rotation_x",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 2.0,
                            "stiffness": 200.0,
                            "damping": 2.0
                        }
                    },
                    {
                        "dof": "rotation_y",
                        "constrained": true,
                        "system_parameters": {
                            "mass": 10.0,
                            "stiffness": 1000.0,
                            "damping": 10.0
                        }
                    }
                ]
            }
        }''')
        Parameters.RecursivelyAddMissingParameters(self.default_parameters)

        simulation = RigidBodySolver(self.model, Parameters)
        
        Constrains = {'displacement_x': False, 'displacement_y': True, 'displacement_z': False, 'rotation_x': False, 'rotation_y': True, 'rotation_z': False}
        M = [[ 40.,  0. ,   0. ,  0.,   0.,   0.],
            [  0. , 50. ,   0. ,  0. ,  0. ,  0.],
            [  0. ,  0. , 100. ,  0. ,  0. ,  0.],
            [  0. ,  0. ,   0. ,  2. ,  0. ,  0.],
            [  0. ,  0. ,   0. ,  0. , 10. ,  0.],
            [  0. ,  0. ,   0. ,  0. ,  0. ,  1.]]
        C = [[ 4.,  0.,  0.,  0.,  0.,  0.],
             [ 0.,  5.,  0.,  0.,  0.,  0.],
             [ 0.,  0., 10.,  0.,  0.,  0.],
             [ 0.,  0.,  0.,  2.,  0.,  0.],
             [ 0.,  0.,  0.,  0., 10.,  0.],
             [ 0.,  0.,  0.,  0.,  0.,  0.]]
        K = [[ 400.,    0.,    0.,    0.,    0.,    0.],
             [   0.,  500.,    0.,    0.,    0.,    0.],
             [   0.,    0., 1000.,    0.,    0.,    0.],
             [   0.,    0.,    0.,  200.,    0.,    0.],
             [   0.,    0.,    0.,    0., 1000.,    0.],
             [   0.,    0.,    0.,    0.,    0.,  1.]]

        for i, dof in enumerate(['displacement_x', 'displacement_y', 'displacement_z', 'rotation_x', 'rotation_y', 'rotation_z']):
            self.assertEqual(simulation.is_constrained[dof], Constrains[dof])
            for j in range(6):
                self.assertEqual(simulation.M[i][j], M[i][j])
                self.assertEqual(simulation.C[i][j], C[i][j])
                self.assertEqual(simulation.K[i][j], K[i][j])

    def test_InitializeGeneralizedAlphaParameters1(self):
        simulation = RigidBodySolver(self.model, self.default_parameters)

        self.assertAlmostEqual(simulation.alpha_f, 0.5)
        self.assertAlmostEqual(simulation.alpha_m, 0.5)
        self.assertAlmostEqual(simulation.beta, 0.25)
        self.assertAlmostEqual(simulation.gamma, 0.5)
        self.assertAlmostEqual(simulation.a1h, 20000.0)
        self.assertAlmostEqual(simulation.a2h, 100.0)
        self.assertAlmostEqual(simulation.a3h, 0.5)
        self.assertAlmostEqual(simulation.a1m, 20000.0)
        self.assertAlmostEqual(simulation.a2m, 200.0)
        self.assertAlmostEqual(simulation.a3m, 0.0)
        self.assertAlmostEqual(simulation.a1b, 100.0)
        self.assertAlmostEqual(simulation.a2b, 0.0)
        self.assertAlmostEqual(simulation.a3b, 0.0)
        self.assertAlmostEqual(simulation.a1k, -0.5)
        self.assertAlmostEqual(simulation.a1v, 200.0)
        self.assertAlmostEqual(simulation.a2v, -1.0)
        self.assertAlmostEqual(simulation.a3v, 0.0)
        self.assertAlmostEqual(simulation.a1a, 40000.0)
        self.assertAlmostEqual(simulation.a2a, -400.0)
        self.assertAlmostEqual(simulation.a3a, -1.0)

        LHS= [[20000.5,     0. ,     0. ,     0. ,     0. ,     0. ],
              [    0. , 20000.5,     0. ,     0. ,     0. ,     0. ],
              [    0. ,     0. , 20000.5,     0. ,     0. ,     0. ],
              [    0. ,     0. ,     0. , 20000.5,     0. ,     0. ],
              [    0. ,     0. ,     0. ,     0. , 20000.5,     0. ],
              [    0. ,     0. ,     0. ,     0. ,     0. , 20000.5]]


        for i in range(6):
            for j in range(6):
                self.assertAlmostEqual(simulation.LHS[i][j], LHS[i][j])

    def test_InitializeGeneralizedAlphaParameters2(self):
        Parameters = KM.Parameters('''{
            "solver_settings": {
                "time_integration_parameters": {
                        "rho_inf": 0.16,
                        "time_step": 0.01
                },
                "active_dofs": [
                    {
                        "dof": "displacement_x",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 40.0,
                            "stiffness": 400.0,
                            "damping": 4.0
                        }
                    },
                    {
                        "dof": "displacement_y",
                        "constrained": true,
                        "system_parameters": {
                            "mass": 50.0,
                            "stiffness": 500.0,
                            "damping": 5.0
                        }
                    },
                    {
                        "dof": "displacement_z",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 100.0,
                            "stiffness": 1000.0,
                            "damping": 10.0
                        }
                    },
                    {
                        "dof": "rotation_x",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 2.0,
                            "stiffness": 200.0,
                            "damping": 2.0
                        }
                    },
                    {
                        "dof": "rotation_y",
                        "constrained": true,
                        "system_parameters": {
                            "mass": 10.0,
                            "stiffness": 1000.0,
                            "damping": 10.0
                        }
                    }
                ]
            }
        }''')
        Parameters.RecursivelyAddMissingParameters(self.default_parameters)
        simulation = RigidBodySolver(self.model, Parameters)

        self.assertAlmostEqual(simulation.alpha_f, 0.13793103448275865)
        self.assertAlmostEqual(simulation.alpha_m, -0.5862068965517241)
        self.assertAlmostEqual(simulation.beta, 0.7431629013079668)
        self.assertAlmostEqual(simulation.gamma, 1.2241379310344829)
        self.assertAlmostEqual(simulation.a1h, 21344.0)
        self.assertAlmostEqual(simulation.a2h, 142.0)
        self.assertAlmostEqual(simulation.a3h, 0.8620689655172413)
        self.assertAlmostEqual(simulation.a1m, 21344.0)
        self.assertAlmostEqual(simulation.a2m, 213.44)
        self.assertAlmostEqual(simulation.a3m, 0.0671999999999999)
        self.assertAlmostEqual(simulation.a1b, 142.0)
        self.assertAlmostEqual(simulation.a2b, 0.41999999999999993)
        self.assertAlmostEqual(simulation.a3b, -0.0015206896551724137)
        self.assertAlmostEqual(simulation.a1k, -0.13793103448275865)
        self.assertAlmostEqual(simulation.a1v, 164.72)
        self.assertAlmostEqual(simulation.a2v, -0.6472)
        self.assertAlmostEqual(simulation.a3v, 0.001764)
        self.assertAlmostEqual(simulation.a1a, 13455.999999999998)
        self.assertAlmostEqual(simulation.a2a, -134.55999999999997)
        self.assertAlmostEqual(simulation.a3a, 0.32720000000000005)


        LHS= [[854672.82758621,     0. ,     0. ,     0. ,     0. ,     0. ],
              [    0. , 1068341.03448276,     0. ,     0. ,     0. ,     0. ],
              [    0. ,     0. , 2136682.06896552,     0. ,     0. ,     0. ],
              [    0. ,     0. ,     0. , 43144.4137931,     0. ,     0. ],
              [    0. ,     0. ,     0. ,     0. , 215722.06896552,     0. ],
              [    0. ,     0. ,     0. ,     0. ,     0. , 21344.86206897]]


        for i in range(6):
            for j in range(6):
                self.assertAlmostEqual(simulation.LHS[i][j], LHS[i][j])

    def test_check_variables(self):
        simulation = RigidBodySolver(self.model, self.default_parameters)
        # # Skip __init__
        # simulation = object.__new__(RigidBodySolver)
        # # Creating model
        # simulation.model = KM.Model()
        # simulation.main_model_part = simulation.model.CreateModelPart("Main")
        # simulation.rigid_body_model_part = simulation.main_model_part.CreateSubModelPart("RigidBody")
        # simulation.root_point_model_part = simulation.main_model_part.CreateSubModelPart("RootPoint")
        # # Adding variables
        # simulation.AddVariables()

        # Kinematic variables (work with both RigidBody and RootPoint model parts)
        self.assertTrue(simulation.main_model_part.HasNodalSolutionStepVariable(KM.DISPLACEMENT))
        self.assertTrue(simulation.main_model_part.HasNodalSolutionStepVariable(KM.ROTATION))
        self.assertTrue(simulation.main_model_part.HasNodalSolutionStepVariable(KM.VELOCITY))
        self.assertTrue(simulation.main_model_part.HasNodalSolutionStepVariable(KM.ANGULAR_VELOCITY))
        self.assertTrue(simulation.main_model_part.HasNodalSolutionStepVariable(KM.ACCELERATION))
        self.assertTrue(simulation.main_model_part.HasNodalSolutionStepVariable(KM.ANGULAR_ACCELERATION))
        
        # Specific variables for the model part RigidBody
        self.assertTrue(simulation.rigid_body_model_part.HasNodalSolutionStepVariable(KM.FORCE))
        self.assertTrue(simulation.rigid_body_model_part.HasNodalSolutionStepVariable(KM.MOMENT))
        self.assertTrue(simulation.rigid_body_model_part.HasNodalSolutionStepVariable(KMC.PRESCRIBED_FORCE))
        self.assertTrue(simulation.rigid_body_model_part.HasNodalSolutionStepVariable(KMC.PRESCRIBED_MOMENT))
        self.assertTrue(simulation.rigid_body_model_part.HasNodalSolutionStepVariable(KMC.EFFECTIVE_FORCE))
        self.assertTrue(simulation.rigid_body_model_part.HasNodalSolutionStepVariable(KMC.EFFECTIVE_MOMENT))
        self.assertTrue(simulation.rigid_body_model_part.HasNodalSolutionStepVariable(KM.BODY_FORCE))
        self.assertTrue(simulation.rigid_body_model_part.HasNodalSolutionStepVariable(KM.BODY_MOMENT))
        
        # Specific variables for the model part RootPoint
        self.assertTrue(simulation.root_point_model_part.HasNodalSolutionStepVariable(KM.REACTION))
        self.assertTrue(simulation.root_point_model_part.HasNodalSolutionStepVariable(KM.REACTION_MOMENT))
        self.assertTrue(simulation.root_point_model_part.HasNodalSolutionStepVariable(KMC.PRESCRIBED_DISPLACEMENT))
        self.assertTrue(simulation.root_point_model_part.HasNodalSolutionStepVariable(KMC.PRESCRIBED_ROTATION))

    def test_SetCompleteVector_rigid_body(self):
        simulation = RigidBodySolver(self.model, self.default_parameters)
        buffer = 0
        random_linear = np.random.rand(simulation.linear_size)
        random_angular = np.random.rand(simulation.angular_size)
        random_vector = np.array(list(random_linear) + list(random_angular))

        simulation._SetCompleteVector("rigid_body", KM.FORCE, KM.MOMENT, random_vector, buffer)
        simulation_linear = simulation.rigid_body_model_part.Nodes[1].GetSolutionStepValue(KM.FORCE, buffer)
        simulation_angular = simulation.rigid_body_model_part.Nodes[1].GetSolutionStepValue(KM.MOMENT, buffer)

        self.assertVectorAlmostEqual(random_linear, simulation_linear)
        self.assertVectorAlmostEqual(random_angular, simulation_angular)

    def test_SetCompleteVector_root_point(self):
        simulation = RigidBodySolver(self.model, self.default_parameters)
        buffer = 0
        random_linear = np.random.rand(simulation.linear_size)
        random_angular = np.random.rand(simulation.angular_size)
        random_vector = np.array(list(random_linear) + list(random_angular))

        simulation._SetCompleteVector("root_point", KM.DISPLACEMENT, KM.ROTATION, random_vector, buffer)
        simulation_linear = simulation.root_point_model_part.Nodes[2].GetSolutionStepValue(KM.DISPLACEMENT, buffer)
        simulation_angular = simulation.root_point_model_part.Nodes[2].GetSolutionStepValue(KM.ROTATION, buffer)

        self.assertVectorAlmostEqual(random_linear, simulation_linear)
        self.assertVectorAlmostEqual(random_angular, simulation_angular)

    def test_SetCompleteVector_false_model_part(self):
        simulation = RigidBodySolver(self.model, self.default_parameters)
        buffer = 0
        random_vector = np.random.rand(simulation.system_size)
        self.assertRaises(Exception, simulation._SetCompleteVector, "false_model_part", KM.DISPLACEMENT, KM.ROTATION, random_vector, buffer)

    def test_GetCompleteVector_rigid_body(self):
        simulation = RigidBodySolver(self.model, self.default_parameters)
        buffer = 0
        random_linear = np.random.rand(simulation.linear_size)
        random_angular = np.random.rand(simulation.angular_size)
        random_vector = np.array(list(random_linear) + list(random_angular))
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.FORCE, buffer, random_linear)
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.MOMENT, buffer, random_angular)

        simulation_vector = simulation._GetCompleteVector("rigid_body", KM.FORCE, KM.MOMENT, buffer)

        self.assertVectorAlmostEqual(random_vector, simulation_vector)
    
    def test_GetCompleteVector_root_point(self):
        simulation = RigidBodySolver(self.model, self.default_parameters)
        buffer = 0
        random_linear = np.random.rand(simulation.linear_size)
        random_angular = np.random.rand(simulation.angular_size)
        random_vector = np.array(list(random_linear) + list(random_angular))
        simulation.root_point_model_part.Nodes[2].SetSolutionStepValue(KM.DISPLACEMENT, buffer, random_linear)
        simulation.root_point_model_part.Nodes[2].SetSolutionStepValue(KM.ROTATION, buffer, random_angular)

        simulation_vector = simulation._GetCompleteVector("root_point", KM.DISPLACEMENT, KM.ROTATION, buffer)

        self.assertVectorAlmostEqual(random_vector, simulation_vector)

    def test_GetCompleteVector_false_model_part(self):
        simulation = RigidBodySolver(self.model, self.default_parameters)
        buffer = 0
        self.assertRaises(Exception, simulation._GetCompleteVector, "false_model_part", KM.DISPLACEMENT, KM.ROTATION, buffer)

    def test_ResetExternalVariables(self):
        simulation = RigidBodySolver(self.model, self.default_parameters)

        random_vector = np.random.rand(simulation.system_size)
        simulation._SetCompleteVector("rigid_body", KM.FORCE, KM.MOMENT, random_vector)
        simulation._SetCompleteVector("rigid_body", KMC.PRESCRIBED_FORCE, KMC.PRESCRIBED_MOMENT, random_vector)
        simulation._SetCompleteVector("root_point", KM.DISPLACEMENT, KM.ROTATION, random_vector)
        simulation._SetCompleteVector("root_point", KMC.PRESCRIBED_DISPLACEMENT, KMC.PRESCRIBED_ROTATION, random_vector)
        
        simulation._ResetExternalVariables()

        zero_vector = np.zeros(simulation.system_size)
        self.assertVectorAlmostEqual(zero_vector, simulation._GetCompleteVector("rigid_body", KM.FORCE, KM.MOMENT))
        self.assertVectorAlmostEqual(zero_vector, simulation._GetCompleteVector("rigid_body", KMC.PRESCRIBED_FORCE, KMC.PRESCRIBED_MOMENT))
        self.assertVectorAlmostEqual(zero_vector, simulation._GetCompleteVector("root_point", KM.DISPLACEMENT, KM.ROTATION))
        self.assertVectorAlmostEqual(zero_vector, simulation._GetCompleteVector("root_point", KMC.PRESCRIBED_DISPLACEMENT, KMC.PRESCRIBED_ROTATION))

    def test_CalculateEffectiveLoad(self):
        # Setting up the model
        Parameters = KM.Parameters('''{
            "solver_settings":{
                "active_dofs": [
                    {
                        "dof": "displacement_x",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 40.0,
                            "stiffness": 400.0,
                            "damping": 4.0
                        }
                    },
                    {
                        "dof": "displacement_y",
                        "constrained": true,
                        "system_parameters": {
                            "mass": 50.0,
                            "stiffness": 500.0,
                            "damping": 5.0
                        }
                    },
                    {
                        "dof": "displacement_z",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 100.0,
                            "stiffness": 1000.0,
                            "damping": 10.0
                        }
                    },
                    {
                        "dof": "rotation_x",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 2.0,
                            "stiffness": 200.0,
                            "damping": 2.0
                        }
                    },
                    {
                        "dof": "rotation_y",
                        "constrained": true,
                        "system_parameters": {
                            "mass": 10.0,
                            "stiffness": 1000.0,
                            "damping": 10.0
                        }
                    }
                ]
            },
            "processes": {
                "gravity": [
                    {
                        "python_module": "process_factory",
                        "kratos_module": "KratosMultiphysics",
                        "process_name": "ApplyConstantScalarValueProcess",
                        "Parameters": {
                            "model_part_name": "Main.RigidBody",
                            "variable_name": "BODY_FORCE_Y",
                            "is_fixed": true,
                            "value": -981
                        }
                    },
                    {
                        "python_module": "process_factory",
                        "kratos_module": "KratosMultiphysics",
                        "process_name": "ApplyConstantScalarValueProcess",
                        "Parameters": {
                            "model_part_name": "Main.RigidBody",
                            "variable_name": "BODY_MOMENT_Y",
                            "is_fixed": true,
                            "value": -98.1
                        }
                    }
                ],
                "initial_conditions_process_list": [
                    {
                        "python_module": "process_factory",
                        "kratos_module": "KratosMultiphysics",
                        "process_name": "ApplyConstantScalarValueProcess",
                        "Parameters": {
                            "model_part_name": "Main.RigidBody",
                            "variable_name": "DISPLACEMENT_X",
                            "value": 1
                        }
                    },
                    {
                        "python_module": "process_factory",
                        "kratos_module": "KratosMultiphysics",
                        "process_name": "ApplyConstantScalarValueProcess",
                        "Parameters": {
                            "model_part_name": "Main.RigidBody",
                            "variable_name": "ACCELERATION_X",
                            "value": -400
                        }
                    },
                    {
                        "python_module": "process_factory",
                        "kratos_module": "KratosMultiphysics",
                        "process_name": "ApplyConstantScalarValueProcess",
                        "Parameters": {
                            "model_part_name": "Main.RigidBody",
                            "variable_name": "ANGULAR_VELOCITY_Z",
                            "value": -10
                        }
                    }
                ],
                "boundary_conditions_process_list": [
                    {
                        "python_module" : "assign_vector_variable_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"  : "AssignVectorVariableProcess",
                        "Parameters": {
                            "model_part_name" : "Main.RigidBody",
                            "variable_name"   : "PRESCRIBED_FORCE",
                            "interval"        : [0, "End"],
                            "constrained"     : [false,false,false],
                            "value"           : ["5*sin((5+2*t)*t)", "3*sin((5+2*t)*t)", "2*sin((5+2*t)*t)"]
                        }
                    },
                    {
                        "python_module" : "assign_vector_variable_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"  : "AssignVectorVariableProcess",
                        "Parameters": {
                            "model_part_name" : "Main.RigidBody",
                            "variable_name"   : "FORCE",
                            "interval"        : [0, "End"],
                            "constrained"     : [false,false,false],
                            "value"           : ["5", "3", "2"]
                        }
                    },
                    {
                        "python_module" : "assign_vector_variable_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"  : "AssignVectorVariableProcess",
                        "Parameters": {
                            "model_part_name" : "Main.RigidBody",
                            "variable_name"   : "MOMENT",
                            "interval"        : [0, "End"],
                            "constrained"     : [false,false,false],
                            "value"           : ["4", "3", "5"]
                        }
                    },
                    {
                        "python_module": "assign_vector_variable_process",
                        "kratos_module": "KratosMultiphysics",
                        "process_name": "AssignVectorVariableProcess",
                        "Parameters": {
                            "model_part_name": "Main.RootPoint",
                            "variable_name": "PRESCRIBED_DISPLACEMENT",
                            "interval": [
                                0,
                                "End"
                            ],
                            "constrained": [
                                false,
                                false,
                                false
                            ],
                            "value": [
                                "0.05*sin(20*t)",
                                "0.05*sin(10*t)",
                                "0.05*sin(5*t)"
                            ]
                        }
                    },
                    {
                        "python_module": "assign_vector_variable_process",
                        "kratos_module": "KratosMultiphysics",
                        "process_name": "AssignVectorVariableProcess",
                        "Parameters": {
                            "model_part_name": "Main.RootPoint",
                            "variable_name": "PRESCRIBED_ROTATION",
                            "interval": [
                                0,
                                "End"
                            ],
                            "constrained": [
                                false,
                                false,
                                false
                            ],
                            "value": [
                                "0.05*sin(6*t)",
                                "0.05*sin(8*t)",
                                "0.05*sin(13*t)"
                            ]
                        }
                    },
                    {
                        "python_module": "assign_vector_variable_process",
                        "kratos_module": "KratosMultiphysics",
                        "process_name": "AssignVectorVariableProcess",
                        "Parameters": {
                            "model_part_name": "Main.RigidBody",
                            "variable_name": "PRESCRIBED_MOMENT",
                            "interval": [
                                0,
                                "End"
                            ],
                            "constrained": [
                                false,
                                false,
                                false
                            ],
                            "value": [
                                "100*sin(4*t)",
                                "100*sin(3*t)",
                                "100*sin(10*t)"
                            ]
                        }
                    }
                ]
            }
        }''')
        Parameters.RecursivelyAddMissingParameters(self.default_parameters)
        simulation = RigidBodySolver(self.model, Parameters)
        simulation.Initialize()

        while simulation.main_model_part.ProcessInfo[KM.TIME] < 0.1:
            simulation.AdvanceInTime(simulation.main_model_part.ProcessInfo[KM.TIME])
            simulation.InitializeSolutionStep()
            simulation.Predict()
            simulation.SolveSolutionStep()
            simulation.FinalizeSolutionStep()
            simulation.OutputSolutionStep()

        simulation.AdvanceInTime(simulation.main_model_part.ProcessInfo[KM.TIME])
        simulation.InitializeSolutionStep()
        simulation.Predict()
        simulation.SolveSolutionStep()

        reference = [14.47730397, -954.52993411, 30.97172112, 56.62280302, -20.61982872, 98.25390568]
        simulation_effective_load = simulation._GetCompleteVector("rigid_body", KMC.EFFECTIVE_FORCE, KMC.EFFECTIVE_MOMENT)

        self.assertVectorAlmostEqual(reference, simulation_effective_load)

    def test_GetKinematics(self):
        simulation = RigidBodySolver(self.model, self.default_parameters)
        buffer = 0
        random_linear = np.random.rand(simulation.linear_size)
        random_angular = np.random.rand(simulation.angular_size)
        random_vector_x = np.array(list(random_linear) + list(random_angular))
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.DISPLACEMENT, buffer, random_linear)
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.ROTATION, buffer, random_angular)

        random_linear = np.random.rand(simulation.linear_size)
        random_angular = np.random.rand(simulation.angular_size)
        random_vector_v = np.array(list(random_linear) + list(random_angular))
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.VELOCITY, buffer, random_linear)
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.ANGULAR_VELOCITY, buffer, random_angular)

        random_linear = np.random.rand(simulation.linear_size)
        random_angular = np.random.rand(simulation.angular_size)
        random_vector_a = np.array(list(random_linear) + list(random_angular))
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.ACCELERATION, buffer, random_linear)
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.ANGULAR_ACCELERATION, buffer, random_angular)

        x, v, a = simulation._GetKinematics("rigid_body", 0)

        self.assertVectorAlmostEqual(random_vector_x, x)
        self.assertVectorAlmostEqual(random_vector_v, v)
        self.assertVectorAlmostEqual(random_vector_a, a)

    def test_UpdateKinematics(self):
        simulation = RigidBodySolver(self.model, self.default_parameters)
        # Sets previous kinematics
        buffer = 1
        prev_linear = np.array([0.3496 , 0.575, 0.7984])
        prev_angular = np.array([0.4725 , 0.6698, 0.3468])
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.DISPLACEMENT, buffer, prev_linear)
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.ROTATION, buffer, prev_angular)

        prev_linear = np.array([1.7849 , 2.2954, 3.3984])
        prev_angular = np.array([0.3989 , 0.8954, 1.1389])
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.VELOCITY, buffer, prev_linear)
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.ANGULAR_VELOCITY, buffer, prev_angular)

        prev_linear = np.array([3.1488 , 0.5647, 0.1417])
        prev_angular = np.array([1.8998 , 2.1987, 0.9874])
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.ACCELERATION, buffer, prev_linear)
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.ANGULAR_ACCELERATION, buffer, prev_angular)

        current_x = np.array([0.4861 , 0.6813, 0.8942, 0.7549 , 1.1698, 0.7918])

        simulation._UpdateKinematics("rigid_body", current_x)

        ref_DISPLACEMENT = np.array([0.4861, 0.6813, 0.8942])
        ref_ROTATION = np.array([0.7549 , 1.1698, 0.7918])
        ref_VELOCITY = np.array([25.5151, 18.9646, 15.7616])
        ref_ANGULAR_VELOCITY = np.array([56.0811, 99.1046, 87.8611])
        ref_ACCELERATION = np.array([4742.8912, 3333.2753, 2472.4983])
        ref_ANGULAR_ACCELERATION = np.array([11134.5402, 19639.6413, 17343.4526])

        buffer = 0
        simulation_DISPLACEMENT = simulation.rigid_body_model_part.Nodes[1].GetSolutionStepValue(KM.DISPLACEMENT, buffer)
        simulation_ROTATION = simulation.rigid_body_model_part.Nodes[1].GetSolutionStepValue(KM.ROTATION, buffer)
        simulation_VELOCITY = simulation.rigid_body_model_part.Nodes[1].GetSolutionStepValue(KM.VELOCITY, buffer)
        simulation_ANGULAR_VELOCITY = simulation.rigid_body_model_part.Nodes[1].GetSolutionStepValue(KM.ANGULAR_VELOCITY, buffer)
        simulation_ACCELERATION = simulation.rigid_body_model_part.Nodes[1].GetSolutionStepValue(KM.ACCELERATION, buffer)
        simulation_ANGULAR_ACCELERATION = simulation.rigid_body_model_part.Nodes[1].GetSolutionStepValue(KM.ANGULAR_ACCELERATION, buffer)

        self.assertVectorAlmostEqual(ref_DISPLACEMENT, simulation_DISPLACEMENT)
        self.assertVectorAlmostEqual(ref_ROTATION, simulation_ROTATION)
        self.assertVectorAlmostEqual(ref_VELOCITY, simulation_VELOCITY)
        self.assertVectorAlmostEqual(ref_ANGULAR_VELOCITY, simulation_ANGULAR_VELOCITY)
        self.assertVectorAlmostEqual(ref_ACCELERATION, simulation_ACCELERATION)
        self.assertVectorAlmostEqual(ref_ANGULAR_ACCELERATION, simulation_ANGULAR_ACCELERATION)

    def test_CalculateReaction(self):
        Parameters = KM.Parameters('''{
            "solver_settings": {
                "time_integration_parameters": {
                        "rho_inf": 0.16,
                        "time_step": 0.01
                },
                "active_dofs": [
                    {
                        "dof": "displacement_x",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 40.0,
                            "stiffness": 400.0,
                            "damping": 4.0
                        }
                    },
                    {
                        "dof": "displacement_z",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 100.0,
                            "stiffness": 1000.0,
                            "damping": 10.0
                        }
                    },
                    {
                        "dof": "rotation_x",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 2.0,
                            "stiffness": 200.0,
                            "damping": 2.0
                        }
                    },
                    {
                        "dof": "rotation_y",
                        "constrained": true,
                        "system_parameters": {
                            "mass": 10.0,
                            "stiffness": 1000.0,
                            "damping": 10.0
                        }
                    }
                ]
            }
        }''')
        
        Parameters.RecursivelyAddMissingParameters(self.default_parameters)
        simulation = RigidBodySolver(self.model, Parameters)

        buffer = 0
        linear = np.array([3.1488 , 0.5647, 0.1417])
        angular = np.array([1.8998 , 2.1987, 0.9874])
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.ACCELERATION, buffer, linear)
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.ANGULAR_ACCELERATION, buffer, angular)

        buffer = 1
        linear = np.array([2.1488 , 1.5647, 0.5417])
        angular = np.array([2.8998 , 3.1987, 1.5874])
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.ACCELERATION, buffer, linear)
        simulation.rigid_body_model_part.Nodes[1].SetSolutionStepValue(KM.ANGULAR_ACCELERATION, buffer, angular)

        simulation.total_load = 1.4648

        ref_buffer_0 = np.array([-124.4872, 0., -12.7052, -2.3348, -20.5222, 0.])
        ref_buffer_1 = np.array([-84.4872, 0., -52.7052, -4.3348, -30.5222, 0.])

        self.assertVectorAlmostEqual(ref_buffer_0, simulation._CalculateReaction(0))
        self.assertVectorAlmostEqual(ref_buffer_1, simulation._CalculateReaction(1))

    def test_CalculateEquivalentForceFromRootPointDisplacement(self):
        Parameters = KM.Parameters('''{
            "solver_settings":{
                "active_dofs": [
                    {
                        "dof": "displacement_x",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 40.0,
                            "stiffness": 400.0,
                            "damping": 4.0
                        }
                    },
                    {
                        "dof": "displacement_y",
                        "constrained": true,
                        "system_parameters": {
                            "mass": 50.0,
                            "stiffness": 500.0,
                            "damping": 5.0
                        }
                    },
                    {
                        "dof": "displacement_z",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 100.0,
                            "stiffness": 1000.0,
                            "damping": 10.0
                        }
                    },
                    {
                        "dof": "rotation_x",
                        "constrained": false,
                        "system_parameters": {
                            "mass": 2.0,
                            "stiffness": 200.0,
                            "damping": 2.0
                        }
                    },
                    {
                        "dof": "rotation_y",
                        "constrained": true,
                        "system_parameters": {
                            "mass": 10.0,
                            "stiffness": 1000.0,
                            "damping": 10.0
                        }
                    }
                ]
            }
        }''')

        Parameters.RecursivelyAddMissingParameters(self.default_parameters)
        simulation = RigidBodySolver(self.model, Parameters)

        buffer = 0
        linear = np.array([0.3496 , 0.575, 0.7984])
        angular = np.array([0.4725 , 0.6698, 0.3468])
        simulation.root_point_model_part.Nodes[2].SetSolutionStepValue(KM.DISPLACEMENT, buffer, linear)
        simulation.root_point_model_part.Nodes[2].SetSolutionStepValue(KM.ROTATION, buffer, angular)

        linear = np.array([1.7849 , 2.2954, 3.3984])
        angular = np.array([0.3989 , 0.8954, 1.1389])
        simulation.root_point_model_part.Nodes[2].SetSolutionStepValue(KM.VELOCITY, buffer, linear)
        simulation.root_point_model_part.Nodes[2].SetSolutionStepValue(KM.ANGULAR_VELOCITY, buffer, angular)

        ref_equivalent_force = np.array([1.469796e+02, 2.989770e+02, 8.323840e+02, 9.529780e+01, 6.787540e+02, 3.468000e-01])
        
        self.assertVectorAlmostEqual(ref_equivalent_force, simulation._CalculateEquivalentForceFromRootPointDisplacement())

        


# class TestSdofSolver(KratosUnittest.TestCase):

#     def setUp(self):
#         self.system_settings = {
#         "system_parameters":{
#             "mass"      : 10.0,
#             "stiffness" : 1579.14,
#             "damping"   : 10.0
#             },
#         "time_integration_parameters":{
#                     "alpha_m"   : -0.3,
#                     "start_time": 0.0,
#                     "time_step" : 0.01
#                 },
#         "output_parameters":{
#             "write_output_file": True,
#             "file_name" : "result.dat"
#             }
#         }
#         #result.dat
#         self.end_time = 1.0
#         self.time = 0.0

#     def tearDown(self):
#         kratos_utilities.DeleteFileIfExisting("result.dat")

#     def __CompareResults(self, reference, result):
#         ref = np.loadtxt(reference, skiprows=1)
#         res = np.loadtxt(result, skiprows=1)
#         self.assertEqual(len(ref), len(res))
#         for line_ref, line_res in zip(ref, res):
#             self.assertEqual(len(line_ref), len(line_res))
#             for entry_ref, entry_res in zip(line_ref, line_res):
#                 self.assertAlmostEqual(entry_ref, entry_res)

#     def __ExecuteTest(self, settings, ref_file_name):
#         settings.update(self.system_settings)
#         system = SDoFSolver(settings)
#         system.Initialize()

#         while(self.time <= self.end_time):
#             self.time =  system.AdvanceInTime(self.time)
#             system.SolveSolutionStep()
#             system.OutputSolutionStep()

#         self.__CompareResults(os.path.join("reference_files", ref_file_name), "result.dat")


#     def test_initial_displacement(self):
#         settings = {
#         "initial_values":{
#             "displacement"  : 1.0,
#             "velocity"      : 0.0,
#             "acceleration"  : 0.0
#             }
#         }
#         self.__ExecuteTest(settings, "ref_sdof_initial_displacement.dat")

#     def test_initial_velocity(self):
#         settings = {
#         "initial_values":{
#             "displacement"  : 0.0,
#             "velocity"      : 1.0,
#             "acceleration"  : 0.0
#             }
#         }
#         self.__ExecuteTest(settings, "ref_sdof_initial_velocity.dat")

#     def test_impulse(self):
#         settings = {
#                 "boundary_conditions":{
#                     "load_impulse" : 1000.0
#                 }
#         }
#         self.__ExecuteTest(settings, "ref_sdof_impulse.dat")

#     def test_force_excitation(self):
#         settings = {
#                 "boundary_conditions":{
#                     "load_impulse" : 0.0,
#                     "omega_force"        : 12.57,
#                     "omega_root_point_displacement"        : 0.0,
#                     "excitation_function_force": "A * sin(omega * t)",
#                     "excitation_function_root_point_displacement": "A * sin(omega * t)",
#                     "amplitude_root_point_displacement": 0.0,
#                     "amplitude_force": 1000.0
#                 }
#         }
#         self.__ExecuteTest(settings, "ref_sdof_force_excitation.dat")

#     def test_root_point_excitation(self):
#         settings = {
#                 "boundary_conditions":{
#                     "load_impulse" : 0.0,
#                     "omega_force"        : 0.0,
#                     "omega_root_point_displacement"        : 12.57,
#                     "excitation_function_force": "A * sin(omega * t)",
#                     "excitation_function_root_point_displacement": "A * sin(omega * t)",
#                     "amplitude_root_point_displacement": 1.0,
#                     "amplitude_force": 0.0
#                 }
#         }
#         self.__ExecuteTest(settings, "ref_sdof_root_point_displacement.dat")

#     def test_root_point_excitation_force_excitation(self):
#         settings = {
#                 "boundary_conditions":{
#                     "load_impulse" : 0.0,
#                     "omega_force"        : 3.0,
#                     "omega_root_point_displacement"        : 12.57,
#                     "excitation_function_force": "A * sin(omega * t)",
#                     "excitation_function_root_point_displacement": "A * sin(omega * t)",
#                     "amplitude_root_point_displacement": 1.0,
#                     "amplitude_force": 1000.0
#                 }
#         }
#         self.__ExecuteTest(settings, "ref_sdof_root_point_displacement_external_force.dat")

#     def test_root_point_excitation_force_excitation_impulse(self):
#         settings = {
#                 "boundary_conditions":{
#                     "load_impulse" : 10.0,
#                     "omega_force"        : 3.0,
#                     "omega_root_point_displacement"        : 12.57,
#                     "excitation_function_force": "A * sin(omega * t)",
#                     "excitation_function_root_point_displacement": "A * sin(omega * t)",
#                     "amplitude_root_point_displacement": 1.0,
#                     "amplitude_force": 1000.0
#                 }
#         }
#         self.__ExecuteTest(settings, "ref_sdof_root_point_displacement_impulse.dat")

#     def test_self_weight_calculation(self):
#         system = SDoFSolver(self.system_settings)
#         system.Initialize()
#         self_weight = system.GetSolutionStepValue("VOLUME_ACCELERATION")
#         self.assertAlmostEqual(98.10000000000001, self_weight)

#     def test_wrapper_variables_check(self):
#         # resusing the sdof fsi parameters
#         wrapper_settings = KM.Parameters("""{
#             "type" : "solver_wrappers.sdof.sdof_solver_wrapper",
#             "solver_wrapper_settings" : {
#                 "input_file"  : "fsi_sdof/ProjectParametersSdof"
#             },
#             "data" : {
#                 "disp" : {
#                     "model_part_name" : "Sdof",
#                     "variable_name" : "DISPLACEMENT",
#                     "location"      : "model_part",
#                     "dimension"     : 1
#                 }
#             }
#         }""")
#         sdof_solver_wrapper = CreateSDofSolverWrapper(wrapper_settings, None, "custom_sdof_solver_wrapper")

#         sdof_solver_wrapper.Initialize()
#         with self.assertRaisesRegex(Exception, 'Variable "DISPLACEMENT" of interface data "disp" of solver "custom_sdof_solver_wrapper" cannot be used for the SDof Solver!'):
#             sdof_solver_wrapper.Check()

if __name__ == '__main__':
    KratosUnittest.main()
