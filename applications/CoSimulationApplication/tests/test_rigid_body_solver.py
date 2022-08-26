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


    def test_InitializeGeneralizedAlphaParameters(self):

        simulation = RigidBodySolver(self.model, self.default_parameters)

        self.assertEqual(simulation.alpha_f, 0.5)
        self.assertEqual(simulation.alpha_m, 0.5)
        self.assertEqual(simulation.beta, 0.25)
        self.assertEqual(simulation.gamma, 0.5)
        self.assertEqual(simulation.a1h, 20000.0)
        self.assertEqual(simulation.a2h, 100.0)
        self.assertEqual(simulation.a3h, 0.5)
        self.assertEqual(simulation.a1m, 20000.0)
        self.assertEqual(simulation.a2m, 200.0)
        self.assertEqual(simulation.a3m, 0.0)
        self.assertEqual(simulation.a1b, 100.0)
        self.assertEqual(simulation.a2b, 0.0)
        self.assertEqual(simulation.a3b, 0.0)
        self.assertEqual(simulation.a1k, -0.5)
        self.assertEqual(simulation.a1v, 200.0)
        self.assertEqual(simulation.a2v, -1.0)
        self.assertEqual(simulation.a3v, 0.0)
        self.assertEqual(simulation.a1a, 40000.0)
        self.assertEqual(simulation.a2a, -400.0)
        self.assertEqual(simulation.a3a, -1.0)

        LHS= [[20000.5,     0. ,     0. ,     0. ,     0. ,     0. ],
              [    0. , 20000.5,     0. ,     0. ,     0. ,     0. ],
              [    0. ,     0. , 20000.5,     0. ,     0. ,     0. ],
              [    0. ,     0. ,     0. , 20000.5,     0. ,     0. ],
              [    0. ,     0. ,     0. ,     0. , 20000.5,     0. ],
              [    0. ,     0. ,     0. ,     0. ,     0. , 20000.5]]


        for i in range(6):
            for j in range(6):
                self.assertEqual(simulation.LHS[i][j], LHS[i][j])

    def test_check_variables(self):
        # HasNodalSolutionStepVariable
        simulation = RigidBodySolver(self.model, self.default_parameters)


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
        self.assertTrue(simulation.rigid_body_model_part.HasNodalSolutionStepVariable(KM.REACTION))
        self.assertTrue(simulation.rigid_body_model_part.HasNodalSolutionStepVariable(KM.REACTION_MOMENT))
        self.assertTrue(simulation.rigid_body_model_part.HasNodalSolutionStepVariable(KMC.PRESCRIBED_DISPLACEMENT))
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
        pass

    def test_GetKinematics(self):
        pass

    def test_UpdateKinematics(self):
        pass

    def test_CalculateReaction(self):
        pass

    def test_CalculateEquivalentForceFromRootPointDisplacement(self):
        pass

        


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
