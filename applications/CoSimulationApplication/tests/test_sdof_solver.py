from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.sdof.sdof_solver import SDoFSolver

import os
import numpy as np

class TestSdofSolver(KratosUnittest.TestCase):

    def setUp(self):
        self.system_settings = {
        "system_parameters":{
            "mass"      : 10.0,
            "stiffness" : 1579.14,
            "damping"   : 10.0
            },
        "time_integration_parameters":{
                    "alpha_m"   : -0.3,
                    "start_time": 0.0,
                    "time_step" : 0.01
                },
        "output_parameters":{
            "write_output_file": True,
            "file_name" : "./result.dat"
            }
        }
        #result.dat
        self.end_time = 1.0
        self.time = 0.0

    def tearDown(self):
        os.remove("result.dat")

    def __CompareResults(self, reference, result):
        ref = np.loadtxt(reference, skiprows=1)
        res = np.loadtxt(result, skiprows=1)
        self.assertEqual(len(ref), len(res))
        for line_ref, line_res in zip(ref, res):
            self.assertEqual(len(line_ref), len(line_res))
            for entry_ref, entry_res in zip(line_ref, line_res):
                self.assertAlmostEqual(entry_ref, entry_res)

    def __ExecuteTest(self, settings, ref_file_name):
        settings.update(self.system_settings)
        system = SDoFSolver(settings)
        system.Initialize()

        while(self.time <= self.end_time):
            self.time =  system.AdvanceInTime(self.time)
            system.SolveSolutionStep()
            system.OutputSolutionStep()

        self.__CompareResults(os.path.join("reference_files", ref_file_name), "result.dat")


    def test_initial_displacement(self):
        settings = {
        "initial_values":{
            "displacement"  : 1.0,
            "velocity"      : 0.0,
            "acceleration"  : 0.0
            }
        }
        self.__ExecuteTest(settings, "ref_sdof_initial_displacement.dat")

    def test_initial_velocity(self):
        settings = {
        "initial_values":{
            "displacement"  : 0.0,
            "velocity"      : 1.0,
            "acceleration"  : 0.0
            }
        }
        self.__ExecuteTest(settings, "ref_sdof_initial_velocity.dat")

    def test_impulse(self):
        settings = {
                "boundary_conditions":{
                    "load_impulse" : 1000.0
                }
        }
        self.__ExecuteTest(settings, "ref_sdof_impulse.dat")

    def test_force_excitation(self):
        settings = {
                "boundary_conditions":{
                    "load_impulse" : 0.0,
                    "omega_force"        : 12.57,
                    "omega_root_point_displacement"        : 0.0,
                    "excitation_function_force": "A * sin(omega * t)",
                    "excitation_function_root_point_displacement": "A * sin(omega * t)",
                    "amplitude_root_point_displacement": 0.0,
                    "amplitude_force": 1000.0
                }
        }
        self.__ExecuteTest(settings, "ref_sdof_force_excitation.dat")

    def test_root_point_excitation(self):
        settings = {
                "boundary_conditions":{
                    "load_impulse" : 0.0,
                    "omega_force"        : 0.0,
                    "omega_root_point_displacement"        : 12.57,
                    "excitation_function_force": "A * sin(omega * t)",
                    "excitation_function_root_point_displacement": "A * sin(omega * t)",
                    "amplitude_root_point_displacement": 1.0,
                    "amplitude_force": 0.0
                }
        }
        self.__ExecuteTest(settings, "ref_sdof_root_point_displacement.dat")

    def test_root_point_excitation_force_excitation(self):
        settings = {
                "boundary_conditions":{
                    "load_impulse" : 0.0,
                    "omega_force"        : 3.0,
                    "omega_root_point_displacement"        : 12.57,
                    "excitation_function_force": "A * sin(omega * t)",
                    "excitation_function_root_point_displacement": "A * sin(omega * t)",
                    "amplitude_root_point_displacement": 1.0,
                    "amplitude_force": 1000.0
                }
        }
        self.__ExecuteTest(settings, "ref_sdof_root_point_displacement_external_force.dat")

    def test_root_point_excitation_force_excitation_impulse(self):
        settings = {
                "boundary_conditions":{
                    "load_impulse" : 10.0,
                    "omega_force"        : 3.0,
                    "omega_root_point_displacement"        : 12.57,
                    "excitation_function_force": "A * sin(omega * t)",
                    "excitation_function_root_point_displacement": "A * sin(omega * t)",
                    "amplitude_root_point_displacement": 1.0,
                    "amplitude_force": 1000.0
                }
        }
        self.__ExecuteTest(settings, "ref_sdof_root_point_displacement_impulse.dat")

if __name__ == '__main__':
    KratosUnittest.main()
