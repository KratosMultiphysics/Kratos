from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

from KratosMultiphysics import kratos_utilities
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.sdof.sdof_static_solver import SDoFStaticSolver

import os
import numpy as np

class TestSdofStaticSolver(KratosUnittest.TestCase):

    def setUp(self):
        self.system_settings = {
        "system_parameters":{
            "stiffness" : 50000.0,
            },
        "output_parameters":{
            "write_output_file": True,
            "file_name" : "result.dat"
            }
        }
        #result.dat
        self.end_time = 1.0
        self.time = 0.0

    @classmethod
    def tearDownClass(self):
        kratos_utilities.DeleteFileIfExisting("result.dat")
        kratos_utilities.DeleteFileIfExisting('fsi_sdof_static/results_final_sdof.dat')

    def __CompareResults(self, reference, result):
        ref = np.loadtxt(reference, skiprows=1)
        res = np.loadtxt(result, skiprows=1)
        self.assertEqual(ref.all(), res.all())

    def __ExecuteTest(self, settings, ref_file_name):
        settings.update(self.system_settings)
        system = SDoFStaticSolver(settings)
        system.Initialize()

        system.SolveSolutionStep()
        system.OutputSolutionStep()
        self.__CompareResults(os.path.join("reference_files", ref_file_name), "result.dat")


    def test_initial_displacement(self):
        settings = {
        "initial_values":{
            "displacement"  : 1.0,
            }
        }
        self.__ExecuteTest(settings, "ref_sdof_static_initial_displacement.dat")

    def test_final_displacement(self):
        import json
        parameter_file_name = "fsi_sdof_static/ProjectParametersSDoF.json"
        with open(parameter_file_name, 'r') as parameter_file:
            settings = json.load(parameter_file)

        settings["output_parameters"]["write_output_file"] = True

        system = SDoFStaticSolver(settings)
        system.Initialize()

        system.SolveSolutionStep()
        system.OutputSolutionStep()

        results_obtained = np.loadtxt('fsi_sdof_static/results_final_sdof.dat', skiprows=1)
        results_reference = np.loadtxt('reference_files/ref_sdof_static_final_displacement.dat', skiprows=1)
        self.assertEqual(results_reference.all(), results_obtained.all())

if __name__ == '__main__':
    KratosUnittest.main()
