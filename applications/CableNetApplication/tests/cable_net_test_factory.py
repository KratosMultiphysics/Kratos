from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils


import KratosMultiphysics.CableNetApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.DEMApplication

from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import cable_net_test_case
import os

try:
    import numpy
    numpy_available = True
except ImportError:
    numpy_available = False


have_dem_fem_dependencies = kratos_utils.CheckIfApplicationsAvailable("DEMApplication", "StructuralMechanicsApplication", "MappingApplication", "ExternalSolversApplication")


def GetFilePath(fileName):
    return os.path.join(os.path.dirname(__file__), fileName)


class TestCableNetCoSimulationCases(cable_net_test_case.CableNetTestCase):

    def _createTest(self, problem_dir_name, parameter_file_name):
        self.problem_dir_name = problem_dir_name

        full_parameter_file_name = os.path.join(problem_dir_name, parameter_file_name + '_parameters.json')

        with open(full_parameter_file_name, 'r') as parameter_file:
            self.cable_net_parameters = KM.Parameters(parameter_file.read())

        # To avoid many prints
        echo_level = self.cable_net_parameters["problem_data"]["echo_level"].GetInt()
        if (echo_level == 0):
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
        else:
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.INFO)

        self.test = CoSimulationAnalysis(self.cable_net_parameters)


    def test_DEMFEMCableNet_SlidingEdges(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_dem_fem_dependencies:
            self.skipTest("DEM FEM dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("sliding_edges_with_friction_dem_fem","cosim_dem_fem_cable_net")
            self._runTest()

    def test_DEMFEMCableNet_RingElements(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_dem_fem_dependencies:
            self.skipTest("DEM FEM dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("ring_set_up_dem_fem","cosim_dem_fem_cable_net")
            self._runTest()

    @classmethod
    def tearDownClass(cls):
        super(TestCableNetCoSimulationCases,cls).tearDownClass()

        # delete superfluous dem files
        dir_name_list = ["sliding_edges_with_friction_dem_fem","ring_set_up_dem_fem"]
        project_name_list = ["slidingDEM", "smallNet"]

        for dir_name, project_name in zip(dir_name_list, project_name_list):
            kratos_utils.DeleteFileIfExisting(GetFilePath(dir_name+"/"+project_name+".post.lst"))
            kratos_utils.DeleteDirectoryIfExisting(GetFilePath(dir_name+"/"+project_name+"_Graphs"))
            kratos_utils.DeleteDirectoryIfExisting(GetFilePath(dir_name+"/"+project_name+"_MPI_results"))
            kratos_utils.DeleteDirectoryIfExisting(GetFilePath(dir_name+"/"+project_name+"_Post_Files"))
            kratos_utils.DeleteDirectoryIfExisting(GetFilePath(dir_name+"/"+project_name+"_Results_and_Data"))

class TestCableNetFEMCases(cable_net_test_case.CableNetTestCase):

        def _createTest(self, problem_dir_name):
            self.problem_dir_name = problem_dir_name

            full_parameter_file_name = os.path.join(problem_dir_name,'ProjectParameters.json')

            with open(full_parameter_file_name, 'r') as parameter_file:
                self.cable_net_parameters = KM.Parameters(parameter_file.read())

            # To avoid many prints
            echo_level = self.cable_net_parameters["problem_data"]["echo_level"].GetInt()
            if (echo_level == 0):
                KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
            else:
                KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.INFO)

            model = KM.Model()
            self.test = StructuralMechanicsAnalysis(model, self.cable_net_parameters)


        def test_SlidingNodes(self):
            if not numpy_available:
                self.skipTest("Numpy not available")

            with KratosUnittest.WorkFolderScope(".", __file__):
                self._createTest("sliding_element")
                self._runTest()

        def test_RingElement(self):
            if not numpy_available:
                self.skipTest("Numpy not available")

            with KratosUnittest.WorkFolderScope(".", __file__):
                self._createTest("ring_element")
                self._runTest()

if __name__ == '__main__':
    KratosUnittest.main()