from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils


import KratosMultiphysics.CableNetApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.DEMApplication

from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

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

    def _runTest(self):
        super(TestCableNetCoSimulationCases,self)._runTest()
        CoSimulationAnalysis(self.cable_net_parameters).Run()

    def test_DEMFEMCableNet(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_dem_fem_dependencies:
            self.skipTest("DEM FEM dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("sliding_edges_with_friction_dem_fem","cosim_dem_fem_cable_net")
            self._runTest()

    @classmethod
    def tearDownClass(cls):
        super(TestCableNetCoSimulationCases,cls).tearDownClass()

        # delete superfluous dem files
        kratos_utils.DeleteFileIfExisting(GetFilePath("sliding_edges_with_friction_dem_fem/slidingDEM.post.lst"))
        kratos_utils.DeleteDirectoryIfExisting(GetFilePath("sliding_edges_with_friction_dem_fem/slidingDEM_Graphs"))
        kratos_utils.DeleteDirectoryIfExisting(GetFilePath("sliding_edges_with_friction_dem_fem/slidingDEM_MPI_results"))
        kratos_utils.DeleteDirectoryIfExisting(GetFilePath("sliding_edges_with_friction_dem_fem/slidingDEM_Post_Files"))
        kratos_utils.DeleteDirectoryIfExisting(GetFilePath("sliding_edges_with_friction_dem_fem/slidingDEM_Results_and_Data"))

if __name__ == '__main__':
    KratosUnittest.main()