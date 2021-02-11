import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

import co_simulation_test_case
import os

try:
    import numpy
    numpy_available = True
except ImportError:
    numpy_available = False

have_fsi_dependencies = kratos_utils.CheckIfApplicationsAvailable("FluidDynamicsApplication", "StructuralMechanicsApplication", "MappingApplication", "MeshMovingApplication", "LinearSolversApplication")
have_potential_fsi_dependencies = kratos_utils.CheckIfApplicationsAvailable("CompressiblePotentialFlowApplication", "StructuralMechanicsApplication", "MappingApplication", "MeshMovingApplication", "LinearSolversApplication")
have_mpm_fem_dependencies = kratos_utils.CheckIfApplicationsAvailable("ParticleMechanicsApplication", "StructuralMechanicsApplication", "MappingApplication", "LinearSolversApplication")
have_dem_fem_dependencies = kratos_utils.CheckIfApplicationsAvailable("DEMApplication", "StructuralMechanicsApplication", "MappingApplication", "LinearSolversApplication")
have_fem_fem_dependencies = kratos_utils.CheckIfApplicationsAvailable("StructuralMechanicsApplication", "MappingApplication")

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestTinyFetiCoSimulationCases(co_simulation_test_case.CoSimulationTestCase):
    '''This class contains "tiny" FETI CoSimulation-Cases, small enough to run in the CI
    '''
    def test_FEM_FEM_small_2d_plate_feti_explict_explicit(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_small_2d_plate_feti_explict_explicit"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/small_2d_plate_feti/explicit_explicit", "cosim_fem_fem_small_2d_plate_feti_explicit_explicit")
            self._runTest()

    def test_FEM_FEM_small_2d_plate_feti_implict_explicit(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_small_2d_plate_feti_implict_explicit"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/small_2d_plate_feti/implicit_explicit", "cosim_fem_fem_small_2d_plate_feti_implicit_explicit")
            self._runTest()

    def test_FEM_FEM_small_2d_plate_feti_implict_implicit(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_small_2d_plate_feti_implict_implicit"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/small_2d_plate_feti/implicit_implicit", "cosim_fem_fem_small_2d_plate_feti_implicit_implicit")
            self._runTest()

    def test_FEM_FEM_small_2d_plate_feti_implict_explicit_mixed(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_small_2d_plate_feti_implict_explicit_mixed"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/small_2d_plate_feti/implicit_explicit_mixed", "cosim_fem_fem_small_2d_plate_feti_implicit_explicit_mixed")
            self._runTest()


class TestSmallCoSimulationCases(co_simulation_test_case.CoSimulationTestCase):
    '''This class contains "small" CoSimulation-Cases, small enough to run in the nightly suite
    '''

    def test_MPM_FEM_beam_penalty(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_mpm_fem_dependencies:
            self.skipTest("MPM-FEM dependencies are not available!")

        self.name = "penalty_beam"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("mpm_fem_beam", "cosim_mpm_fem_beam")
            self._runTest()

    def test_FEM_FEM_small_2d_plate_dual_mortar(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_small_2d_plate_dual_mortar"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/small_2d_plate", "cosim_fem_fem_small_2d_plate_dual_mortar")
            self._runTest()

    def test_FEM_FEM_small_2d_plate_full_mortar(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_small_2d_plate_full_mortar"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/small_2d_plate", "cosim_fem_fem_small_2d_plate_full_mortar")
            self._runTest()

    #def test_FEM_FEM_dynamic_2d_cantilever_implicit_implicit(self):
    #    if not numpy_available:
    #        self.skipTest("Numpy not available")
    #    if not have_fem_fem_dependencies:
    #        self.skipTest("FEM-FEM dependencies are not available!")
    #
    #    self.name = "test_FEM_FEM_dynamic_2d_cantilever_implicit_implicit"
    #    with KratosUnittest.WorkFolderScope(".", __file__):
    #        self._createTest("fem_fem/dynamic_2d_cantilever/implicit_implicit", "fem_fem_dynamic_2d_cantilever")
    #        self._runTest()
    #
    #def test_FEM_FEM_dynamic_2d_cantilever_implicit_implicit_nonconforming(self):
    #    if not numpy_available:
    #        self.skipTest("Numpy not available")
    #    if not have_fem_fem_dependencies:
    #        self.skipTest("FEM-FEM dependencies are not available!")
    #
    #    self.name = "test_FEM_FEM_dynamic_2d_cantilever_implicit_implicit_nonconforming"
    #    with KratosUnittest.WorkFolderScope(".", __file__):
    #        self._createTest("fem_fem/dynamic_2d_cantilever/implicit_implicit", "fem_fem_dynamic_2d_cantilever_nonconforming")
    #        self._runTest()
    #
    #def test_FEM_FEM_dynamic_2d_cantilever_implicit_implicit_mixed_timestep(self):
    #    if not numpy_available:
    #        self.skipTest("Numpy not available")
    #    if not have_fem_fem_dependencies:
    #        self.skipTest("FEM-FEM dependencies are not available!")
    #
    #    self.name = "test_FEM_FEM_dynamic_2d_cantilever_implicit_implicit_mixed_timestep"
    #    with KratosUnittest.WorkFolderScope(".", __file__):
    #        self._createTest("fem_fem/dynamic_2d_cantilever/implicit_implicit", "fem_fem_dynamic_2d_cantilever_mixed_timestep")
    #        self._runTest()
    #
    #def test_FEM_FEM_dynamic_2d_cantilever_implicit_explicit_mixed_nonconforming(self):
    #    if not numpy_available:
    #        self.skipTest("Numpy not available")
    #    if not have_fem_fem_dependencies:
    #        self.skipTest("FEM-FEM dependencies are not available!")
    #
    #    self.name = "test_FEM_FEM_dynamic_2d_cantilever_implicit_explicit_mixed_nonconforming"
    #    with KratosUnittest.WorkFolderScope(".", __file__):
    #        self._createTest("fem_fem/dynamic_2d_cantilever/implicit_explicit", "fem_fem_dynamic_2d_cantilever_mixed_nonconforming")
    #        self._runTest()

    def test_sdof_static_fsi(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_potential_fsi_dependencies:
            self.skipTest("FSI dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_sdof_static", "project_cosim_naca0012_small_fsi")
            # self.__AddVtkOutputToCFD() # uncomment to get output
            self._runTest()

class TestCoSimulationCases(co_simulation_test_case.CoSimulationTestCase):
    '''This class contains "full" CoSimulation-Cases, too large for the nightly suite and therefore
    have to be in the validation-suite
    '''
    def test_WallFSI(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_fsi_dependencies:
            self.skipTest("FSI dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_wall", "cosim_wall_weak_coupling_fsi")
            self._runTest()

    def test_DEMFEMCableNet(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_dem_fem_dependencies:
            self.skipTest("DEM FEM dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("dem_fem_cable_net","cosim_dem_fem_cable_net")
            self._runTest()

    def test_sdof_fsi(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if not have_fsi_dependencies:
            self.skipTest("FSI dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_sdof", "cosim_sdof_fsi")
            # self.__AddVtkOutputToCFD() # uncomment to get output
            self._runTest()

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()

        # delete superfluous dem files
        kratos_utils.DeleteFileIfExisting(GetFilePath("dem_fem_cable_net/cableNet.post.lst"))
        kratos_utils.DeleteDirectoryIfExisting(GetFilePath("dem_fem_cable_net/cableNet_Graphs"))
        kratos_utils.DeleteDirectoryIfExisting(GetFilePath("dem_fem_cable_net/cableNet_MPI_results"))
        kratos_utils.DeleteDirectoryIfExisting(GetFilePath("dem_fem_cable_net/cableNet_Post_Files"))
        kratos_utils.DeleteDirectoryIfExisting(GetFilePath("dem_fem_cable_net/cableNet_Results_and_Data"))

if __name__ == '__main__':
    KratosUnittest.main()
