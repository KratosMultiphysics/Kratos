import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

data_comm = KM.Testing.GetDefaultDataCommunicator()

import co_simulation_test_case
import os

have_fsi_dependencies = kratos_utils.CheckIfApplicationsAvailable("FluidDynamicsApplication", "StructuralMechanicsApplication", "MappingApplication", "MeshMovingApplication", "LinearSolversApplication")
have_potential_fsi_dependencies = kratos_utils.CheckIfApplicationsAvailable("CompressiblePotentialFlowApplication", "StructuralMechanicsApplication", "MappingApplication", "MeshMovingApplication", "LinearSolversApplication")
have_mpm_fem_dependencies = kratos_utils.CheckIfApplicationsAvailable("MPMApplication", "StructuralMechanicsApplication", "MappingApplication", "LinearSolversApplication", "ConstitutiveLawsApplication")
have_dem_fem_dependencies = kratos_utils.CheckIfApplicationsAvailable("DEMApplication", "StructuralMechanicsApplication", "MappingApplication", "LinearSolversApplication")
have_mpm_dem_dependencies = kratos_utils.CheckIfApplicationsAvailable("DEMApplication", "MPMApplication", "MappingApplication", "LinearSolversApplication")
have_fem_fem_dependencies = kratos_utils.CheckIfApplicationsAvailable("StructuralMechanicsApplication", "MappingApplication")
have_pfem_fem_dependencies = kratos_utils.CheckIfApplicationsAvailable("PfemFluidDynamicsApplication", "StructuralMechanicsApplication", "MappingApplication", "LinearSolversApplication", "ConstitutiveLawsApplication")

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestTinyFetiCoSimulationCases(co_simulation_test_case.CoSimulationTestCase):
    '''This class contains "tiny" FETI CoSimulation-Cases, small enough to run in the CI
    '''
    def test_FEM_FEM_small_2d_plate_feti_explict_explicit(self):
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_small_2d_plate_feti_explict_explicit"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/small_2d_plate_feti/explicit_explicit", "cosim_fem_fem_small_2d_plate_feti_explicit_explicit")
            self._runTest()

    def test_FEM_FEM_small_2d_plate_feti_implict_explicit(self):
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_small_2d_plate_feti_implict_explicit"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/small_2d_plate_feti/implicit_explicit", "cosim_fem_fem_small_2d_plate_feti_implicit_explicit")
            self._runTest()

    def test_FEM_FEM_small_2d_plate_feti_implict_implicit(self):
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_small_2d_plate_feti_implict_implicit"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/small_2d_plate_feti/implicit_implicit", "cosim_fem_fem_small_2d_plate_feti_implicit_implicit")
            self._runTest()

    def test_FEM_FEM_small_2d_plate_feti_implict_explicit_mixed(self):
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_small_2d_plate_feti_implict_explicit_mixed"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/small_2d_plate_feti/implicit_explicit_mixed", "cosim_fem_fem_small_2d_plate_feti_implicit_explicit_mixed")
            self._runTest()

    def test_MPMDEMCoupling(self):
        if not have_mpm_dem_dependencies:
            self.skipTest("MPM DEM dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("mpm_dem","cosim_mpm_dem")
            self._runTest()

        # removing superfluous dem files after test
        self.addCleanup(kratos_utils.DeleteFileIfExisting, GetFilePath("mpm_dem/dempart.post.lst"))
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, GetFilePath("mpm_dem/dempart_Graphs"))
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, GetFilePath("mpm_dem/dempart_MPI_results"))
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, GetFilePath("mpm_dem/dempart_Post_Files"))
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, GetFilePath("mpm_dem/dempart_Results_and_Data"))


class TestSmallCoSimulationCases(co_simulation_test_case.CoSimulationTestCase):
    '''This class contains "small" CoSimulation-Cases, small enough to run in the nightly suite
    '''

    def test_FEM_FEM_small_2d_plate_dual_mortar(self):
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_small_2d_plate_dual_mortar"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/small_2d_plate", "cosim_fem_fem_small_2d_plate_dual_mortar")
            self._runTest()

    def test_FEM_FEM_small_2d_plate_full_mortar(self):
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_small_2d_plate_full_mortar"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/small_2d_plate", "cosim_fem_fem_small_2d_plate_full_mortar")
            self._runTest()

    def test_FEM_FEM_Neumann_Neumann_Jacobi_Solver(self):
        if not have_fem_fem_dependencies:
            self.skipTest("FEM-FEM dependencies are not available!")

        self.name = "test_FEM_FEM_neumann_neumann_jacobi_solver"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fem_fem/static_2d_cantilever", "cosim_fem_fem_neumann_neumann_jacobi_solver")
            self._runTest()

    #def test_FEM_FEM_dynamic_2d_cantilever_implicit_implicit(self):
    #    if not have_fem_fem_dependencies:
    #        self.skipTest("FEM-FEM dependencies are not available!")
    #
    #    self.name = "test_FEM_FEM_dynamic_2d_cantilever_implicit_implicit"
    #    with KratosUnittest.WorkFolderScope(".", __file__):
    #        self._createTest("fem_fem/dynamic_2d_cantilever/implicit_implicit", "fem_fem_dynamic_2d_cantilever")
    #        self._runTest()
    #
    #def test_FEM_FEM_dynamic_2d_cantilever_implicit_implicit_nonconforming(self):
    #    if not have_fem_fem_dependencies:
    #        self.skipTest("FEM-FEM dependencies are not available!")
    #
    #    self.name = "test_FEM_FEM_dynamic_2d_cantilever_implicit_implicit_nonconforming"
    #    with KratosUnittest.WorkFolderScope(".", __file__):
    #        self._createTest("fem_fem/dynamic_2d_cantilever/implicit_implicit", "fem_fem_dynamic_2d_cantilever_nonconforming")
    #        self._runTest()
    #
    #def test_FEM_FEM_dynamic_2d_cantilever_implicit_implicit_mixed_timestep(self):
    #    if not have_fem_fem_dependencies:
    #        self.skipTest("FEM-FEM dependencies are not available!")
    #
    #    self.name = "test_FEM_FEM_dynamic_2d_cantilever_implicit_implicit_mixed_timestep"
    #    with KratosUnittest.WorkFolderScope(".", __file__):
    #        self._createTest("fem_fem/dynamic_2d_cantilever/implicit_implicit", "fem_fem_dynamic_2d_cantilever_mixed_timestep")
    #        self._runTest()
    #
    #def test_FEM_FEM_dynamic_2d_cantilever_implicit_explicit_mixed_nonconforming(self):
    #    if not have_fem_fem_dependencies:
    #        self.skipTest("FEM-FEM dependencies are not available!")
    #
    #    self.name = "test_FEM_FEM_dynamic_2d_cantilever_implicit_explicit_mixed_nonconforming"
    #    with KratosUnittest.WorkFolderScope(".", __file__):
    #        self._createTest("fem_fem/dynamic_2d_cantilever/implicit_explicit", "fem_fem_dynamic_2d_cantilever_mixed_nonconforming")
    #        self._runTest()

    def test_sdof_static_fsi(self):
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

    def test_MPM_FEM_beam_penalty(self):
        if not have_mpm_fem_dependencies:
            self.skipTest("MPM-FEM dependencies are not available!")

        self.name = "penalty_beam"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("mpm_fem_beam", "cosim_mpm_fem_beam")
            self._runTest()

    def test_WallFSI(self):
        if not have_fsi_dependencies:
            self.skipTest("FSI dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_wall", "cosim_wall_weak_coupling_fsi")
            self._runTest()

    def test_DEMFEMCableNet(self):
        if not have_dem_fem_dependencies:
            self.skipTest("DEM FEM dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("dem_fem_cable_net","cosim_dem_fem_cable_net")
            self._runTest()

        # removing superfluous dem files after test
        self.addCleanup(kratos_utils.DeleteFileIfExisting, GetFilePath("dem_fem_cable_net/cableNet.post.lst"))
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, GetFilePath("dem_fem_cable_net/cableNet_Graphs"))
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, GetFilePath("dem_fem_cable_net/cableNet_MPI_results"))
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, GetFilePath("dem_fem_cable_net/cableNet_Post_Files"))
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, GetFilePath("dem_fem_cable_net/cableNet_Results_and_Data"))

    def test_sdof_fsi(self):
        if not have_fsi_dependencies:
            self.skipTest("FSI dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_sdof", "cosim_sdof_fsi")
            if data_comm.IsDistributed():
                allowed_num_processes = [3,4,5] # problem is small and needs a very specific number of processes, otherwise the linear solver gives slightly different results and the test fails
                if data_comm.Size() not in allowed_num_processes:
                    self.skipTest("This test runs only with {} processes".format(allowed_num_processes))

                self.cosim_parameters["problem_data"]["parallel_type"].SetString("MPI")
                fluid_solver_settings = self.cosim_parameters["solver_settings"]["solvers"]["fluid"]["solver_wrapper_settings"]
                fluid_solver_settings["input_file"].SetString("fsi_sdof/ProjectParametersCFD_mpi") # TODO refactor such that serial file can be reused. Requires to update and dump new CFD settings (similar to mok test)

            self._runTest()

    @KratosUnittest.skipUnless(False, "this test result is not evaluated")
    def test_PFEM_FEM_water_slide_2d(self):
        if not have_pfem_fem_dependencies:
            self.skipTest("PFEM FEM dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("pfem_fem_waterslide2d","cosim_pfem_fem_waterslide2d")
            self._runTest()


if __name__ == '__main__':
    KratosUnittest.main()
