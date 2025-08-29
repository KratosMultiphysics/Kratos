import KratosMultiphysics as KM

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.testing.utilities import GetPython3Command

import co_simulation_test_case
import os

have_fsi_dependencies = kratos_utils.CheckIfApplicationsAvailable("FluidDynamicsApplication", "StructuralMechanicsApplication", "MappingApplication", "MeshMovingApplication", "LinearSolversApplication")

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestMokFSI(co_simulation_test_case.CoSimulationTestCase):
    cfd_tes_file_name = "fsi_mok/ProjectParametersCFD_for_test.json"

    def setUp(self):
        self.err_tol = "1e-6"
        if not have_fsi_dependencies:
            self.skipTest("FSI dependencies are not available!")


    def test_mok_fsi_mvqn(self):
        self.accelerator_type = "mvqn"

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_mok", "cosim_mok_fsi")
            self.__ManipulateCFDSettings()
            self.__RemoveOutputFromCFD() # comment to get output
            self.__AddTestingToCFD()
            self.__DumpUpdatedCFDSettings()
            self._runTest()

    def test_mok_fsi_aitken(self):
        self.accelerator_type = "aitken"

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_mok", "cosim_mok_fsi")
            self.__ManipulateCFDSettings()
            self.__RemoveOutputFromCFD() # comment to get output
            self.__AddTestingToCFD()
            self.__DumpUpdatedCFDSettings()
            self._runTest()

    def test_mok_fsi_iqnimvj_explicit(self):
        self.accelerator_type = "iqnimvj_explicit"

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_mok", "cosim_mok_fsi")
            self.__ManipulateCFDSettings()
            self.__RemoveOutputFromCFD() # comment to get output
            self.__AddTestingToCFD()
            self.__DumpUpdatedCFDSettings()
            self._runTest()

    def test_mok_fsi_iqnimvj_implicit(self):
        self.accelerator_type = "iqnimvj_implicit"

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_mok", "cosim_mok_fsi")
            self.__ManipulateCFDSettings()
            self.__RemoveOutputFromCFD() # comment to get output
            self.__AddTestingToCFD()
            self.__DumpUpdatedCFDSettings()
            self._runTest()

    def test_mok_fsi_block_mvqn(self):
        self.accelerator_type = "block_mvqn"

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_mok", "cosim_mok_fsi_block")
            self.__ManipulateCFDSettings()
            self.__RemoveOutputFromCFD() # comment to get output
            self.__AddTestingToCFD()
            self.__DumpUpdatedCFDSettings()
            self._runTest()

    def test_mok_fsi_block_ibqnls(self):
        self.accelerator_type = "block_ibqnls"
        self.err_tol = "6e-5"

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_mok", "cosim_mok_fsi_block")
            self.__ManipulateCFDSettings()
            additional_accel_settings = KM.Parameters("""2""")
            self.cosim_parameters["solver_settings"]["convergence_accelerators"][0].AddValue("timestep_horizon", additional_accel_settings)
            self.__RemoveOutputFromCFD() # comment to get output
            self.__AddTestingToCFD()
            self.__DumpUpdatedCFDSettings()
            self._runTest()

    def test_mok_fsi_mvqn_external_structure(self):
        self.accelerator_type = "mvqn"

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_mok", "cosim_mok_fsi")
            ext_parameter_file_name = os.path.join(self.problem_dir_name, "ProjectParametersCSM.json")
            self.__ManipulateCSMSettingsToRunExternally("solver_wrappers.external.external_solver_wrapper")
            self.__ManipulateCFDSettings()
            self.__RemoveOutputFromCFD() # comment to get output
            self.__AddTestingToCFD()
            self.__DumpUpdatedCFDSettings()
            self._runTestWithExternal([GetPython3Command(), "testing_structural_mechanics_analysis_with_co_sim_io.py", ext_parameter_file_name])

    def test_mok_fsi_mvqn_external_structure_remote_controlled(self):
        self.accelerator_type = "mvqn"

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_mok", "cosim_mok_fsi")
            ext_parameter_file_name = os.path.join(self.problem_dir_name, "ProjectParametersCSM.json")
            self.__ManipulateCSMSettingsToRunExternally("solver_wrappers.external.remote_controlled_solver_wrapper")
            self.__ManipulateCFDSettings()
            self.__RemoveOutputFromCFD() # comment to get output
            self.__AddTestingToCFD()
            self.__DumpUpdatedCFDSettings()
            self._runTestWithExternal([GetPython3Command(), "structural_mechanics_analysis_remote_controlled.py", ext_parameter_file_name])

    @KratosUnittest.skipUnless(KM.IsDistributedRun(), "this test requires MPI")
    def test_mok_fsi_mvqn_external_structure_mpi(self):
        self.accelerator_type = "mvqn"

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_mok", "cosim_mok_fsi")
            ext_parameter_file_name = os.path.join(self.problem_dir_name, "ProjectParametersCSM.json")
            self.__ManipulateCSMSettingsToRunExternally("solver_wrappers.external.external_solver_wrapper")
            self.__ManipulateCFDSettings()
            self.__RemoveOutputFromCFD() # comment to get output
            self.__AddTestingToCFD()
            self.cosim_parameters["problem_data"]["parallel_type"].SetString("MPI")
            self.cfd_parameters["problem_data"]["parallel_type"].SetString("MPI")
            self.__DumpUpdatedCFDSettings()
            num_procs = KM.Testing.GetDefaultDataCommunicator().Size()
            self._runTestWithExternal(["mpiexec", "-np", str(num_procs), GetPython3Command(), "testing_structural_mechanics_analysis_with_co_sim_io.py", "--using-mpi", ext_parameter_file_name])


    def __ManipulateCFDSettings(self):
        self.cosim_parameters["solver_settings"]["convergence_accelerators"][0]["type"].SetString(self.accelerator_type)
        self.cosim_parameters["solver_settings"]["solvers"]["fluid"]["solver_wrapper_settings"]["input_file"].SetString(self.cfd_tes_file_name)

        with open("fsi_mok/ProjectParametersCFD.json",'r') as parameter_file:
            self.cfd_parameters = KM.Parameters(parameter_file.read())

    def __ManipulateCSMSettingsToRunExternally(self, solver_wrapper_name):
        structure_settings = self.cosim_parameters["solver_settings"]["solvers"]["structure"]
        structure_settings.RemoveValue("solver_wrapper_settings")

        structure_settings["type"].SetString(solver_wrapper_name)
        solver_wrapper_settings = KM.Parameters("""{
            "import_meshes" : ["Structure.GENERIC_FSI"],
            "export_data"   : ["load"],
            "import_data"   : ["disp"]
        }""")
        io_settings = KM.Parameters("""{
            "type" : "kratos_co_sim_io",
            "echo_level" : 3,
            "connect_to" : "ext_structure"
        }""")
        structure_settings.AddValue("solver_wrapper_settings", solver_wrapper_settings)
        structure_settings.AddValue("io_settings", io_settings)

    def __AddTestingToCFD(self):
        disp_ref_file_name = "fsi_mok/fsi_mok_cfd_results_disp_ref_{}.dat".format(self.accelerator_type)
        fluid_ref_file_name = "fsi_mok/fsi_mok_cfd_results_fluid_ref_{}.dat".format(self.accelerator_type)

        self.cfd_parameters["processes"].AddValue("testing_processes", KM.Parameters("""[{
            "kratos_module"   : "KratosMultiphysics",
            "python_module"   : "point_output_process",
            "process_name"    : "PointOutputProcess",
            "Parameters" : {
                "position"         : [0.5, 0.25, 0.0],
                "entity_type"      : "node",
                "model_part_name"  : "FluidModelPart",
                "output_file_settings": {
                    "file_name"  : "fsi_mok_cfd_results_disp.dat",
                    "output_path": "fsi_mok"
                },
                "output_variables" : [
                    "MESH_DISPLACEMENT_X",
                    "MESH_DISPLACEMENT_Y",
                    "MESH_VELOCITY_X",
                    "MESH_VELOCITY_Y"]
                }
            },{
            "python_module"   : "compare_two_files_check_process",
            "kratos_module"   : "KratosMultiphysics",
            "process_name"    : "CompareTwoFilesCheckProcess",
            "Parameters" :{
                "output_file_name"    : "fsi_mok/fsi_mok_cfd_results_disp.dat",
                "reference_file_name" : \""""+disp_ref_file_name.replace("\\", "\\\\")+"""\",
                "comparison_type"     : "dat_file_variables_time_history",
                "tolerance"      : """+self.err_tol+"""
                }
            },{
            "kratos_module"   : "KratosMultiphysics",
            "python_module"   : "point_output_process",
            "help"            : "",
            "process_name"    : "PointOutputProcess",
            "Parameters" : {
                "position"         : [0.57, 0.27, 0.0],
                "model_part_name"  : "FluidModelPart",
                "output_file_settings": {
                    "file_name"  : "fsi_mok_cfd_results_fluid.dat",
                    "output_path": "fsi_mok"
                },
                "output_variables" : [
                    "VELOCITY_X",
                    "VELOCITY_Y",
                    "MESH_DISPLACEMENT_X",
                    "MESH_DISPLACEMENT_Y",
                    "MESH_VELOCITY_X",
                    "MESH_VELOCITY_Y"]
                }
            },{
            "python_module"   : "compare_two_files_check_process",
            "kratos_module"   : "KratosMultiphysics",
            "process_name"    : "CompareTwoFilesCheckProcess",
            "Parameters" :{
                "output_file_name"    : "fsi_mok/fsi_mok_cfd_results_fluid.dat",
                "reference_file_name" : \""""+fluid_ref_file_name.replace("\\", "\\\\")+"""\",
                "comparison_type"     : "dat_file_variables_time_history",
                "tolerance"      : """+self.err_tol+"""
                }
            }]"""))

    def __DumpUpdatedCFDSettings(self):
        with open(self.cfd_tes_file_name, 'w') as parameter_output_file:
            parameter_output_file.write(self.cfd_parameters.PrettyPrintJsonString())

    def __RemoveOutputFromCFD(self):
        # by default output is written (in case the same files are used outside the testing environment)
        # for the tests however output is disabled by default
        self.cfd_parameters.RemoveValue("output_processes")

    @classmethod
    def tearDownClass(cls):
        super(TestMokFSI,cls).tearDownClass()
        kratos_utils.DeleteFileIfExisting(GetFilePath(cls.cfd_tes_file_name))


if __name__ == '__main__':
    KratosUnittest.main()
