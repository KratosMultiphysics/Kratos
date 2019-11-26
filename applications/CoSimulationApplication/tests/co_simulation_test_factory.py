from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import UsingPyKratos
import co_simulation_test_case
import os

try:
    import numpy
    numpy_available = True
except ImportError:
    numpy_available = False

have_fsi_dependencies = kratos_utils.CheckIfApplicationsAvailable("FluidDynamicsApplication", "StructuralMechanicsApplication", "MappingApplication", "MeshMovingApplication", "ExternalSolversApplication")
have_mpm_fem_dependencies = kratos_utils.CheckIfApplicationsAvailable("ParticleMechanicsApplication", "StructuralMechanicsApplication", "MappingApplication", "ExternalSolversApplication")
have_dem_fem_dependencies = kratos_utils.CheckIfApplicationsAvailable("DEMApplication", "StructuralMechanicsApplication", "MappingApplication", "ExternalSolversApplication")

using_pykratos = UsingPyKratos()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestSmallCoSimulationCases(co_simulation_test_case.CoSimulationTestCase):
    '''This class contains "small" CoSimulation-Cases, small enough to run in the nightly suite
    '''

    def test_MokFSI_mvqn(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if using_pykratos:
            self.skipTest("This test cannot be run with pyKratos!")
        if not have_fsi_dependencies:
            self.skipTest("FSI dependencies are not available!")

        self.name = "mvqn"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_mok", "cosim_mok_fsi")
            self.__ManipulateSettings()
            # self.__AddVtkOutputToCFD() # uncomment to get output
            self.__DumpUpdatedCFDSettings()
            self._runTest()

    def test_MokFSI_aitken(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if using_pykratos:
            self.skipTest("This test cannot be run with pyKratos!")
        if not have_fsi_dependencies:
            self.skipTest("FSI dependencies are not available!")

        self.name = "aitken"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_mok", "cosim_mok_fsi")
            self.__ManipulateSettings()
            # self.__AddVtkOutputToCFD() # uncomment to get output
            self.__DumpUpdatedCFDSettings()
            self._runTest()

    def test_MPM_FEM_beam_penalty(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if using_pykratos:
            self.skipTest("This test cannot be run with pyKratos!")
        if not have_mpm_fem_dependencies:
            self.skipTest("MPM-FEM dependencies are not available!")

        self.name = "penalty_beam"
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("mpm_fem_beam", "cosim_mpm_fem_beam")
            self._runTest()

    def __ManipulateSettings(self):

        self.cosim_parameters["solver_settings"]["convergence_accelerators"][0]["type"].SetString(self.name)

        with open("fsi_mok/ProjectParametersCFD_ref.json",'r') as parameter_file:
            self.cfd_parameters = KM.Parameters(parameter_file.read())

        aux_processes = self.cfd_parameters["processes"]["auxiliar_process_list"]
        ref_name = aux_processes[1]["Parameters"]["reference_file_name"].GetString()
        aux_processes[1]["Parameters"]["reference_file_name"].SetString(ref_name[:-4]+"_"+self.name+".dat")
        ref_name = aux_processes[3]["Parameters"]["reference_file_name"].GetString()
        aux_processes[3]["Parameters"]["reference_file_name"].SetString(ref_name[:-4]+"_"+self.name+".dat")


    def __DumpUpdatedCFDSettings(self):
        with open("fsi_mok/ProjectParametersCFD.json", 'w') as parameter_output_file:
            parameter_output_file.write(self.cfd_parameters.PrettyPrintJsonString())

    def __AddVtkOutputToCFD(self):
        self.cfd_parameters["output_processes"].AddValue("vtk_output", KM.Parameters(R'''[{
            "python_module" : "vtk_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "VtkOutputProcess",
            "help"          : "This process writes postprocessing files for Paraview",
            "Parameters"    : {
                "model_part_name"                    : "FluidModelPart.Parts_Fluid",
                "output_control_type"                : "step",
                "output_frequency"                   : 1,
                "file_format"                        : "binary",
                "output_precision"                   : 7,
                "output_sub_model_parts"             : false,
                "folder_name"                        : "fsi_mok/vtk_output_mok_fsi_''' + self.name + '''\",
                "save_output_files_in_folder"        : true,
                "nodal_solution_step_data_variables" : ["VELOCITY","PRESSURE","MESH_DISPLACEMENT","MESH_VELOCITY"],
                "nodal_data_value_variables"         : [],
                "element_data_value_variables"       : [],
                "condition_data_value_variables"     : []
            }
        }]'''))

    @classmethod
    def tearDownClass(cls):
        super(TestSmallCoSimulationCases,cls).tearDownClass()
        kratos_utils.DeleteFileIfExisting(GetFilePath("fsi_mok/ProjectParametersCFD.json"))


class TestCoSimulationCases(co_simulation_test_case.CoSimulationTestCase):
    '''This class contains "full" CoSimulation-Cases, too large for the nightly suite and therefore
    have to be in the validation-suite
    '''
    def test_WallFSI(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if using_pykratos:
            self.skipTest("This test cannot be run with pyKratos!")
        if not have_fsi_dependencies:
            self.skipTest("FSI dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_wall", "cosim_wall_weak_coupling_fsi")
            self._runTest()

    def test_DEMFEMCableNet(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if using_pykratos:
            self.skipTest("This test cannot be run with pyKratos!")
        if not have_dem_fem_dependencies:
            self.skipTest("DEM FEM dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("dem_fem_cable_net","cosim_dem_fem_cable_net")
            self._runTest()

    def test_sdof_fsi(self):
        if not numpy_available:
            self.skipTest("Numpy not available")
        if using_pykratos:
            self.skipTest("This test cannot be run with pyKratos!")
        if not have_fsi_dependencies:
            self.skipTest("FSI dependencies are not available!")

        with KratosUnittest.WorkFolderScope(".", __file__):
            self._createTest("fsi_sdof", "cosim_sdof_fsi")
            # self.__AddVtkOutputToCFD() # uncomment to get output
            self._runTest()

    @classmethod
    def tearDownClass(cls):
        super(TestCoSimulationCases,cls).tearDownClass()

        # delete superfluous dem files
        kratos_utils.DeleteFileIfExisting(GetFilePath("dem_fem_cable_net/cableNet.post.lst"))
        kratos_utils.DeleteDirectoryIfExisting(GetFilePath("dem_fem_cable_net/cableNet_Graphs"))
        kratos_utils.DeleteDirectoryIfExisting(GetFilePath("dem_fem_cable_net/cableNet_MPI_results"))
        kratos_utils.DeleteDirectoryIfExisting(GetFilePath("dem_fem_cable_net/cableNet_Post_Files"))
        kratos_utils.DeleteDirectoryIfExisting(GetFilePath("dem_fem_cable_net/cableNet_Results_and_Data"))

if __name__ == '__main__':
    KratosUnittest.main()
