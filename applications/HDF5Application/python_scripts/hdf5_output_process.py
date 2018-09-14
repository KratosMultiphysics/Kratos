from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5

# Other imports
import os
try:
    # in case the h5py-module is not installed (e.g. on clusters) we don't want it to crash the simulation!
    # => in such a case the xdmf can be created manually afterwards locally
    import create_xdmf_file
    have_xdmf = True
except ImportError:
    have_xdmf = False


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return HDF5OutputProcess(Model, settings["Parameters"])


class HDF5OutputProcess(KratosMultiphysics.Process):
    """This process writes output in hdf5-files
    Additionally it can create xdmf-files based on the hdf5-files which can be used
    to visualize the files in postptocessing-tools
    """

    def __init__(self, Model, settings):
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"                    : "PLEASE_SPECIFY_MOEL_PART_NAME",
            "file_name"                          : "",
            "save_h5_files_in_folder"            : true,
            "backwards_in_time"                  : false,
            "create_xdmf_file_level"             : 1,
            "nodal_solution_step_data_variables" : [],
            "nodal_data_value_variables"         : [],
            "element_data_value_variables"       : []
        }
        """)

        self.model = Model

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part_name = self.settings["model_part_name"].GetString()
        self.file_name = self.settings["file_name"].GetString()
        if self.file_name == "": # in case not file name is specified then the name of the ModelPart is used
            self.file_name = self.model_part_name

        self.raw_path, self.file_name = os.path.split(self.file_name) # maybe add os.getcwd to raw_path?
        self.raw_path = os.path.join(os.getcwd(), self.raw_path)

        self.folder_name = self.file_name + "__h5_files"
        self.save_h5_files_in_folder = self.settings["save_h5_files_in_folder"].GetBool()

        self.create_xdmf_file_level = self.settings["create_xdmf_file_level"].GetInt()
        if self.create_xdmf_file_level > 0 and have_xdmf is False:
            KratosMultiphysics.Logger.PrintWarning("XDMF-Writing is not available!")
            self.create_xdmf_file_level = 0
        # self.xdmf_file_name = self.__GetFullFilePath() + ".xdmf"
        self.num_output_files = 0 # force rewritting xdmf

    def ExecuteInitialize(self):
        # create folder if necessary and adapt paths
        hfd5_writer_process_parameters = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name" : "%s",
                "nodal_solution_step_data_settings" : { },
                "nodal_data_value_settings"         : { },
                "element_data_value_settings"       : { },
                "file_settings" : {
                    "file_access_mode" : "truncate"
                },
                "output_time_settings" : {
                    "output_step_frequency" : 20,
                    "file_name" : "%s"
                }
            }
        }""" % (self.model_part_name, self.__GetFullFilePath(self.file_name)))

        hdf5_params = hfd5_writer_process_parameters["Parameters"]

        hdf5_params["nodal_solution_step_data_settings"].AddValue(
            "list_of_variables", self.settings["nodal_solution_step_data_variables"])
        hdf5_params["nodal_data_value_settings"].AddValue(
            "list_of_variables", self.settings["nodal_data_value_variables"])
        hdf5_params["element_data_value_settings"].AddValue(
            "list_of_variables", self.settings["element_data_value_variables"])


        model_part_comm = self.model[self.model_part_name].GetCommunicator()
        is_mpi_execution = (model_part_comm.TotalProcesses() > 1)

        if self.save_h5_files_in_folder:
            folder_path = self.__GetFolderPathSave()
            if not os.path.isdir(folder_path) and model_part_comm.MyPID() == 0:
                os.makedirs(folder_path)
            model_part_comm.Barrier()

        if is_mpi_execution:
            import partitioned_single_mesh_temporal_output_process as hdf5_process
            hdf5_params["file_settings"].AddEmptyValue("file_driver").SetString("mpio") # TODO needed?
        else:
            import single_mesh_temporal_output_process as hdf5_process

        print(hfd5_writer_process_parameters.PrettyPrintJsonString())

        # TODO delete the old files!

        self.hfd5_writer_process = hdf5_process.Factory(hfd5_writer_process_parameters, self.model)
        self.hfd5_writer_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.hfd5_writer_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.hfd5_writer_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.hfd5_writer_process.ExecuteFinalizeSolutionStep()

        # Create xdmf-file
        if self.create_xdmf_file_level > 1:
            current_num_out_files = self.__GetNumberOfOutputFiles()
            # the number of files has changed, the xdmf is rewritten
            if self.num_output_files != current_num_out_files:
                self.num_output_files = current_num_out_files
                self.__WriteXdmfFile()
                print("Recreating XDMF")

    def ExecuteBeforeOutputStep(self):
        self.hfd5_writer_process.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        self.hfd5_writer_process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        self.hfd5_writer_process.ExecuteFinalize()

        if self.create_xdmf_file_level == 1:
            # if create_xdmf_file_level is larger, then the xdmf-file will
            # already have been created in "ExecuteFinalizeSolutionStep"
            self.__WriteXdmfFile()

    def __WriteXdmfFile(self):
        current_time = self.model[self.model_part_name].ProcessInfo[KratosMultiphysics.TIME]
        create_xdmf_file.WriteXdmfFile(self.file_name + ".h5",
                                       self.__GetFolderPathSave(),
                                       current_time)

    def __GetNumberOfOutputFiles(self):
        h5_files_path = os.path.join(self.__GetFolderPathSave(), self.file_name + ".h5")
        current_time = self.model[self.model_part_name].ProcessInfo[KratosMultiphysics.TIME]
        list_time_labels = create_xdmf_file.GetListOfTimeLabels(h5_files_path)

        num_output_files = sum(float(time_label) <= current_time for time_label in list_time_labels)

        return num_output_files

    def __GetFolderPathSave(self):
        if self.save_h5_files_in_folder:
            return os.path.join(self.raw_path, self.folder_name)
        else:
            return self.raw_path

    def __GetFullFilePath(self, file_name):
        return os.path.join(self.__GetFolderPathSave(), file_name)