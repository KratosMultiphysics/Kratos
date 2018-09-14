from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosMultiphysics.kratos_utilities as kratos_utils

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
            "output_control_type"                : "step",
            "output_frequency"                   : 1.0,
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

        self.raw_path, self.file_name = os.path.split(self.file_name)
        # self.raw_path = os.path.join(os.getcwd(), self.raw_path) # maybe add os.getcwd to raw_path?

        self.folder_name = self.file_name + "__h5_files"
        self.save_h5_files_in_folder = self.settings["save_h5_files_in_folder"].GetBool()

        self.create_xdmf_file_level = self.settings["create_xdmf_file_level"].GetInt()
        if self.create_xdmf_file_level > 0 and have_xdmf is False:
            KratosMultiphysics.Logger.PrintWarning("XDMF-Writing is not available!")
            self.create_xdmf_file_level = 0
        self.num_output_files = 0 # force rewritting xdmf


    def ExecuteInitialize(self):
        model_part_comm = self.model[self.model_part_name].GetCommunicator()
        is_mpi_execution = (model_part_comm.TotalProcesses() > 1)

        if self.save_h5_files_in_folder:
            folder_path = self.__GetFolderPathSave()
            if not os.path.isdir(folder_path) and model_part_comm.MyPID() == 0:
                os.makedirs(folder_path)
            model_part_comm.Barrier()

        self.__DeleteOldH5Files()

        hfd5_writer_process_parameters = self.__CreateHDF5WriterProcessParams()

        if is_mpi_execution:
            import partitioned_single_mesh_temporal_output_process as hdf5_process
            file_settings = hfd5_writer_process_parameters["Parameters"]["file_settings"]
            file_settings.AddEmptyValue("file_driver").SetString("mpio") # TODO needed?
        else:
            import single_mesh_temporal_output_process as hdf5_process

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
        create_xdmf_file.WriteXdmfFile(self.file_name + ".h5", self.__GetFolderPathSave())

    def __GetNumberOfOutputFiles(self):
        h5_files_path = os.path.join(self.__GetFolderPathSave(), self.file_name + ".h5")
        list_time_labels = create_xdmf_file.GetListOfTimeLabels(h5_files_path)

        return len(create_xdmf_file.GetListOfTimeLabels(h5_files_path))

    def __GetFolderPathSave(self):
        if self.save_h5_files_in_folder:
            return os.path.join(self.raw_path, self.folder_name)
        else:
            return self.raw_path

    def __GetFullFilePath(self, file_name):
        return os.path.join(self.__GetFolderPathSave(), file_name)

    def __DeleteOldH5Files(self):
        """This function deletes the old h5-files containing the time-data
        """
        kratos_utils.DeleteFileIfExisting(self.file_name+".xdmf")

        current_time = self.model[self.model_part_name].ProcessInfo[KratosMultiphysics.TIME]
        time_prefix = self.file_name+"-"
        backwards_in_time = self.settings["backwards_in_time"].GetBool() # e.g. for adjoint problems
        for name in os.listdir(self.__GetFolderPathSave()):
            if name.find(time_prefix) == 0:
                file_time = float(name.replace(".h5", "")[len(time_prefix):])
                if backwards_in_time:
                    if file_time < current_time:
                        kratos_utils.DeleteFileIfExisting(self.__GetFullFilePath(name))
                else:
                    if file_time > current_time:
                        kratos_utils.DeleteFileIfExisting(self.__GetFullFilePath(name))

    def __CreateHDF5WriterProcessParams(self):
        """This function translates the input for the HDF-output-processes
        """
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

        output_control_type = self.settings["output_control_type"].GetString()
        if output_control_type == "time":
            hdf5_params["output_time_settings"].AddValue("output_time_frequency", self.settings["output_frequency"])
        elif output_control_type == "step":
            # "output_step_frequency" is an int, therefore explicit conversion necessary!
            output_frequency = int(self.settings["output_frequency"].GetDouble())
            hdf5_params["output_time_settings"].AddEmptyValue("output_step_frequency").SetInt(output_frequency)
        else:
            err_msg =  'The requested output_control_type "' + output_control_type + '" is not available!\n'
            err_msg += 'Available options are: "time", "step"'
            raise Exception(err_msg)

        return hfd5_writer_process_parameters