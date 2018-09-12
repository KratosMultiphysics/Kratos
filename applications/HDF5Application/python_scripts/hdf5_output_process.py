from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5

# Other imports
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
            "create_xdmf_file_level"         : 1,
            "save_h5_files_in_folder"        : true,
            "hdf5_writer_process_parameters" : { }
        }
        """)

        self.model = Model

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)
        self.model_part_name = self.settings["hdf5_writer_process_parameters"]["model_part_name"].GetString()
        self.create_xdmf_file_level = self.settings["create_xdmf_file_level"].GetInt()
        if self.create_xdmf_file_level > 0 and have_xdmf is False:
            KratosMultiphysics.Logger.PrintWarning("XDMF-Writing is not available!")
            self.create_xdmf_file_level = 0

    def ExecuteInitialize(self):
        model_part = self.model[self.model_part_name]
        is_mpi_execution = (model_part.GetCommunicator().TotalProcesses() > 1)

        if is_mpi_execution:
            import partitioned_single_mesh_temporal_output_process as hdf5_process
            # TODO set mpiio
        else:
            import single_mesh_temporal_output_process as hdf5_process

        # create folder if necessary and adapt paths

        hfd5_writer_process_parameters = KratosMultiphysics.Parameters("{}")
        hfd5_writer_process_parameters.AddValue("Parameters", self.settings["hdf5_writer_process_parameters"])

        self.hfd5_writer_process = hdf5_process.Factory(hfd5_writer_process_parameters, self.model)
        self.hfd5_writer_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.hfd5_writer_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.hfd5_writer_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.hfd5_writer_process.ExecuteFinalizeSolutionStep()

        # maybe create an xdmf-file every time to be able to constantly visualize it ...?
        # And leave one as backup...?
        # => should be selectable
        # => use the same setting as the one for the hdf5-file itself

        # Create xdmf-file
        if self.create_xdmf_file_level > 1:
            ## TODO also check if a new file was written and only then do this WriteXDMF!
            # => this could be done by checking if len(create_xdmf_file.GetListOfTimeLabels()) has changed
            # in case the h5py-module is not installed (e.g. on clusters) we don't want it to crash the simulation!
            # => in such a case the xdmf can be created manually afterwards locall
            ## Before copying check if the file exists!
            import shutil
            shutil.copyfile(str(self.model_part_name) + ".xdmf", str(self.model_part_name) + "_backup.xdmf")
            ## TODO remove the file again...
            ## for higher lvl keep them and make a counter (Here I could use the counter from above!)
            self.__WriteXdmfFile()

    def ExecuteBeforeOutputStep(self):
        self.hfd5_writer_process.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        self.hfd5_writer_process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        self.hfd5_writer_process.ExecuteFinalize()

        if self.create_xdmf_file_level == 1:
            # if create_xdmf_file_level is larger then the xdmf-file will
            # already have been created in "ExecuteFinalizeSolutionStep"
            self.__WriteXdmfFile()

    def __WriteXdmfFile(self):
        create_xdmf_file.WriteXdmfFile(str(self.model_part_name) + ".h5")
