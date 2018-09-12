from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return HDF5OutputProcess(Model, settings["Parameters"])

class HDF5OutputProcess(KratosMultiphysics.Process):
    """This class is used in order to compute some pre and post process on the SPRISM solid shell elements

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing the settings.
    """

    def __init__(self, Model, settings):
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "create_xdmf_file_level"            : 1,
            "save_restart_files_in_folder"      : true,
            "hdf5_writer_process_parameters"    : { }
        }
        """)

        self.model = Model

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)
        self.model_part_name = self.settings["hdf5_writer_process_parameters"]["model_part_name"].GetString()
        self.create_xdmf_file_level = self.settings["create_xdmf_file_level"].GetInt()

    def ExecuteInitialize(self):
        model_part = self.model[self.model_part_name]
        is_mpi_execution = (model_part.GetCommunicator().TotalProcesses() > 1)

        if is_mpi_execution:
            import partitioned_single_mesh_temporal_output_process as hdf5_process
            # todo set mpiio
        else:
            import single_mesh_temporal_output_process as hdf5_process

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
            # in case the h5py-module is not installed (e.g. on clusters) we don't want it to crash the simulation!
            # => in such a case the xdmf can be created manually afterwards locall
            try:
                import create_xdmf_file
            except ImportError:
                KratosMultiphysics.Logger.PrintWarning("HDF5OutputProcess", "xdmf-file could not be created!")
            create_xdmf_file.main(str(self.model_part_name) + ".h5")

    def ExecuteBeforeOutputStep(self):
        self.hfd5_writer_process.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        self.hfd5_writer_process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        self.hfd5_writer_process.ExecuteFinalize()

        # Create xdmf-file
        if self.create_xdmf_file_level == 1: # if it it larger then it will already have been created in "ExecuteFinalizeSolutionStep"
            # in case the h5py-module is not installed (e.g. on clusters) we don't want it to crash the simulation!
            # => in such a case the xdmf can be created manually afterwards locall
            try:
                import create_xdmf_file
            except ImportError:
                KratosMultiphysics.Logger.PrintWarning("HDF5OutputProcess", "xdmf-file could not be created!")
            create_xdmf_file.main(str(self.model_part_name) + ".h5")
