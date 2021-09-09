# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Other imports
from KratosMultiphysics.restart_utility import RestartUtility

class DistributedRestartUtility(RestartUtility):
    """
    This class overwrites the methods that are different
    in MPI parallel execution

    See restart_utility.py for more information.
    """
    def __init__(self, model_part, settings):
        # Construct the base class
        super().__init__(model_part, settings)

        self.set_mpi_communicator = settings["set_mpi_communicator"].GetBool()
        # the mpi-comm is not set yet, maybe change this once the communicator can be splitted
        self.rank = KratosMultiphysics.DataCommunicator.GetDefault().Rank()

    #### Protected functions ####

    def _GetFileLabelLoad(self):
        return str(self.rank) + '_' + self.input_file_label

    def _GetFileLabelSave(self, file_label):
        return str(self.rank) + '_' + str(file_label)

    def _ExecuteAfterLoad(self):
        if self.set_mpi_communicator:
            KratosMPI.ParallelFillCommunicator(self.main_model_part.GetRootModelPart()).Execute()

    def _GetFileNamePattern( self ):
        """Return the pattern of flags in the file name for FileNameDataCollector."""
        file_name_pattern = "<model_part_name>_<rank>"
        if self.restart_control_type_is_time:
            file_name_pattern += "_<time>"
        else:
            file_name_pattern += "_<step>"
        file_name_pattern += ".rest"
        return file_name_pattern

    def _UpdateRestartFilesMap(self, restart_files_map, step_id, file_name_data):
        if (self.rank == file_name_data.GetRank()):
            restart_files_map[step_id] = file_name_data.GetFileName()