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

        data_comm_name = settings["data_communicator_name"].GetString()
        if data_comm_name != "":
            self.comm = KratosMultiphysics.ParallelEnvironment.GetDataCommunicator(data_comm_name)
        else:
            self.comm = KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator()

        self.rank = self.comm.Rank()

    #### Protected functions ####

    def _GetFileLabelLoad(self):
        return str(self.rank) + '_' + self.input_file_label

    def _GetFileLabelSave(self, file_label):
        return str(self.rank) + '_' + str(file_label)

    def _ExecuteAfterLoad(self):
        if self.set_mpi_communicator:
            KratosMPI.ParallelFillCommunicator(self.main_model_part.GetRootModelPart(), self.comm).Execute()

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

    def _GetSerializerFlags(self):
        return KratosMultiphysics.Serializer.MPI | KratosMultiphysics.Serializer.SHALLOW_GLOBAL_POINTERS_SERIALIZATION

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "set_mpi_communicator"   : true,
            "data_communicator_name" : ""
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
