from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

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
        super(DistributedRestartUtility, self).__init__(model_part, settings)

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

    def __ExtractFileLabel(self, file_name):
        label_begin = file_name.find(self.raw_file_name + '_') + len(self.raw_file_name + '_')
        label_end   = file_name.find('.rest')
        labels      = file_name[label_begin:label_end].split('_')
        return labels