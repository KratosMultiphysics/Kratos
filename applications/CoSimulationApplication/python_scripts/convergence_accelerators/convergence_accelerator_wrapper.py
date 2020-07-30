from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.factories.convergence_accelerator_factory import CreateConvergenceAccelerator

# Other imports
import numpy as np

class ConvergenceAcceleratorWrapper(object):
    """This class wraps the convergence accelerators such that they can be used "automized"
    => this class stores the residual and updates the solutions, such that the
    convergence accelerator can be configured through json
    In case of distributed data, it is checked whether the convergence accelerator supports it.
    If not, the data is gathered / scattered and the accelerator is executed on only one rank
    """
    def __init__(self, settings, solver_wrapper):
        self.interface_data = solver_wrapper.GetInterfaceData(settings["data_name"].GetString())
        settings.RemoveValue("data_name")
        settings.RemoveValue("solver")

        self.conv_acc = CreateConvergenceAccelerator(settings)

    def Initialize(self):
        self.conv_acc.Initialize()

    def Finalize(self):
        self.conv_acc.Finalize()

    def InitializeSolutionStep(self):
        self.conv_acc.InitializeSolutionStep()

        # MPI related - TODO might be better to do one in Initialize, but the InterfaceData is not yet initialized there yet (might be possible in the future)
        # However if this is done in initialize, then we would have to Clear or sth in order to make it work with Remeshing (or if the sizes change for other reasons) ...
        conv_acc_supports_dist_data = self.conv_acc.SupportsDistributedData()
        self.executing_rank = conv_acc_supports_dist_data or (self.interface_data.GetModelPart().GetCommunicator().MyPID() == 0)
        self.gather_scatter_required = self.interface_data.IsDistributed() and not conv_acc_supports_dist_data
        if self.gather_scatter_required:
            self.data_comm = self.interface_data.GetModelPart().GetCommunicator().GetDataCommunicator()
            self.sizes_from_ranks = np.cumsum(self.data_comm.GatherInts([self.interface_data.Size()], 0))

    def FinalizeSolutionStep(self):
        self.conv_acc.FinalizeSolutionStep()

    def InitializeNonLinearIteration(self):
        # Saving the previous data for the computation of the residual
        # and the computation of the solution update
        self.input_data = self.interface_data.GetData()

        self.conv_acc.InitializeNonLinearIteration()

    def FinalizeNonLinearIteration(self):
        self.conv_acc.FinalizeNonLinearIteration()

    def ComputeAndApplyUpdate(self):
        current_data = self.interface_data.GetData()
        residual = current_data - self.input_data
        input_data_for_acc = self.input_data

        if self.gather_scatter_required:
            residual = np.array(np.concatenate(self.data_comm.GathervDoubles(residual, 0)))
            input_data_for_acc = np.array(np.concatenate(self.data_comm.GathervDoubles(input_data_for_acc, 0)))

        if self.executing_rank:
            updated_data = input_data_for_acc + self.conv_acc.UpdateSolution(residual, input_data_for_acc)

        if self.gather_scatter_required:
            if self.executing_rank:
                data_to_scatter = np.split(updated_data, self.sizes_from_ranks[:-1])
            else:
                data_to_scatter = []

            updated_data = self.data_comm.ScattervDoubles(data_to_scatter, 0)

        self.interface_data.SetData(updated_data)

    def PrintInfo(self):
        self.conv_acc.PrintInfo()

    def Check(self):
        self.conv_acc.Check()
