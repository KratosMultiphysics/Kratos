from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.factories.convergence_accelerator_factory import CreateConvergenceAccelerator

class ConvergenceAcceleratorWrapper(object):
    """This class wraps the convergence accelerators such that they can be used "automized"
    => this class stores the residual and updates the solutions, such that the
    convergence accelerator can be configured through json
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
        executing_rank = self.conv_acc.SupportsDistributedData() or (self.interface_data.GetModelPart().GetCommunicator().MyPID() == 0)
        gather_scatter_required = self.interface_data.IsDistributed() and not self.conv_acc.SupportsDistributedData()

        current_data = self.interface_data.GetData()
        residual = current_data - self.input_data

        if gather_scatter_required:
            raise NotImplementedError
            # TODO use the DataComm here
            residual_to_use = gather(residual)
            input_data_to_use = gather(self.input_data)
        else:
            residual_to_use = residual
            input_data_to_use = self.input_data

        if executing_rank:
            updated_data = input_data_to_use + self.conv_acc.UpdateSolution(residual_to_use, input_data_to_use)

        if gather_scatter_required:
            raise NotImplementedError
            # TODO use the DataComm here
            updated_data_to_use = scatter(updated_data)
        else:
            updated_data_to_use = updated_data

        self.interface_data.SetData(updated_data_to_use)

    def PrintInfo(self):
        self.conv_acc.PrintInfo()

    def Check(self):
        self.conv_acc.Check()
