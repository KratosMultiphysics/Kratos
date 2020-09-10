from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.factories.convergence_criterion_factory import CreateConvergenceCriterion
import KratosMultiphysics.CoSimulationApplication.colors as colors

# Other imports
import numpy as np

class ConvergenceCriteriaWrapper(object):
    """This class wraps the convergence criteria such that they can be used "automized"
    => this class stores the residual and updates the solutions, such that the
    convergence criteria can be configured through json
    In case of distributed data, the data is gathered on one rank, the convergence checked and the result broadcasted to the other ranks
    """
    def __init__(self, settings, solver_wrapper):
        self.interface_data = solver_wrapper.GetInterfaceData(settings["data_name"].GetString())
        settings.RemoveValue("data_name")
        settings.RemoveValue("solver")

        if not settings.Has("label"):
            settings.AddEmptyValue("label").SetString(colors.bold('{}.{}'.format(self.interface_data.solver_name, self.interface_data.name)))

        self.conv_crit = CreateConvergenceCriterion(settings)

    def Initialize(self):
        self.conv_crit.Initialize()

    def Finalize(self):
        self.conv_crit.Finalize()

    def InitializeSolutionStep(self):
        self.conv_crit.InitializeSolutionStep()

        # MPI related - TODO might be better to do one in Initialize, but the InterfaceData is not yet initialized there yet (might be possible in the future)
        self.executing_rank = (self.interface_data.GetModelPart().GetCommunicator().MyPID() == 0)
        if self.interface_data.IsDistributed():
            self.data_comm = self.interface_data.GetModelPart().GetCommunicator().GetDataCommunicator()

    def FinalizeSolutionStep(self):
        self.conv_crit.FinalizeSolutionStep()

    def InitializeNonLinearIteration(self):
        # Saving the previous data for the computation of the residual
        # and the computation of the solution update
        self.input_data = self.interface_data.GetData()

        self.conv_crit.InitializeNonLinearIteration()

    def FinalizeNonLinearIteration(self):
        self.conv_crit.FinalizeNonLinearIteration()

    def IsConverged(self):
        current_data = self.interface_data.GetData()
        residual = current_data - self.input_data

        if self.interface_data.IsDistributed():
            residual = np.array(np.concatenate(self.data_comm.GathervDoubles(residual, 0)))
            current_data = np.array(np.concatenate(self.data_comm.GathervDoubles(current_data, 0)))

        is_converged = 0
        if self.executing_rank:
            is_converged = self.conv_crit.IsConverged(residual, current_data)

        if self.interface_data.IsDistributed():
            is_converged = self.data_comm.Broadcast(is_converged, 0)

        return is_converged

    def PrintInfo(self):
        self.conv_crit.PrintInfo()

    def Check(self):
        self.conv_crit.Check()
