# Core imports
import KratosMultiphysics

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.factories.convergence_criterion_factory import CreateConvergenceCriterion
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
import KratosMultiphysics.CoSimulationApplication.colors as colors

# Other imports
import numpy as np

class ConvergenceCriteriaWrapper:
    """ @brief This class wraps the convergence criteria such that they can be used "automated".
        @details This class stores the residual and updates the solutions, such that the
                 convergence criteria can be configured through JSON.
                 In case of distributed data, the data is gathered on one rank, the convergence
                 checked and the result broadcast to the other ranks.
    """
    def __init__(self,
                 settings: KratosMultiphysics.Parameters,
                 interface_data: CouplingInterfaceData,
                 parent_coupled_solver_data_communicator: KratosMultiphysics.DataCommunicator):
        self.interface_data = interface_data

        for key in ("data_name", "solver"):
            if settings.Has(key):
                settings.RemoveValue(key)

        if not settings.Has("label"):
            settings.AddEmptyValue("label").SetString(colors.bold('{}.{}'.format(self.interface_data.solver_name, self.interface_data.name)))

        self.conv_crit = CreateConvergenceCriterion(settings)
        self.data_communicator = parent_coupled_solver_data_communicator

        self.executing_rank = False
        if self.interface_data.IsDefinedOnThisRank():
            self.data_comm = self.interface_data.GetModelPart().GetCommunicator().GetDataCommunicator()
            self.executing_rank = (self.data_comm.Rank() == 0)

    def Initialize(self):
        self.conv_crit.Initialize()

    def Finalize(self):
        self.conv_crit.Finalize()

    def InitializeSolutionStep(self):
        self.conv_crit.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.conv_crit.FinalizeSolutionStep()

    def InitializeNonLinearIteration(self):
        if self.interface_data.IsDefinedOnThisRank():
            # Saving the previous data for the computation of the residual
            # and the computation of the solution update
            self.input_data = self.interface_data.GetData()

        self.conv_crit.InitializeNonLinearIteration()

    def FinalizeNonLinearIteration(self):
        self.conv_crit.FinalizeNonLinearIteration()

    def IsConverged(self):
        if self.interface_data.IsDefinedOnThisRank():
            current_data = self.interface_data.GetData()
            residual = current_data - self.input_data

            if self.interface_data.IsDistributed():
                residual = np.array(np.concatenate(self.data_comm.GathervDoubles(residual, 0)))
                current_data = np.array(np.concatenate(self.data_comm.GathervDoubles(current_data, 0)))

        is_converged = 0
        if self.executing_rank:
            is_converged = self.conv_crit.IsConverged(residual, current_data)

        # all ranks of the coupled solver need to know the convergence information
        is_converged = bool(self.data_communicator.Broadcast(bool(is_converged), 0))

        return is_converged

    def PrintInfo(self):
        self.conv_crit.PrintInfo()

    def Check(self):
        self.conv_crit.Check()
