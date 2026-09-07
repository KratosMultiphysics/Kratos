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
                 interface_data,
                 parent_coupled_solver_data_communicator: KratosMultiphysics.DataCommunicator):
        self.interface_data = interface_data if isinstance(interface_data, list) else [interface_data]

        self.data_combination = "none"
        if settings.Has("data_combination"):
            self.data_combination = settings["data_combination"].GetString()
            settings.RemoveValue("data_combination")

        if self.data_combination == "complex":
            if len(self.interface_data) != 2:
                raise Exception('The "complex" data combination requires exactly two entries in "data_name": real first, imaginary second.')
        elif self.data_combination != "none":
            raise Exception('Unsupported "data_combination": "{}"'.format(self.data_combination))
        elif len(self.interface_data) != 1:
            raise Exception('Multiple entries in "data_name" require a "data_combination".')

        for key in ("data_name", "solver"):
            if settings.Has(key):
                settings.RemoveValue(key)

        if not settings.Has("label"):
            data_labels = [interface_data.name for interface_data in self.interface_data]
            settings.AddEmptyValue("label").SetString(colors.bold('{}.{}'.format(self.interface_data[0].solver_name, '+'.join(data_labels))))

        self.conv_crit = CreateConvergenceCriterion(settings)
        self.data_communicator = parent_coupled_solver_data_communicator

        self.executing_rank = False
        if self.interface_data[0].IsDefinedOnThisRank():
            self.data_comm = self.interface_data[0].GetModelPart().GetCommunicator().GetDataCommunicator()
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
        if self.interface_data[0].IsDefinedOnThisRank():
            # Saving the previous data for the computation of the residual
            # and the computation of the solution update
            self.input_data = self.__GetData()

        self.conv_crit.InitializeNonLinearIteration()

    def FinalizeNonLinearIteration(self):
        self.conv_crit.FinalizeNonLinearIteration()

    def IsConverged(self):
        if self.interface_data[0].IsDefinedOnThisRank():
            current_data = self.__GetData()
            residual = current_data - self.input_data

            if self.interface_data[0].IsDistributed():
                if self.data_combination == "complex":
                    residual = self.__GatherComplexData(residual)
                    current_data = self.__GatherComplexData(current_data)
                else:
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

    def __GetData(self):
        if self.data_combination == "complex":
            real_data = self.interface_data[0].GetData()
            imaginary_data = self.interface_data[1].GetData()
            if real_data.size != imaginary_data.size:
                raise Exception("Real and imaginary interface data must have the same size.")
            return real_data + 1j * imaginary_data

        return self.interface_data[0].GetData()

    def __GatherComplexData(self, data):
        real_data = np.concatenate(self.data_comm.GathervDoubles(data.real, 0))
        imaginary_data = np.concatenate(self.data_comm.GathervDoubles(data.imag, 0))
        return np.array(real_data) + 1j * np.array(imaginary_data)
