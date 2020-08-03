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
        self.criteria_composition = settings["criteria_composition"].GetString()


        #self.interface_data = solver_wrapper.GetInterfaceData(settings["data_name"].GetString())
        self.interface_data = [solver_wrapper.GetInterfaceData(settings["data_name"].GetString())]
        if self.criteria_composition == "energy_conjugate":
            self.interface_data.append(solver_wrapper.GetInterfaceData(settings["additional_data_name"].GetString()))
            print("\n\n\nconjugate\n\n\n")


        settings.RemoveValue("data_name")
        settings.RemoveValue("solver")

        if not settings.Has("label"):
            settings.AddEmptyValue("label").SetString(colors.bold('{}.{}'.format(self.interface_data[0].solver_name, self.interface_data[0].name)))

        self.conv_crit = CreateConvergenceCriterion(settings)

    def Initialize(self):
        self.conv_crit.Initialize()

    def Finalize(self):
        self.conv_crit.Finalize()

    def InitializeSolutionStep(self):
        self.conv_crit.InitializeSolutionStep()

        # MPI related - TODO might be better to do one in Initialize, but the InterfaceData is not yet initialized there yet (might be possible in the future)
        self.executing_rank = (self.interface_data[0].GetModelPart().GetCommunicator().MyPID() == 0)
        if self.interface_data[0].IsDistributed():
            self.data_comm = self.interface_data[0].GetModelPart().GetCommunicator().GetDataCommunicator()

    def FinalizeSolutionStep(self):
        self.conv_crit.FinalizeSolutionStep()

    def InitializeNonLinearIteration(self):
        # Saving the previous data for the computation of the residual
        # and the computation of the solution update

        self.input_data = self.GetData()
        #self.input_data = self.interface_data[0].GetData()

        self.conv_crit.InitializeNonLinearIteration()

    def FinalizeNonLinearIteration(self):
        self.conv_crit.FinalizeNonLinearIteration()

    def IsConverged(self):
        current_data = self.GetData()
        #current_data = self.interface_data.GetData()
        residual = current_data - self.input_data

        if self.interface_data[0].IsDistributed():
            residual = np.array(np.concatenate(self.data_comm.GathervDoubles(residual, 0)))
            current_data = np.array(np.concatenate(self.data_comm.GathervDoubles(current_data, 0)))

        is_converged = 0
        if self.executing_rank:
            is_converged = self.conv_crit.IsConverged(residual, current_data)

        if self.interface_data[0].IsDistributed():
            is_converged = self.data_comm.Broadcast(is_converged, 0)

        return is_converged

    def PrintInfo(self):
        self.conv_crit.PrintInfo()

    def Check(self):
        self.conv_crit.Check()

    def GetData(self):
        if self.criteria_composition == "primal":
            return self.interface_data[0].GetData()
        elif self.criteria_composition == "energy_conjugate":
            #check length of data vectors are the same
            data_1 = self.interface_data[0].GetData()
            data_2 = self.interface_data[1].GetData()
            if len(data_1) != len(data_2):
                self.__RaiseException('Data vector lengths for conjugate criteria composition must be identical but they are different!')
            else:
                scalar_data = 0.0
                for i in range(0,len(data_1)):
                    scalar_data += data_1[i]*data_2[i]
                return scalar_data
        else:
            self.__RaiseException('Invalid criteria_composition specified in cosim parameters json')
