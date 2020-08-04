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
        is_dummy_solver_wrapper = False
        try:
            self.interface_data = [None]*len(solver_wrapper)
        except TypeError:
            self.interface_data = [None]
            is_dummy_solver_wrapper = True
            pass

        for solver_index in range(0,len(self.interface_data)):
            if is_dummy_solver_wrapper:
                self.interface_data[solver_index] = [solver_wrapper.GetInterfaceData(settings["data_name"].GetString())]
            else:
                self.interface_data[solver_index] = [solver_wrapper[solver_index].GetInterfaceData(settings["data_name"].GetString())]
            if settings.Has("criteria_composition"):
                self.criteria_composition =settings["criteria_composition"].GetString()
                if self.criteria_composition == "energy_conjugate":
                    self.interface_data[solver_index].append(solver_wrapper[solver_index].GetInterfaceData(settings["conjugate_data_name"].GetString()))
            else:
                self.criteria_composition = "primal" # default composition is just a single data source

        self.interface_sign = 1.0
        if settings.Has("criteria_options"):
            if "swap_second_domain_data_sign" in settings["criteria_options"].GetStringArray():
                self.interface_sign = -1.0

        settings.RemoveValue("data_name")
        settings.RemoveValue("solver")

        if not settings.Has("label"):
            settings.AddEmptyValue("label").SetString(colors.bold('{}.{}'.format(self.interface_data[0][0].solver_name,
                                                                                 self.interface_data[0][0].name)))

        self.conv_crit = CreateConvergenceCriterion(settings)

    def Initialize(self):
        self.conv_crit.Initialize()

    def Finalize(self):
        self.conv_crit.Finalize()

    def InitializeSolutionStep(self):
        self.conv_crit.InitializeSolutionStep()

        #TODO not sure how to treat this with potentially multiple solvers

        # MPI related - TODO might be better to do one in Initialize, but the InterfaceData is not yet initialized there yet (might be possible in the future)
        self.executing_rank = (self.interface_data[0][0].GetModelPart().GetCommunicator().MyPID() == 0)
        if self.interface_data[0][0].IsDistributed():
            self.data_comm = self.interface_data[0][0].GetModelPart().GetCommunicator().GetDataCommunicator()

    def FinalizeSolutionStep(self):
        self.conv_crit.FinalizeSolutionStep()

    def InitializeNonLinearIteration(self):
        # Saving the previous data for the computation of the residual
        # and the computation of the solution update

        self.input_data = self.GetInterfaceData()
        #self.input_data = self.interface_data[0].GetData()

        self.conv_crit.InitializeNonLinearIteration()

    def FinalizeNonLinearIteration(self):
        self.conv_crit.FinalizeNonLinearIteration()

    def IsConverged(self):
        current_data = self.GetInterfaceData()
        residual = current_data - self.input_data

        #TODO not sure how this handles different solvers
        if self.interface_data[0][0].IsDistributed():
            residual = np.array(np.concatenate(self.data_comm.GathervDoubles(residual, 0)))
            current_data = np.array(np.concatenate(self.data_comm.GathervDoubles(current_data, 0)))

        is_converged = 0
        if self.executing_rank:
            is_converged = self.conv_crit.IsConverged(residual, current_data)

        #TODO not sure how this handles different solvers
        if self.interface_data[0][0].IsDistributed():
            is_converged = self.data_comm.Broadcast(is_converged, 0)

        return is_converged

    def PrintInfo(self):
        self.conv_crit.PrintInfo()

    def Check(self):
        self.conv_crit.Check()

    def GetInterfaceData(self):
        result = 0.0
        for solver_index in range(0,len(self.interface_data)):
            if self.criteria_composition == "primal":
                if solver_index == 0:
                    result = self.interface_data[solver_index][0].GetData()
                else:
                    result -= self.interface_sign*self.interface_data[solver_index][0].GetData() #assumes domain_difference
            elif self.criteria_composition == "energy_conjugate":
                #check length of data vectors are the same
                interface_energy = 0.0;
                data_1 = self.interface_data[solver_index][0].GetData()
                data_2 = self.interface_data[solver_index][1].GetData()
                if len(data_1) != len(data_2):
                    self.__RaiseException('Data vector lengths for conjugate criteria composition must be identical, but they are different!')
                else:
                    for i in range(0,len(data_1)):
                        interface_energy += data_1[i]*data_2[i]
                if solver_index == 0:
                    result = interface_energy
                else:
                    result -= self.interface_sign*interface_energy #assumes domain_difference
            else:
                self.__RaiseException('Invalid criteria_composition specified in cosim parameters json')
        return result
