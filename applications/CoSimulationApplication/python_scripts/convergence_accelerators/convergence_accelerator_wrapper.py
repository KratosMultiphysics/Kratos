# Core imports
import KratosMultiphysics

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.factories.convergence_accelerator_factory import CreateConvergenceAccelerator
from ..coupling_interface_data import CouplingInterfaceData

# Other imports
import numpy as np
from abc import ABCMeta, abstractmethod

class ConvergenceAcceleratorWrapper:
    """This class wraps the convergence accelerators such that they can be used "automized"
    => this class stores the residual and updates the solutions, such that the
    convergence accelerator can be configured through json
    In case of distributed data, it is checked whether the convergence accelerator supports it.
    If not, the data is gathered / scattered and the accelerator is executed on only one rank
    """
    def __init__(self,
                 settings: KratosMultiphysics.Parameters,
                 interface_data_dict: "dict[str,CouplingInterfaceData]",
                 parent_coupled_solver_data_communicator: KratosMultiphysics.DataCommunicator):
        # self.interface_data = interface_data_dict[settings["data_name"].GetString()]
        # self.residual_computation = CreateResidualComputation(settings, interface_data_dict)
        self.is_two_field_accelerator = settings.Has("data_names")

        if self.is_two_field_accelerator:
            if settings["data_names"].size() != 2:
                raise Exception('"data_names" must contain exactly two entries!')

            self.data_name_1 = settings["data_names"][0].GetString()
            self.data_name_2 = settings["data_names"][1].GetString()

            self.interface_data_1 = interface_data_dict[self.data_name_1]
            self.interface_data_2 = interface_data_dict[self.data_name_2]

            if self.interface_data_1.Size() != self.interface_data_2.Size():
                raise Exception(
                    f'Two-field accelerator requires equal data sizes. '
                    f'{self.data_name_1}: {self.interface_data_1.Size()}, '
                    f'{self.data_name_2}: {self.interface_data_2.Size()}'
                )

            self.interface_data = self.interface_data_1

        else:
            self.interface_data = interface_data_dict[settings["data_name"].GetString()]
            self.residual_computation = CreateResidualComputation(settings, interface_data_dict)

        # Remove extra entries from accelerator parameters
        for key in ("data_name", "data_names", "solver", "residual_computation"):
            if settings.Has(key):
                settings.RemoveValue(key)

        self.conv_acc = CreateConvergenceAccelerator(settings)
        self.data_communicator = parent_coupled_solver_data_communicator

        self.executing_rank = False
        self.gather_scatter_required = False

        if self.interface_data.IsDefinedOnThisRank():
            conv_acc_supports_dist_data = self.conv_acc.SupportsDistributedData()
            self.executing_rank = conv_acc_supports_dist_data or (
                self.interface_data.GetModelPart().GetCommunicator().MyPID() == 0
            )
            self.gather_scatter_required = (
                self.interface_data.IsDistributed() and not conv_acc_supports_dist_data
            )

            if self.gather_scatter_required:
                self.data_comm = self.interface_data.GetModelPart().GetCommunicator().GetDataCommunicator()
                self.sizes_from_ranks = np.cumsum(
                    self.data_comm.GatherInts([self.interface_data.Size()], 0)
                )

    def Initialize(self):
        self.conv_acc.Initialize()

    def Finalize(self):
        self.conv_acc.Finalize()

    def InitializeSolutionStep(self):
        self.conv_acc.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.conv_acc.FinalizeSolutionStep()

    def InitializeNonLinearIteration(self):
        # if self.interface_data.IsDefinedOnThisRank():
        #     # Saving the previous data for the computation of the residual
        #     # and the computation of the solution update
        #     self.input_data = self.interface_data.GetData()
        if self.interface_data.IsDefinedOnThisRank():
            if self.is_two_field_accelerator:
                self.input_data_1 = self.interface_data_1.GetData()
                self.input_data_2 = self.interface_data_2.GetData()
            else:
                self.input_data = self.interface_data.GetData()

        self.conv_acc.InitializeNonLinearIteration()

    def FinalizeNonLinearIteration(self):
        self.conv_acc.FinalizeNonLinearIteration()

    def ComputeAndApplyUpdate(self):
        if self.is_two_field_accelerator:
            if not self.interface_data_1.IsDefinedOnThisRank():
                return

            residual_1 = self.interface_data_1.GetData() - self.input_data_1
            residual_2 = self.interface_data_2.GetData() - self.input_data_2

            input_data_1_for_acc = self.input_data_1
            input_data_2_for_acc = self.input_data_2

            if self.gather_scatter_required:
                residual_1 = np.array(np.concatenate(self.data_comm.GathervDoubles(residual_1, 0)))
                residual_2 = np.array(np.concatenate(self.data_comm.GathervDoubles(residual_2, 0)))

                input_data_1_for_acc = np.array(np.concatenate(self.data_comm.GathervDoubles(input_data_1_for_acc, 0)))
                input_data_2_for_acc = np.array(np.concatenate(self.data_comm.GathervDoubles(input_data_2_for_acc, 0)))

            if self.executing_rank:
                delta_1, delta_2 = self.conv_acc.UpdateSolution(
                    (residual_1, residual_2),
                    (input_data_1_for_acc, input_data_2_for_acc)
                )

                updated_data_1 = input_data_1_for_acc + delta_1
                updated_data_2 = input_data_2_for_acc + delta_2

            if self.gather_scatter_required:
                if self.executing_rank:
                    data_to_scatter_1 = np.split(updated_data_1, self.sizes_from_ranks[:-1])
                    data_to_scatter_2 = np.split(updated_data_2, self.sizes_from_ranks[:-1])
                else:
                    data_to_scatter_1 = []
                    data_to_scatter_2 = []

                updated_data_1 = self.data_comm.ScattervDoubles(data_to_scatter_1, 0)
                updated_data_2 = self.data_comm.ScattervDoubles(data_to_scatter_2, 0)

            self.interface_data_1.SetData(updated_data_1)
            self.interface_data_2.SetData(updated_data_2)

            return

        # original single-field implementation continues below
        if not self.interface_data.IsDefinedOnThisRank(): return

        residual = self.residual_computation.ComputeResidual(self.input_data)
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

class BlockConvergenceAcceleratorWrapper:
    """Similar class for wrapping convergence accelerators for block coupling solvers
    """
    def __init__(self,
                 settings: KratosMultiphysics.Parameters,
                 interface_data_dict: "dict[str,CouplingInterfaceData]",
                 parent_coupled_solver_data_communicator: KratosMultiphysics.DataCommunicator):
        self.interface_data_dict = interface_data_dict
        self.current_solver_id = None
        self.coupl_data_names = {}
        self.input_data = {}
        self.is_block_accelerator = True
        for solverData in settings["solver_sequence"].values():
            self.coupl_data_names[solverData["data_name"].GetString()] = solverData["coupled_data_name"].GetString()
        self.residual_computation = CreateBlockResidualComputation(settings, interface_data_dict)

        # Remove extra entries from accelerator parameters
        for key in ("data_name", "solver", "residual_computation"):
            if settings.Has(key):
                settings.RemoveValue(key)

        self.conv_acc = CreateConvergenceAccelerator(settings)
        self.data_communicator = parent_coupled_solver_data_communicator

        self.executing_rank = {}
        self.gather_scatter_required = {}
        self.data_comm = {}
        self.sizes_from_ranks = {}
        for data_name in interface_data_dict:
            if self.interface_data_dict[data_name].IsDefinedOnThisRank():
                conv_acc_supports_dist_data = self.conv_acc.SupportsDistributedData()
                self.executing_rank[data_name] = conv_acc_supports_dist_data or (self.interface_data_dict[data_name].GetModelPart().GetCommunicator().MyPID() == 0)
                self.gather_scatter_required[data_name] = self.interface_data_dict[data_name].IsDistributed() and not conv_acc_supports_dist_data
                if self.gather_scatter_required[data_name]:
                    self.data_comm[data_name] = self.interface_data_dict[data_name].GetModelPart().GetCommunicator().GetDataCommunicator()
                    self.sizes_from_ranks[data_name] = np.cumsum(self.data_comm[data_name].GatherInts([self.interface_data_dict[data_name].Size()], 0))

    def Initialize(self):
        self.conv_acc.Initialize()
        self.current_solver_id = 0
        self.prev_input_data = {}
        self.output_data = {}
        for data_name in self.interface_data_dict:
            if self.interface_data_dict[data_name].IsDefinedOnThisRank():
                self.prev_input_data[data_name] = self.interface_data_dict[data_name].GetData()
                self.output_data[data_name] = self.interface_data_dict[data_name].GetData()

    def Finalize(self):
        self.conv_acc.Finalize()

    def InitializeSolutionStep(self):
        self.conv_acc.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.conv_acc.FinalizeSolutionStep()

    def InitializeNonLinearIteration(self):
        # Saving the previous data for the computation of the residual
        # and the computation of the solution update
        for data_name in self.interface_data_dict:
            interface_data = self.interface_data_dict[data_name]
            if interface_data.IsDefinedOnThisRank():
                self.input_data[data_name] = interface_data.GetData()

        self.conv_acc.InitializeNonLinearIteration()

    def FinalizeNonLinearIteration(self):
        self.conv_acc.FinalizeNonLinearIteration()
        self.current_solver_id = 0

    def ComputeAndApplyUpdate(self):
        # Retrieving solver and data names
        solver_id = self.current_solver_id
        data_name, interface_data = list(self.interface_data_dict.items())[solver_id]
        interface_data = self.interface_data_dict[data_name]
        coupled_data_name = self.coupl_data_names[data_name]

        if not interface_data.IsDefinedOnThisRank(): return

        # Data Comm
        executing_rank = self.executing_rank[data_name]
        gather_scatter_required = self.gather_scatter_required[data_name]

        # Retrieving data
        input_data = self.input_data[data_name]
        self.prev_input_data[data_name] = interface_data.GetData()
        prev_input_data = self.prev_input_data[coupled_data_name]

        residual = self.residual_computation.ComputeResidual(input_data, data_name) #  = x~ - x
        yResidual = self.residual_computation.ComputeResidual(prev_input_data, coupled_data_name) # = y - y~

        input_data_for_acc = input_data
        input_other_solver = self.output_data[coupled_data_name]

        if gather_scatter_required:
            residual = np.array(np.concatenate(self.data_comm[data_name].GathervDoubles(residual, 0)))
            yResidual = np.array(np.concatenate(self.data_comm[coupled_data_name].GathervDoubles(yResidual, 0)))
            input_data_for_acc = np.array(np.concatenate(self.data_comm[data_name].GathervDoubles(input_data_for_acc, 0)))
            input_other_solver = np.array(np.concatenate(self.data_comm[coupled_data_name].GathervDoubles(input_other_solver, 0)))

        if executing_rank:
            updated_data = input_data_for_acc + self.conv_acc.UpdateSolution(residual, input_data_for_acc, input_other_solver, data_name, yResidual)

        if gather_scatter_required:
            if executing_rank:
                data_to_scatter = np.split(updated_data, self.sizes_from_ranks[data_name][:-1])
            else:
                data_to_scatter = []

            updated_data = self.data_comm[data_name].ScattervDoubles(data_to_scatter, 0)

        interface_data.SetData(updated_data)
        self.output_data[data_name] = updated_data
        self.current_solver_id += 1

    def PrintInfo(self):
        self.conv_acc.PrintInfo()

    def Check(self):
        self.conv_acc.Check()

class ConvergenceAcceleratorResidual(metaclass=ABCMeta):
    @abstractmethod
    def ComputeResidual(self, input_data, data_name=None): pass

class DataDifferenceBlockResidual(ConvergenceAcceleratorResidual):
    def __init__(self,
                 settings: KratosMultiphysics.Parameters,
                 interface_data_dict: "dict[str,CouplingInterfaceData]"):

        self.interface_data_dict = interface_data_dict
        self.is_block_residual_computer = True

    def ComputeResidual(self, input_data, data_name=None):
        return self.interface_data_dict[data_name].GetData() - input_data

class DataDifferenceResidual(ConvergenceAcceleratorResidual):
    def __init__(self,
                 settings: KratosMultiphysics.Parameters,
                 interface_data_dict: "dict[str,CouplingInterfaceData]"):
        self.is_block_residual_computer = False
        self.interface_data = interface_data_dict[settings["data_name"].GetString()]

    def ComputeResidual(self, input_data, data_name=None):
        return self.interface_data.GetData() - input_data

class DifferentDataDifferenceResidual(ConvergenceAcceleratorResidual):
    def __init__(self,
                 settings: KratosMultiphysics.Parameters,
                 interface_data_dict: "dict[str,CouplingInterfaceData]"):
        self.interface_data =  interface_data_dict[settings["data_name"].GetString()]
        self.interface_data1 = interface_data_dict[settings["residual_computation"]["data_name1"].GetString()]
        self.interface_data2 = interface_data_dict[settings["residual_computation"]["data_name2"].GetString()]

    def ComputeResidual(self, input_data):
        return self.interface_data1.GetData() - self.interface_data2.GetData()

def CreateBlockResidualComputation(settings: KratosMultiphysics.Parameters,
                                interface_data_dict: "dict[str,CouplingInterfaceData]"):
    residual_computation_type = "data_difference"
    if settings.Has("residual_computation"):
        residual_computation_type = settings["residual_computation"]["type"].GetString()

    if residual_computation_type == "data_difference":
        return DataDifferenceBlockResidual(settings, interface_data_dict)
    else:
        raise Exception(f'The specified residual computation "{residual_computation_type}" is not available!')

def CreateResidualComputation(settings: KratosMultiphysics.Parameters,
                              interface_data_dict: "dict[str,CouplingInterfaceData]"):
    residual_computation_type = "data_difference"
    if settings.Has("residual_computation"):
        residual_computation_type = settings["residual_computation"]["type"].GetString()

    if residual_computation_type == "data_difference":
        return DataDifferenceResidual(settings, interface_data_dict)
    elif residual_computation_type == "different_data_difference":
        return DifferentDataDifferenceResidual(settings, interface_data_dict)
    else:
        raise Exception(f'The specified residual computation "{residual_computation_type}" is not available!')
