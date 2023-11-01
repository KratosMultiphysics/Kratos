"""!@package HDF5Application

HDF5 core processes.

This module only contains the most general HDF5 IO processes which should not
change frequently.

license: HDF5Application/license.txt
"""
# Kratos imports
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application.core.operations.aggregated_operations import AggregatedControlledOperations

##!@addtogroup HDF5Application
##!@{
##!@name Kratos classes
##!@{

class HDF5Process(KratosMultiphysics.Process):
    __process_local_counter = 0
    __process_mpi_counter = 0

    def __init__(self, model_part: KratosMultiphysics.ModelPart) -> None:
        KratosMultiphysics.Process.__init__(self)
        self.Clear()

        if KratosMultiphysics.IsDistributedRun():
            if model_part.IsDistributed():
                # the run and model parts are distributed.
                HDF5Process.__process_mpi_counter += 1
                self.__process_id = [-1, HDF5Process.__process_mpi_counter]
            else:
                # the run is distributed, but the model part is serial
                HDF5Process.__process_local_counter += 1
                data_comm: KratosMultiphysics.DataCommunicator = model_part.GetCommunicator().GetDataCommunicator()
                self.__process_id = [data_comm.Rank(), HDF5Process.__process_local_counter]
        else:
            # the run is serial, and the model part should be always serial
            HDF5Process.__process_local_counter += 1
            self.__process_id = [-1, HDF5Process.__process_local_counter]

    def GetProcessId(self) -> 'list[int]':
        return self.__process_id

    def Check(self) -> int:
        list(map(lambda x: x.Check(), *self.__aggregated_operations_dict.values()))
        return 0

    def Clear(self) -> None:
        self.__aggregated_operations_dict: 'dict[str, list[AggregatedControlledOperations]]' = {}

    def _AddAggregatedControlledOperations(self, execution_point: str, aggregated_controlled_operation: AggregatedControlledOperations) -> None:
        if not execution_point in self.__aggregated_operations_dict.keys():
            self.__aggregated_operations_dict[execution_point] = []
        self.__aggregated_operations_dict[execution_point].append(aggregated_controlled_operation)

    def _ExecuteAggregatedControlledOperations(self, execution_point: str) -> None:
        # map is evaluated lazily, hence need to do something with return to evaluate the map.
        list(map(lambda x: x.Execute(), self._GetAggregatedControlledOperations(execution_point)))

    def _UpdateAggregatedControlledOperations(self, execution_point: str) -> None:
        # map is evaluated lazily, hence need to do something with return to evaluate the map.
        list(map(lambda x: x.Update(), self._GetAggregatedControlledOperations(execution_point)))

    def _GetAggregatedControlledOperations(self, execution_point: str) -> 'list[AggregatedControlledOperations]':
        if execution_point in self.__aggregated_operations_dict.keys():
            return self.__aggregated_operations_dict[execution_point]
        else:
            return []

    def _GetValidatedParameters(self, sub_parameter_name: str, parameters: KratosMultiphysics.Parameters) -> KratosMultiphysics.Parameters:
        default_sub_parameters = self.GetDefaultParameters()[sub_parameter_name]
        sub_parameters = parameters[sub_parameter_name]
        sub_parameters.ValidateAndAssignDefaults(default_sub_parameters)
        return sub_parameters

    def _GetOperationParameters(self, sub_parameter_name: str, parameters: KratosMultiphysics.Parameters) -> KratosMultiphysics.Parameters:
        sub_parameters = parameters[sub_parameter_name]

        # we do the custom attributes addition first then the validation
        # to make sure that the operation defaults has custom attributes and is capable of writing
        # them.
        if not sub_parameters.Has("custom_attributes"):
            sub_parameters.AddEmptyValue("custom_attributes")

        KratosHDF5.AddProcessId(sub_parameters["custom_attributes"], self.__process_id[0], self.__process_id[1])
        return self._GetValidatedParameters(sub_parameter_name, parameters)

    def AddInitialize(self, aggregated_controlled_operation: AggregatedControlledOperations) -> None:
        self._AddAggregatedControlledOperations("ExecuteInitialize", aggregated_controlled_operation)

    def AddBeforeSolutionLoop(self, aggregated_controlled_operation: AggregatedControlledOperations) -> None:
        self._AddAggregatedControlledOperations("ExecuteBeforeSolutionLoop", aggregated_controlled_operation)

    def AddInitializeSolutionStep(self, aggregated_controlled_operation: AggregatedControlledOperations) -> None:
        self._AddAggregatedControlledOperations("ExecuteInitializeSolutionStep", aggregated_controlled_operation)

    def AddFinalizeSolutionStep(self, aggregated_controlled_operation: AggregatedControlledOperations) -> None:
        self._AddAggregatedControlledOperations("ExecuteFinalizeSolutionStep", aggregated_controlled_operation)

    def AddBeforeOutputStep(self, aggregated_controlled_operation: AggregatedControlledOperations) -> None:
        self._AddAggregatedControlledOperations("ExecuteBeforeOutputStep", aggregated_controlled_operation)

    def AddAfterOutputStep(self, aggregated_controlled_operation: AggregatedControlledOperations) -> None:
        self._AddAggregatedControlledOperations("ExecuteAfterOutputStep", aggregated_controlled_operation)

    def AddFinalize(self, aggregated_controlled_operation: AggregatedControlledOperations) -> None:
        self._AddAggregatedControlledOperations("ExecuteFinalize", aggregated_controlled_operation)

    def ExecuteInitialize(self) -> None:
        self._ExecuteAggregatedControlledOperations("ExecuteInitialize")
        self._UpdateAggregatedControlledOperations("ExecuteInitialize")

    def ExecuteBeforeSolutionLoop(self) -> None:
        self._ExecuteAggregatedControlledOperations("ExecuteBeforeSolutionLoop")
        self._UpdateAggregatedControlledOperations("ExecuteBeforeSolutionLoop")

    def ExecuteInitializeSolutionStep(self) -> None:
        self._ExecuteAggregatedControlledOperations("ExecuteInitializeSolutionStep")
        self._UpdateAggregatedControlledOperations("ExecuteInitializeSolutionStep")

    def ExecuteFinalizeSolutionStep(self) -> None:
        self._ExecuteAggregatedControlledOperations("ExecuteFinalizeSolutionStep")
        self._UpdateAggregatedControlledOperations("ExecuteFinalizeSolutionStep")

    def ExecuteBeforeOutputStep(self) -> None:
        self._ExecuteAggregatedControlledOperations("ExecuteBeforeOutputStep")
        self._UpdateAggregatedControlledOperations("ExecuteBeforeOutputStep")

    def ExecuteAfterOutputStep(self) -> None:
        self._ExecuteAggregatedControlledOperations("ExecuteAfterOutputStep")
        self._UpdateAggregatedControlledOperations("ExecuteAfterOutputStep")

    def ExecuteFinalize(self) -> None:
        self._ExecuteAggregatedControlledOperations("ExecuteFinalize")
        self._UpdateAggregatedControlledOperations("ExecuteFinalize")

class HDF5OutputProcess(KratosMultiphysics.OutputProcess, HDF5Process):
    def __init__(self, model_part: KratosMultiphysics.ModelPart) -> None:
        KratosMultiphysics.OutputProcess.__init__(self)
        HDF5Process.__init__(self, model_part)

    def AddPrintOutput(self, aggregated_controlled_operation: AggregatedControlledOperations) -> None:
        self._AddAggregatedControlledOperations("PrintOutput", aggregated_controlled_operation)

    def IsOutputStep(self) -> bool:
        return any([aggregated_output_operation.Evaluate() for aggregated_output_operation in self._GetAggregatedControlledOperations("PrintOutput")])

    def PrintOutput(self) -> None:
        self._ExecuteAggregatedControlledOperations("PrintOutput")
        self._UpdateAggregatedControlledOperations("PrintOutput")

##!@}
##!@}
