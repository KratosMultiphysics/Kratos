from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ContainerVariableDataHolderUnion
from KratosMultiphysics.OptimizationApplication.output.output_utilities import GetOptimizationInfoAvailableKeysForType
from KratosMultiphysics.OptimizationApplication.output.output_utilities import GetPatternMatchedOptimizationKeys
from KratosMultiphysics.OptimizationApplication.output.output_utilities import GetPrefixWithoutDataName
from KratosMultiphysics.OptimizationApplication.output.output_utilities import GetPlaceHolderWithValues

class ContainerVariableDataHolderInfo(ABC):
    @classmethod
    @abstractmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        pass

    @classmethod
    @abstractmethod
    def GetPattern(cls) -> str:
        pass

    def __init__(self, optimization_info: OptimizationInfo, settings: Kratos.Parameters):
        default_parameters = self.GetDefaultParameters()
        settings.ValidateAndAssignDefaults(default_parameters)

        self.optimization_info = optimization_info
        self.common_prefix = GetPrefixWithoutDataName(self.GetPattern(), settings)
        self.place_holder_values = GetPlaceHolderWithValues(self.GetPattern(), settings)

        self.data_list: 'list[dict[str, any]]' = []
        for data in settings["data_list"]:
            data.ValidateAndAssignDefaults(default_parameters["data_list"][0])
            self.data_list.append({
                "location" : f'{self.common_prefix}/{data["data_name"].GetString()}',
                "variable" : Kratos.KratosGlobals.GetVariable(data["variable_name"].GetString()),
                "data_name": data["data_name"].GetString(),
                "container_type": settings["container_type"].GetString()
            })

    def Check(self):
        available_keys = GetOptimizationInfoAvailableKeysForType(self.optimization_info, ContainerVariableDataHolderUnion)

        # check for model part and response name availability
        if not self.optimization_info.HasValue(self.common_prefix):
            available_options = GetPatternMatchedOptimizationKeys(available_keys, self.GetPattern(), False)
            raise RuntimeError(f"The provided {self.place_holder_values} pair not found. Followings are available pairs:\n\t" + "\n\t".join(available_options))

        # now check for available data names
        for data in self.data_list:
            data_name = data["data_name"]
            location = data["location"]
            variable = data["variable"]
            container_type = data["container_type"]

            if not self.optimization_info.HasValue(location):
                available_options = GetPatternMatchedOptimizationKeys(available_keys, self.common_prefix, True)
                raise RuntimeError(f"The provided \"data_name\" = \"{data_name}\" for  {self.place_holder_values} not found. Followings are available data names:\n\t" + "\n\t".join(available_options))
            else:
                container_data: ContainerVariableDataHolderUnion = self.optimization_info.GetValue(location)
                variable_data_dimension = KratosOA.OptimizationUtils.GetVariableDimension(variable, container_data.GetModelPart().ProcessInfo[Kratos.DOMAIN_SIZE])
                if container_data.GetDataDimension() != variable_data_dimension:
                    raise RuntimeError(f"The requested \"data_name\" = \"{data_name}\" with {self.place_holder_values} is having data with dimensions {container_data.GetDataDimension()} which is not matching with data dimension {variable_data_dimension} of {variable.Name()}. Followings are the details of the data container:\n\t{container_data}")

                current_container_type = ContainerVariableDataHolderInfo.GetContainerType(container_data)
                if container_type != current_container_type:
                    raise RuntimeError(f"The requested \"data_name\" = {data_name} with {self.place_holder_values} is having a different data container than requested. Followings are the details:\n\tRequested data type: {container_type}\n\tFound data type: {current_container_type}")

    def ApplyDataConainer(self):
        self.Check()
        for data in self.data_list:
            location: str = data["location"]
            variable: any = data["variable"]
            container_data: ContainerVariableDataHolderUnion = self.optimization_info.GetValue(location)

            # now do the conversion from historical to non historical
            converted_data = container_data
            if isinstance(converted_data, KratosOA.HistoricalContainerVariableDataHolder):
                converted_data = KratosOA.NodalContainerVariableDataHolder(container_data)
            elif isinstance(converted_data, KratosOA.ConditionPropertiesContainerVariableDataHolder):
                converted_data = KratosOA.ConditionContainerVariableDataHolder(container_data)
            elif isinstance(converted_data, KratosOA.ElementPropertiesContainerVariableDataHolder):
                converted_data = KratosOA.ElementContainerVariableDataHolder(container_data)

            # store the existing data to temp
            temp = type(converted_data)(converted_data.GetModelPart())
            temp.ReadDataFromContainerVariable(variable)

            # assign data container
            converted_data.AssignDataToContainerVariable(variable)

            # store temp data in the dict
            data["temp"] = temp

    def ResetData(self):
        for data in self.data_list:
            data["temp"].AssignDataToContainerVariable(data["variable"])
            del data["temp"]

    @staticmethod
    def GetContainerType(container: ContainerVariableDataHolderUnion) -> str:
        if isinstance(container, KratosOA.NodalContainerVariableDataHolderBase):
            return "nodes"
        elif isinstance(container, KratosOA.ConditionContainerVariableDataHolderBase):
            return "condition"
        elif isinstance(container, KratosOA.ElementContainerVariableDataHolderBase):
            return "elements"
        else:
            raise RuntimeError(f"The container type cannot be determined. [ Container = {str(container)}]")

class ResponseContainerVariableDataHolderInfo(ContainerVariableDataHolderInfo):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"           : "response_data",
            "model_part_name": "PLEASE_PROVIDE_MODEL_PART_NAME",
            "container_type" : "PLEASE_PROVIDE_CONTAINER_TYPE",
            "response_name"  : "PLEASE_PROVIDE_RESPONSE_FUNCTION_NAME",
            "data_list"      : [
                {
                    "data_name"     : "PLEASE_PROVIDE_DATA_NAME",
                    "variable_name" : "PLEASE_PROVIDE_VARIABLE_NAME"
                }
            ]
        }""")

    @classmethod
    def GetPattern(cls) -> str:
        return "problem_data/response_data/<model_part_name>/<response_name>/sensitivities"

class ControlContainerVariableDataHolderInfo(ContainerVariableDataHolderInfo):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"           : "control_data",
            "model_part_name": "PLEASE_PROVIDE_MODEL_PART_NAME",
            "container_type" : "PLEASE_PROVIDE_CONTAINER_TYPE",
            "control_name"   : "PLEASE_PROVIDE_CONTROL_NAME",
            "data_list"      : [
                {
                    "data_name"     : "PLEASE_PROVIDE_DATA_NAME",
                    "variable_name" : "PLEASE_PROVIDE_VARIABLE_NAME"
                }
            ]
        }""")

    @classmethod
    def GetPattern(cls) -> str:
        return "problem_data/control_data/<model_part_name>/<control_name>"

class ContainerVariableDataHolderOutputDecorator(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optmization_info: OptimizationInfo):
        super().__init__()

        default_parameters = Kratos.Parameters("""{
            "module"                           : "KratosMultiphysics.OptimizationApplication.utilities",
            "type"                             : "ContainerVariableDataHolderOutputDecorator",
            "output_processes"                 : [],
            "optimization_info_output_settings": []
        }""")
#
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__list_of_data_sets = []
        optimization_info_output_settings: Kratos.Parameters
        for optimization_info_output_settings in parameters["optimization_info_output_settings"]:
            if not optimization_info_output_settings.Has("type"):
                raise RuntimeError(f"\"type\" is not defined for optimization_info_output_settings:\n" + str(optimization_info_output_settings))

            output_type = optimization_info_output_settings["type"].GetString()
            if output_type == "response_data":
                self.__list_of_data_sets.append(ResponseContainerVariableDataHolderInfo(optmization_info, optimization_info_output_settings))
            elif output_type == "control_data":
                self.__list_of_data_sets.append(ControlContainerVariableDataHolderInfo(optmization_info, optimization_info_output_settings))
            else:
                raise RuntimeError(f"Unsupported \"{output_type}\" found. Followings are available options:\n\tresponse_data\n\tcontrol_data")

        self.__output_processes_list: 'list[Kratos.OutputProcess]' = KratosProcessFactory(model).ConstructListOfProcesses(parameters["output_processes"])
        for output_process in self.__output_processes_list:
            if not isinstance(output_process, Kratos.OutputProcess):
                raise RuntimeError(f"The provided output process is not of the type Kratos.OutputProcess. Details of the process: " + str(output_process))

    def ExecuteInitialize(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteInitialize)

    def ExecuteBeforeSolutionLoop(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteBeforeSolutionLoop)

    def ExecuteInitializeSolutionStep(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteInitializeSolutionStep)

    def ExecuteBeforeOutputStep(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteBeforeOutputStep)

    def ExecuteAfterOutputStep(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteAfterOutputStep)

    def ExecuteFinalizeSolutionStep(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteFinalizeSolutionStep)

    def ExecuteFinalize(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteFinalize)

    def IsOutputStep(self) -> bool:
        for output_process in self.__output_processes_list:
            if output_process.IsOutputStep():
                return True
        return True

    def PrintOutput(self):
        # write containers to respective model parts
        for container_variable_data_holder_info in self.__list_of_data_sets:
            container_variable_data_holder_info.ApplyDataConainer()

        # write output
        for output_process in self.__output_processes_list:
            output_process.PrintOutput()

        for container_variable_data_holder_info in self.__list_of_data_sets:
            container_variable_data_holder_info.ResetData()





