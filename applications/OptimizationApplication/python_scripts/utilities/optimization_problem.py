from __future__ import annotations
from typing import Any, Union

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control

class OptimizationProblem:
    """This is the main data holder for optimization problems

    This class holds one private @ref BufferedDict container
    which is used to hold components of the optimization problem being solved,
    and the problem data generated while solving the optimization problem.

    """
    def __init__(self, echo_level: int = 0) -> None:
        """Creates an instance of Optimization info

        Creates an instance of optimization info with most basic structure
        for the BufferedDict container and components container.

        """
        self.__echo_level = echo_level

        # create the optimization components data container with basic types
        self.__components: 'dict[Any, dict[str, Any]]' = {
            ResponseFunction: {},
            ExecutionPolicy: {},
            Control: {}
        }

        # now set the processes dict with different categories of processes
        # such as initialization, output processes
        self.__proceses: 'dict[str, list[Kratos.Process]]' = {}

        # create the unbufferd optimization problem data container
        self.__problem_data = BufferedDict(1)

        # initialize the step
        self.__problem_data["step"] = 0

    def GetComponentType(self, component: Union[ExecutionPolicy, ResponseFunction, Control]) -> Any:
        for k in self.__components.keys():
            if isinstance(component, k):
                return k

        return None

    def GetComponentName(self, component: Union[ExecutionPolicy, ResponseFunction, Control]) -> str:
        component_type = self.GetComponentType(component)
        if not component_type is None:
            components = self.__components[component_type]
        else:
            raise RuntimeError(f"Unsupported component type {component} provided.")

        for comp_name, current_component in components.items():
            if component == current_component:
                return comp_name

        raise RuntimeError(f"The given {component} not found in components of type {component_type}.")

    def AddComponent(self, component: Union[ExecutionPolicy, ResponseFunction, Control]) -> None:
        for k, v in self.__components.items():
            if isinstance(component, k):
                added_component_type = k.__name__
                v[component.GetName()] = component
                break

        if added_component_type is None:
            raise RuntimeError(f"The given compoennt is of not supported types. Supported types:\n\tKratos.Process\n\t" + "\n\t".join([k.__name__ for k in self.__components.keys()]))
        else:
            if self.__echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Added \"{component.GetName()}\" to \"{added_component_type}\".")

    def GetResponse(self, name: str) -> ResponseFunction:
        return self.GetComponent(name, ResponseFunction)

    def GetListOfResponses(self) -> 'list[ResponseFunction]':
        return self.__components[ResponseFunction].values()

    def RemoveResponse(self, name: str) -> None:
        self.RemoveComponent(name, ResponseFunction)

    def GetExecutionPolicy(self, name: str) -> ExecutionPolicy:
        return self.GetComponent(name, ExecutionPolicy)

    def GetListOfExecutionPolicies(self) -> 'list[ExecutionPolicy]':
        return self.__components[ExecutionPolicy].values()

    def RemoveExecutionPolicy(self, name: str) -> None:
        self.RemoveComponent(name, ExecutionPolicy)

    def GetControl(self, name: str) -> Control:
        return self.GetComponent(name, Control)

    def GetListOfControls(self) -> 'list[Control]':
        return self.__components[Control].values()

    def RemoveControl(self, name: str) -> None:
        self.RemoveComponent(name, Control)

    def AddProcessType(self, process_type: str) -> None:
        self.__proceses[process_type] = []

    def AddProcess(self, process_type: str, process: Kratos.Process) -> None:
        if process_type not in self.__proceses.keys():
            self.__proceses[process_type]: 'list[Kratos.Process]' = []

        if not isinstance(process, Kratos.Process):
            raise RuntimeError(f"The provided process is not of type Kratos.Process. [ process = {process}]")

        if self.__echo_level > 0:
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Added \"{process.__class__.__name__}\" to \"{process_type}\".")
        self.__proceses[process_type].append(process)

    def GetAvailableProcessTypes(self) -> 'list[str]':
        return self.__proceses.keys()

    def GetListOfProcesses(self, process_type: str) -> 'list[Kratos.Process]':
        if process_type not in self.__proceses.keys():
            raise RuntimeError(f"The process type = \"{process_type}\" not found. Followings are available process types:\n\t" + "\n\t".join([k for k in self.__proceses.keys()]))

        return self.__proceses[process_type]

    def GetStep(self) -> int:
        """Gets the current step of the optimization info

        Returns:
            int: Current step of the optimization info.
        """
        return self.__problem_data["step"]

    def AdvanceStep(self) -> None:
        """Advances the problem data by one step.

        This method advances problem data by one step and
        clears all the data in the advanced step.

        """
        current_step = self.__problem_data["step"]
        self.__problem_data.AdvanceStep()
        self.__problem_data["step"] = current_step + 1

    def GetProblemDataContainer(self) -> BufferedDict:
        """Gets the global problem data container.

        Returns:
            BufferedDict: Global problem data container.
        """
        return self.__problem_data

    def GetComponentContainer(self) -> 'dict[Any, dict[str, Any]]':
        return self.__components

    def GetComponent(self, name: str, component_type: Any) -> Any:
        if component_type not in self.__components.keys():
            raise RuntimeError(f"Unsupported component type requested from GetComponent. [ requested component type = \"{component_type}\" ]. Followings are supported component types:\n\t" + "\n\t".join([component_type.__name__ for component_type in list(self.__components.keys())]))

        components = self.__components[component_type]
        if name in components.keys():
            return components[name]
        else:
            raise RuntimeError(f"No {component_type.__name__} with name = \"{name}\" exists. Followings are the available {component_type.__name__}:\n" + "\n\t".join(components.keys()))

    def RemoveComponent(self, name: str, component_type: Any) -> None:
        if component_type not in self.__components.keys():
            raise RuntimeError(f"Unsupported component type requested from RemoveComponent. [ requested component type = \"{str(component_type)}\" ]. Followings are supported component types:\n\t" + "\n\t".join([component_type.__name__ for component_type in list(self.__components.keys())]))

        components = self.__components[component_type]
        if name in components.keys():
            del components[name]
        else:
            raise RuntimeError(f"No {component_type.__name__} with name = \"{name}\" exists. Followings are the available {component_type.__name__}:\n" + "\n\t".join(components.keys()))
