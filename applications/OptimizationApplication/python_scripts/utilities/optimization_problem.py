from __future__ import annotations
from typing import Any

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
    def __init__(self) -> None:
        """Creates an instance of Optimization info

        Creates an instance of optimization info with most basic structure
        for the BufferedDict container and components container.

        """
        # create the optimization components data container with basic types
        self.__components: 'dict[Any, dict[str, Any]]' = {
            ResponseFunction: {},
            ExecutionPolicy: {},
            Control: {}
        }

        # create the unbufferd optimization problem data container
        self.__problem_data = BufferedDict(1)

        # initialize the step
        self.__problem_data["step"] = 0

    def GetComponentType(self, component: Any) -> Any:
        for k in self.__components.keys():
            if isinstance(component, k):
                return k

        return None

    def GetComponentName(self, component: Any) -> str:
        component_type = self.GetComponentType(component)
        if not component_type is None:
            components = self.__components[component_type]
        else:
            raise RuntimeError(f"Unsupported component type {component} provided.")

        for comp_name, current_component in components.items():
            if component == current_component:
                return comp_name

        raise RuntimeError(f"The given {component} not found in components of type {component_type}.")

    def AddResponse(self, name: str, response_function: ResponseFunction) -> None:
        self.__AddComponent(name, response_function, ResponseFunction)

    def GetResponse(self, name: str) -> ResponseFunction:
        return self.__GetComponent(name, ResponseFunction)

    def GetReponseData(self, name: str) -> BufferedDict:
        return self.__GetComponentDataContainer(name, ResponseFunction)

    def RemoveResponse(self, name: str) -> None:
        self.__RemoveComponent(name, ResponseFunction)
        self.__RemoveComponentDataContainer(name, ResponseFunction)

    def AddExecutionPolicy(self, name: str, execution_policy: ExecutionPolicy) -> None:
        self.__AddComponent(name, execution_policy, ExecutionPolicy)

    def GetExecutionPolicy(self, name: str) -> ExecutionPolicy:
        return self.__GetComponent(name, ExecutionPolicy)

    def GetExecutionPolicyData(self, name: str) -> BufferedDict:
        return self.__GetComponentDataContainer(name, ExecutionPolicy)

    def RemoveExecutionPolicy(self, name: str) -> None:
        self.__RemoveComponent(name, ExecutionPolicy)
        self.__RemoveComponentDataContainer(name, ExecutionPolicy)

    def AddControl(self, name: str, control: Control) -> None:
        self.__AddComponent(name, control, Control)

    def GetControl(self, name: str) -> Control:
        return self.__GetComponent(name, Control)

    def RemoveControl(self, name: str) -> None:
        self.__RemoveComponent(name, Control)
        self.__RemoveComponentDataContainer(name, Control)

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

    def __AddComponent(self, name: str, component: Any, component_type: Any) -> None:
        if not isinstance(component, component_type):
            raise RuntimeError(f"Trying to add \"{name}\" of type {component.__class__.__name__} which is not derrived from {component_type.__name__}.")

        components = self.__components[component_type]
        if name in components.keys():
            raise RuntimeError(f"Trying to add \"{name}\" of type {component_type.__name__} when there exists an item with the same name.")

        components[name] = component

    def __GetComponent(self, name: str, component_type: Any) -> Any:
        components = self.__components[component_type]
        if name in components.keys():
            return components[name]
        else:
            raise RuntimeError(f"No {component_type.__name__} with name = \"{name}\" exists. Followings are the available {component_type.__name__}:\n" + "\n\t".join(components.keys()))

    def __RemoveComponent(self, name: str, component_type: Any) -> None:
        components = self.__components[component_type]
        if name in components.keys():
            del components[name]
        else:
            raise RuntimeError(f"No {component_type.__name__} with name = \"{name}\" exists. Followings are the available {component_type.__name__}:\n" + "\n\t".join(components.keys()))

    def __GetComponentDataContainer(self, name: str, component_type: Any) -> BufferedDict:
        components = self.__components[component_type]

        if not name in components.keys():
            raise RuntimeError(f"No {component_type.__name__} with name = \"{name}\" exists. Followings are the available {component_type.__name__}:\n" + "\n\t".join(components.keys()))

        data_name = f"{component_type.__name__}/{name}"
        if not self.__problem_data.HasValue(data_name):
            self.__problem_data[data_name] = {}

        return self.__problem_data[data_name]

    def __RemoveComponentDataContainer(self, name: str, component_type: Any) -> None:
        data_name = f"{component_type.__name__}/{name}"
        if self.__problem_data.HasValue(data_name):
            del self.__problem_data[data_name]

