from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo

class OptimizationComponentCommunicator:
    """Communicator class to unify communication between components and OptimizationInfo

    Instances of this class unifies communication between OptimizationInfo
    and optimization components. Otherwise, a component of optimization info
    can create the string "key" with arbitary names which is difficult for tracking.

    """
    def __init__(self, optimization_info: OptimizationInfo) -> None:
        """Creates an instance of the OptimizationComponentCommunicator

        This creates an instance of OptimizationComponentCommunicator, which
        is connected to existing components. If no components are found, then
        an entry is created to store components.

        provided optimization_info can be an instance of OptimizationInfo
        which is either the root or a child at any level.

        Args:
            optimization_info (OptimizationInfo): OptimizationInfo root or child at any level.
        """
        root_optimization_info = optimization_info.GetRoot()
        if not root_optimization_info.HasValue("components"):
            # creates a components sub item with buffer size 1.
            # The buffer size is 1 because, this will only be used
            # to hold components required for optimization.
            root_optimization_info["components"] = OptimizationInfo(1)

        if not root_optimization_info.HasValue("problem_data"):
            # creates a components sub item with buffer size 1.
            # The buffer size is 1 because, it is the responsibility of other
            # communicators to create their own sub items with required
            # buffer sizes.
            root_optimization_info["problem_data"] = OptimizationInfo(1)
            # now set the step variable
            root_optimization_info["problem_data/step"] = 0

        self.__components: OptimizationInfo = root_optimization_info["components"]
        self.__problem_data: OptimizationInfo = root_optimization_info["problem_data"]

    def GetProblemData(self) -> OptimizationInfo:
        """Get the problem data from the OptimizationInfo

        Returns:
            OptimizationInfo: Problem data optimization info
        """
        return self.__problem_data

    def AdvanceStep(self) -> None:
        """Advances problem data to the next step.

        This method advances the problem data to next step index and
        clears all the data on that step index recursively.

        This also increases the step value by one.

        """
        self.__problem_data.AdvanceStep()
        self.__problem_data.ClearStep()
        self.__problem_data["step"] = self.__problem_data["step", 1] + 1

    def GetStep(self) -> int:
        """Get the current step

        Returns:
            int: Current step of the optimization info
        """
        return self.__problem_data["step"]

    def AddResponseFunction(self, name: str, response_function: ResponseFunction) -> None:
        """Adds response function under the name

        This method adds the provided response function under the name. The name must be
        unique.

        Args:
            name (str): Name of the response function.
            response_function (ResponseFunction): Response function to assign to the name.

        Raises:
            RuntimeError: If an existing response function is found with the given name.
        """
        self.__components[f"responses/{name}"] = response_function

    def GetResponseFunction(self, name: str) -> ResponseFunction:
        """Gets the response function for given name.

        This method returns a response function matching the given name.

        Args:
            name (str): Name of the response function.

        Raises:
            RuntimeError: If a reesponse function is not found with the given name.

        Returns:
            ResponseFunction: Response function matching the given name.
        """
        return self.__components[f"responses/{name}"]

    def AddExecutionPolicyDecorator(self, execution_policy_decorator: ExecutionPolicyDecorator) -> None:
        """Adds execution policy decorator under the name

        This method adds the provided execution policy decorator under the name. The name must be
        unique.

        Args:
            name (str): Name of the execution policy decorator.
            execution_policy_decorator (ExecutionPolicyDecorator): execution policy decorator to assign to the name.

        Raises:
            RuntimeError: If an existing execution policy decorator is found with the given name.
        """
        self.__components[f"execution_policy_decorators/{execution_policy_decorator.GetExecutionPolicyName()}"] = execution_policy_decorator

    def GetExecutionPolicyDecorator(self, name: str) -> ExecutionPolicyDecorator:
        """Gets the execution policy decorator for given name.

        This method returns a execution policy decorator matching the given name.

        Args:
            name (str): Name of the execution policy decorator.

        Raises:
            RuntimeError: If a execution policy decorator is not found with the given name.

        Returns:
            ResponseFunction: execution policy decorator matching the given name.
        """
        return self.__components[f"execution_policy_decorators/{name}"]