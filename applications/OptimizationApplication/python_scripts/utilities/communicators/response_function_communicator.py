import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.communicators.optimization_component_communicator import OptimizationComponentCommunicator
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseFunctionCommunicator:
    """Communicator class to unify communication between responses and OptimizationInfo

    Instances of this class unifies communication between OptimizationInfo
    and responses. Otherwise, a response functioncan create the string "key" with arbitary names
    which is difficult for tracking. This class creates all the necessary basic
    data structure within proper places of OptimizationInfo and gives access to them
    via defined methods such as GetBufferedDataContainer and GetUnbufferedDataContainer

    """
    def __init__(self, name: str, optimization_info: OptimizationInfo, required_buffer_size: int = 2):
        """Creates an instance of ResponseFunctionCommunicator.

        This creates an instance of response function communicator. If it finds an already existing
        response function problem data, then this communicator links to it. Otherwise, it creates
        respective response function problem data containers in OptimizationInfo

        In the event, this creates the response function problem data, then a buffered data container with
        required_buffer_size is created (given by GetBufferedDataContainer). Another unbuffered data container is
        created for sensitivitiy data given by GetUnbufferedDataContainer

        The response function should be present in the optimization_info.

        Args:
            name (str): Name of the response function.
            optimization_info (OptimizationInfo): Root or child instance of type OptimizationInfo
            required_buffer_size (int, optional): Required buffer size for this response. Defaults to 2.

        Raises:
            RuntimeError: If the given name is not found in responses.
            RuntimeError: If the required buffer size does not match the existing buffer size.
        """
        # minimum buffer size of 2 is required for the communicator
        # algorithm can request higher buffer sizes if required.
        required_buffer_size = max(required_buffer_size, 2)

        self.__optimization_component_communicator = OptimizationComponentCommunicator(optimization_info)
        problem_data = self.__optimization_component_communicator.GetProblemData()

        if not problem_data.HasValue(f"responses/{name}"):
            # first create the buffered data containers to store light objects
            # such as primitive variables
            self.__response_buffered_data = OptimizationInfo(required_buffer_size)
            problem_data[f"responses/{name}"] = self.__response_buffered_data

            # now create the sensitivities data container without buffer because
            # they store heavy objects containing sensitivities, and buffer
            # is not required.
            self.__response_unbuffered_data = OptimizationInfo(1)
            problem_data[f"responses/{name}/sensitivities"] = self.__response_unbuffered_data
        else:
            self.__response_buffered_data: OptimizationInfo = problem_data[f"responses/{name}"]
            if self.__response_buffered_data.GetBufferSize() < required_buffer_size:
                raise RuntimeError(f"The required buffer size is not satisfied with the existing problem data container. [ response data container buffer size = {self.__response_buffered_data.GetBufferSize()}, required buffer size = {required_buffer_size}")

            self.__response_unbuffered_data: OptimizationInfo = problem_data[f"responses/{name}/sensitivities"]

        self.__response_function = self.__optimization_component_communicator.GetResponseFunction(name)
        self.__name = name

    def GetName(self) -> str:
        """Get the name of the response function

        Returns:
            str: Response function name.
        """
        return self.__name

    def GetResponseFunction(self) -> ResponseFunction:
        """Get the response function

        Returns:
            ResponseFunction: Response function.
        """
        return self.__response_function

    def GetBufferedDataContainer(self) -> OptimizationInfo:
        """Get the buffered data container for response function.

        Returns:
            OptimizationInfo: Buffered data container for response function
        """
        return self.__response_buffered_data

    def GetUnbufferedDataContainer(self) -> OptimizationInfo:
        """Get the unbuffered data container for the response function.

        This method returns the unbuffered data container for response functions which
        is used to store heavy objects such as response senstivity data.

        Returns:
            OptimizationInfo: Unbuffered data container for response function.
        """
        return self.__response_unbuffered_data

    def GetScaledValue(self, step_index: int = 0, scaling: float = 1.0) -> float:
        """Gets or calculates scaled response value.

        If step__index == 0, then this method first checks whether the response value is
        alread calculated, if not calculates and gets the scaled response function value

        Otherwise, it only returns the response function value by scaling it.

        Args:
            step_index (int, optional): Step index to look the response value. Defaults to 0.
            scaling (float, optional): Scaling to be applied to response value. Defaults to 1.0.

        Returns:
            float: Scaled response value for specified step_index.
        """
        response_problem_data = self.GetBufferedDataContainer()
        if step_index == 0 and not response_problem_data.HasValue("value"):
            response_problem_data["value"] = self.__response_function.CalculateValue()

        return response_problem_data["value"] * scaling

    def GetRelativeChange(self) -> float:
        """Get the relative change ratio of the response value.

        This method returns the relative change ratio of the response function.
        If it is the first step, then this returns 0.0
        If the previous step response value is almost zero, then returns the current response value.

        Returns:
            float: Relative change ratio of the response value.
        """
        if self.__optimization_component_communicator.GetStep() > 1:
            return self.GetScaledValue() / self.GetScaledValue(1) - 1.0 if abs(self.GetScaledValue(1)) > 1e-12 else self.GetScaledValue()
        else:
            return 0.0

    def GetAbsoluteChange(self, reference_value: float) -> float:
        """Get the absolute change ratio of the response value w.r.t. reference_value.

        This method returns the absolute change ratio of the response value w.r.t. given reference value.
        If the reference value is almost zero, then this returns the response value.
        Otherwise it returns absolute change ration.

        Args:
            reference_value (float): Reference value of the response to compute the absolute change ratio.

        Returns:
            float: Absolute change ratio of the response value w.r.t. given reference value.
        """
        return self.GetScaledValue() / reference_value - 1.0 if abs(reference_value) > 1e-12 else self.GetScaledValue()


    def CalculateScaledSensitivity(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]', scaling: float = 1.0) -> None:
        """Calculates the sensitivities requested via sensitivity_variable_collective_expression_info.

        This method calculates sensitivities on the model parts w.r.t. sensitivity variables given in sensitivity_variable_collective_expression_info.
        It first checks whether for given sensitivity and model parts, the sensitivities are already
        computed, if not they are computed and the scaled response sensitivities are put inside the sensitivity_variable_collective_expression_info
        CollectiveExpressions.

        Args:
            sensitivity_variable_collective_expression_info (dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]): Sensitivity variable and CollectiveExpressions pairs
            scaling (float, optional): Scaling to be used on the computed response sensitivities. Defaults to 1.0.
        """
        response_sensitivity_data = self.GetUnbufferedDataContainer()

        # check whether the same sensitivity is computed for each sensitivity model part w.r.t. same sensitivity variable for the same response function.
        required_sensitivity_model_part_variable_info: 'dict[SupportedSensitivityFieldVariableTypes, list[Kratos.ModelPart]]' = {}
        for sensitivity_variable, collective_expression in sensitivity_variable_collective_expression_info.items():
            is_already_computed = True
            for container_expression in collective_expression.GetContainerExpressions():
                key = f"{container_expression.GetModelPart().FullName()}/raw_{self.GetName()}_{sensitivity_variable.Name()}"
                is_already_computed = is_already_computed and response_sensitivity_data.HasValue(key)

            if not is_already_computed:
                # translates container expressions to model parts which is known by responses.
                required_sensitivity_model_part_variable_info[sensitivity_variable] = [container_expression.GetModelPart() for container_expression in collective_expression.GetContainerExpressions()]

        self.__response_function.CalculateSensitivity(required_sensitivity_model_part_variable_info)

        # now store the computed sensitivities from resepective containers to container expressions
        for sensitivity_variable in required_sensitivity_model_part_variable_info.keys():
            for container_expression in sensitivity_variable_collective_expression_info[sensitivity_variable].GetContainerExpressions():
                key = f"{container_expression.GetModelPart().FullName()}/raw_{self.GetName()}_{sensitivity_variable.Name()}"
                container_expression.Read(sensitivity_variable)
                # this cloning is cloning of an expression pointer, hence it is not expensive.
                # no cloning of underlying sensitivity data is involved.
                response_sensitivity_data[key] = container_expression.Clone()

        # now assign the sensitivities for all the collective expressions correctly
        for sensitivity_variable, collective_expression in sensitivity_variable_collective_expression_info.items():
            for container_expression in collective_expression.GetContainerExpressions():
                key = f"{container_expression.GetModelPart().FullName()}/raw_{self.GetName()}_{sensitivity_variable.Name()}"
                # this copying is copying of an expression pointer, hence it is not expensive.
                # no copying of underlying sensitivity data is involved.
                container_expression.CopyFrom(response_sensitivity_data[key])
                container_expression *= scaling
