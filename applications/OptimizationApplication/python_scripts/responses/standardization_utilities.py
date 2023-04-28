import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class StandardizationUtilities:
    """Communicator class to simplify standardization between responses and OptimizationProblem

    Instances of this class simplifies interaction between OptimizationProblem
    and responses.

    """
    def __init__(self, name: str, optimization_info: OptimizationProblem, required_buffer_size: int = 2):
        """Creates an instance of StandardizationUtilities.

        This creates an instance of response function utility. If it finds an already existing
        response function problem data, then this communicator links to it. Otherwise, it creates
        respective response function problem data containers in OptimizationProblem

        In the event, this creates the response function problem data, then a buffered data container with
        required_buffer_size is created (given by GetBufferedData). Another unbuffered data container is
        created for sensitivitiy data given by GetUnbufferedData

        The response function should be present in the optimization_info.

        Args:
            name (str): Name of the response function.
            optimization_info (OptimizationProblem): Root or child instance of type OptimizationProblem
            required_buffer_size (int, optional): Required buffer size for this response. Defaults to 2.

        Raises:
            RuntimeError: If the given name is not found in responses.
            RuntimeError: If the required buffer size does not match the existing buffer size.
        """
        self.__optimization_info = optimization_info

        # minimum buffer size of 2 is required for the communicator
        # algorithm can request higher buffer sizes if required.
        required_buffer_size = max(required_buffer_size, 2)
        response_data = optimization_info.GetReponseData(name)
        self.__response_function: ResponseFunction = optimization_info.GetResponse(name)
        self.__name = name

        if not response_data.HasValue("buffered"):
            # first create the buffered data containers to store light objects
            # such as primitive variables
            self.__response_buffered_data = BufferedDict(required_buffer_size)
            response_data["buffered"] = self.__response_buffered_data

            # now create the sensitivities data container without buffer because
            # they store heavy objects containing sensitivities, and buffer
            # is not required.
            self.__response_unbuffered_data = BufferedDict(1)
            response_data["unbuffered"] = self.__response_unbuffered_data
        else:
            self.__response_buffered_data: BufferedDict = response_data["buffered"]
            if self.__response_buffered_data.GetBufferSize() < required_buffer_size:
                raise RuntimeError(f"The required buffer size is not satisfied with the existing problem data container. [ response data container buffer size = {self.__response_buffered_data.GetBufferSize()}, required buffer size = {required_buffer_size}")

            self.__response_unbuffered_data: BufferedDict = response_data["unbuffered"]

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

    def GetBufferedData(self) -> OptimizationProblem:
        """Get the buffered data container for response function.

        Returns:
            OptimizationProblem: Buffered data container for response function
        """
        return self.__response_buffered_data

    def GetUnbufferedData(self) -> OptimizationProblem:
        """Get the unbuffered data container for the response function.

        This method returns the unbuffered data container for response functions which
        is used to store heavy objects such as response senstivity data.

        Returns:
            OptimizationProblem: Unbuffered data container for response function.
        """
        return self.__response_unbuffered_data

    def GetScaledValue(self, step_index: int = 0, scaling: float = 1.0) -> float:
        """Gets or calculates scaled response value.

        If step__index == 0, then this method first checks whether the response value is
        already calculated, if not calculates and gets the scaled response function value

        Otherwise, it only returns the response function value by scaling it.

        Args:
            step_index (int, optional): Step index to look the response value. Defaults to 0.
            scaling (float, optional): Scaling to be applied to response value. Defaults to 1.0.

        Returns:
            float: Scaled response value for specified step_index.
        """
        response_problem_data = self.GetBufferedData()
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
        if self.__optimization_info.GetStep() > 1:
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
        response_sensitivity_data = self.GetUnbufferedData()

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
