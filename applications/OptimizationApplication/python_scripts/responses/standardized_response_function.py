import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class StandardizedResponseFunction:
    """Communicator class to simplify standardization between responses and OptimizationProblem

    Instances of this class simplifies interaction between OptimizationProblem
    and responses.

    """
    def __init__(self, master_control: MasterControl, response_name: str, optimization_info: OptimizationProblem, required_buffer_size: int = 2):
        """Creates an instance of StandardizedResponseFunction.

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
        response_data = optimization_info.GetReponseData(response_name)

        # create the response routine.
        self.__response_routine = ResponseRoutine(master_control, response_name, optimization_info)

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
        return self.__response_routine.GetReponseName()

    def GetResponseRoutine(self) -> ResponseRoutine:
        """Get the response function

        Returns:
            ResponseRoutine: Response function.
        """
        return self.__response_routine

    def GetBufferedData(self) -> BufferedDict:
        """Get the buffered data container for response function.

        Returns:
            BufferedDict: Buffered data container for response function
        """
        return self.__response_buffered_data

    def GetUnbufferedData(self) -> BufferedDict:
        """Get the unbuffered data container for the response function.

        This method returns the unbuffered data container for response functions which
        is used to store heavy objects such as response senstivity data.

        Returns:
            BufferedDict: Unbuffered data container for response function.
        """
        return self.__response_unbuffered_data
    
    def CalculateScaledValue(self, control_space_updates: KratosOA.ContainerExpression.CollectiveExpressions, scaling: float = 1.0) -> float:
        """Calculates the scaled response value from the routine.

        Args:
            control_space_updates (KratosOA.ContainerExpression.CollectiveExpressions): Control space final design.
            scaling (float, optional): Scaling to be used. Defaults to 1.0.

        Returns:
            float: Calculated scaled response value.
        """
        return self.__response_routine.CalculateValue(control_space_updates) * scaling
    
    def GetReponseValue(self, step_index: int = 0) -> float:
        return self.GetBufferedData().GetValue("value", step_index)

    def SetResponseValue(self, value: float, step_index: int = 0) -> None:
        self.GetBufferedData().SetValue("value", value, step_index)

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

    def CalculateScaledGradient(self, scaling: float = 1.0) -> KratosOA.ContainerExpression.CollectiveExpressions:
        """Returns the scaled control space collective expressions of gradients.

        Args:
            scaling (float, optional): Scaling. Defaults to 1.0.

        Returns:
            KratosOA.ContainerExpression.CollectiveExpressions: Scaled control space gradients.
        """
        return self.__response_routine.CalculateGradient() * scaling
