import os
import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.convergence_criteria.convergence_criterion import ConvergenceCriterion
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentHavingDataByFullName
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ConvergenceCriterion:
    if not parameters.Has("settings"):
        raise RuntimeError(f"MaxSensorErrorCriterion instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return MaxSensorErrorCriterion(parameters["settings"], optimization_problem)

class MaxSensorErrorCriterion(ConvergenceCriterion):
    """
    MaxSensorErrorCriterion is a convergence criterion based on the maximum sensor error in a system identification problem.

    This class checks whether the maximum sensor error falls below a specified tolerance,
    indicating convergence. The criterion is configurable via parameters.
    """
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "component_name": "algorithm",
            "sensor_group_name": "sensors",
            "tolerance"     : 1e-9,
            "output_to_file" : true,
            "sensor_error_output_file_name" : "max_sensor_error.csv"
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.__component_name = parameters["component_name"].GetString()
        self.__sensor_group_name = parameters["sensor_group_name"].GetString()
        self.__tolerance = parameters["tolerance"].GetDouble()
        self.__optimization_problem = optimization_problem
        self.__output_to_file = parameters["output_to_file"].GetBool()
        self.__sensor_error_output_file_name = parameters["sensor_error_output_file_name"].GetString()

    def Initialize(self):
        component = GetComponentHavingDataByFullName(self.__component_name, self.__optimization_problem)
        self.__component_data_view = ComponentDataView(component, self.__optimization_problem)
        if self.__output_to_file:
            with open(self.__sensor_error_output_file_name, "w") as file_output:
                file_output.write("iteration,max_sensor_error\n")

    @time_decorator()
    def IsConverged(self) -> bool:
        """
        Check if the convergence criterion based on sensor error is satisfied.
        This method evaluates convergence by comparing the maximum sensor error across
        all sensors in the sensor group against a predefined tolerance threshold.
        The method performs the following steps:
        1. Retrieves sensor group data from the optimization problem
        2. Extracts all sensors from the sensor group
        3. Computes the maximum absolute error across all sensors
        4. Compares the maximum error against the tolerance threshold
        5. Writes the maximum sensor error to a file for monitoring
        Returns:
            bool: True if the maximum sensor error is less than or equal to the 
                  tolerance threshold, False otherwise.
        """
        sensor_group_data = ComponentDataView(self.__sensor_group_name, self.__optimization_problem)
        self.list_of_sensors = GetSensors(sensor_group_data)
        
        max_sensor_error = 0.0
        for sensor in self.list_of_sensors:
            max_sensor_error = max(max_sensor_error, abs(sensor.GetNode().GetValue(KratosSI.SENSOR_ERROR)))
        
        self.__max_sensor_error = max_sensor_error
        self.__conv = self.__max_sensor_error <= self.__tolerance
        self.WriteMaxSensorErrorToFile()
        return self.__conv

    def GetInfo(self) -> 'list[tuple[str, typing.Union[int, float, str]]]':
        info = [
                    ('type'          , 'sensor_error'),
                    ('max_sensor_error', self.__max_sensor_error),
                    ('tolerance'     , self.__tolerance),
                    ('component_name', self.__component_name),
                    ('sensor_group_name', self.__sensor_group_name),
                    ('status'        , str("converged" if self.__conv else "not converged"))
               ]
        return info

    def Finalize(self):
        return super().Finalize()

    def WriteMaxSensorErrorToFile(self) -> None:
        if self.__output_to_file:
            iteration = self.__optimization_problem.GetStep()
            file_exists = os.path.isfile(self.__sensor_error_output_file_name)
            with open(self.__sensor_error_output_file_name, "a") as file_output:
                if not file_exists:
                    file_output.write("iteration,max_sensor_error\n")
                file_output.write(f"{iteration},{self.__max_sensor_error}\n")
