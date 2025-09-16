import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAlgorithmTimeLogger
from KratosMultiphysics.SystemIdentificationApplication.responses.damage_detection_response import DamageDetectionResponse
from KratosMultiphysics.SystemIdentificationApplication.controls.material_properties_control import MaterialPropertiesControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_steepest_descent import AlgorithmSteepestDescent


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmMaskBasedSteepestDescent(model, parameters, optimization_problem)

class AlgorithmMaskBasedSteepestDescent(AlgorithmSteepestDescent):
    """
        A classical steepest descent algorithm to solve unconstrained optimization problems.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "objective"         : {},
            "controls"          : [],
            "echo_level"        : 0,
            "settings"          : {
                "echo_level"      : 0,
                "mask_name"       : "",
                "line_search"     : {},
                "conv_settings"   : {}
            }
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(model, parameters, optimization_problem)
        self.mask_name = self.parameters["settings"]["mask_name"].GetString()


    @time_decorator()
    def Solve(self):
        while not self.converged:
            with OptimizationAlgorithmTimeLogger("AlgorithmSteepestDescent",self._optimization_problem.GetStep()):
                self._InitializeIteration()

                response: DamageDetectionResponse = self._objective.GetReponse()
                control: MaterialPropertiesControl = self.master_control.GetListOfControls()[0]

                control.skip_forward_filter = True

                # clear all the sensors
                response.damage_response_function.Clear()
                for sensor in response.list_of_sensors:
                    # now we raw sensitivities for each sensor
                    response.damage_response_function.AddSensor(sensor)

                self._obj_val = self._objective.CalculateStandardizedValue(self._control_field)
                obj_info = self._objective.GetInfo()
                self.algorithm_data.GetBufferedData()["std_obj_value"] = obj_info["std_value"]
                self.algorithm_data.GetBufferedData()["rel_obj[%]"] = obj_info["rel_change [%]"]
                if "abs_change [%]" in obj_info:
                    self.algorithm_data.GetBufferedData()["abs_obj[%]"] = obj_info["abs_change [%]"]

                obj_grad = self._control_field.Clone() * 0.0
                for sensor in response.list_of_sensors:
                    # mask = Kratos.Expression.ElementExpression(self._control_field.GetContainerExpressions()[0].GetModelPart())
                    # Kratos.Expression.LiteralExpressionIO.SetData(mask, 1.0)
                    mask = sensor.GetContainerExpression(self.mask_name)

                    # clear all the sensors
                    response.damage_response_function.Clear()
                    # now we raw sensitivities for each sensor
                    response.damage_response_function.AddSensor(sensor)
                    raw_sensor_grad = self._control_field.Clone()

                    control.GetEmptyField().GetModelPart().ProcessInfo[KratosSI.SENSOR_NAME] = sensor.GetName()

                    response.CalculateGradient({Kratos.YOUNG_MODULUS: raw_sensor_grad})
                    mapped_sensor_grad = control.MapGradient({Kratos.YOUNG_MODULUS: raw_sensor_grad.GetContainerExpressions()[0]})
                    modified_sensor_grad = control.filter.ForwardFilterField(mapped_sensor_grad) * mask
                    obj_grad += KratosOA.CollectiveExpression([modified_sensor_grad])

                    self.algorithm_data.GetBufferedData()[f"{sensor.GetName()}_mask"] = mask.Clone()
                    self.algorithm_data.GetBufferedData()[f"{sensor.GetName()}_raw_gradient"] = raw_sensor_grad.Clone()
                    self.algorithm_data.GetBufferedData()[f"{sensor.GetName()}_mapped_gradient"] = mapped_sensor_grad.Clone()
                    self.algorithm_data.GetBufferedData()[f"{sensor.GetName()}_modified_sensor_gradient"] = modified_sensor_grad.Clone()

                self.ComputeSearchDirection(obj_grad)

                alpha = self._line_search_method.ComputeStep()

                self.ComputeControlUpdate(alpha)

                self._FinalizeIteration()

                self.Output()

                self.UpdateControl()

                self.converged = self._convergence_criteria.IsConverged()

                self._optimization_problem.AdvanceStep()

        return self.converged