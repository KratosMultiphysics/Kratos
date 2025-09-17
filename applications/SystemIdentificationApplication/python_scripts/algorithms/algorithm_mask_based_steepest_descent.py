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
    def ComputeControlUpdate(self, alpha, list_of_gradients):
        search_direction = self.algorithm_data.GetBufferedData()["search_direction"].Clone()
        response: DamageDetectionResponse = self._objective.GetReponse()
        control: MaterialPropertiesControl = self.master_control.GetListOfControls()[0]

        # now we compute the modified search direction
        update = search_direction * 0.0
        for i, (filtered_mask, mapped_sensor_grad) in enumerate(list_of_gradients):
            search_direction = -KratosOA.CollectiveExpression([control.filter.ForwardFilterField(mapped_sensor_grad) * filtered_mask])
            sensor_update = KratosOA.ExpressionUtils.Scale(search_direction, alpha)
            sensor = response.list_of_sensors[i]
            self.algorithm_data.GetBufferedData()[f"{response.list_of_sensors[i].GetName()}_search_direction"] = search_direction.Clone()
            self.algorithm_data.GetBufferedData()[f"{response.list_of_sensors[i].GetName()}_update"] = sensor_update.Clone()
            sensor.AddContainerExpression("search_direction", search_direction.GetContainerExpressions()[0])
            sensor.AddContainerExpression("update", update.GetContainerExpressions()[0])
            update = KratosOA.ExpressionUtils.Collapse(update + sensor_update)

        self.algorithm_data.GetBufferedData()["control_field_update"] = update.Clone()

    @time_decorator()
    def Solve(self):
        response: DamageDetectionResponse = self._objective.GetReponse()
        control: MaterialPropertiesControl = self.master_control.GetListOfControls()[0]
        smooth_masks_list = KratosSI.MaskUtils.SmoothenMasks([sensor.GetContainerExpression(self.mask_name) for sensor in response.list_of_sensors], 3e-2)
        while not self.converged:
            with OptimizationAlgorithmTimeLogger("AlgorithmSteepestDescent",self._optimization_problem.GetStep()):
                self._InitializeIteration()

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
                list_of_masked_sensor_gradients = []
                for i, sensor in enumerate(response.list_of_sensors):
                    mask = sensor.GetContainerExpression(self.mask_name)
                    mask = Kratos.Expression.ElementExpression(self._control_field.GetContainerExpressions()[0].GetModelPart())
                    Kratos.Expression.LiteralExpressionIO.SetData(mask, 1.0)
                    filtered_mask = mask.Clone()

                    # clear all the sensors
                    response.damage_response_function.Clear()
                    # now we raw sensitivities for each sensor
                    response.damage_response_function.AddSensor(sensor)
                    raw_sensor_grad = self._control_field.Clone()

                    control.GetEmptyField().GetModelPart().ProcessInfo[KratosSI.SENSOR_NAME] = sensor.GetName()

                    response.CalculateGradient({Kratos.YOUNG_MODULUS: raw_sensor_grad})
                    masked_sensor_grad = raw_sensor_grad * KratosOA.CollectiveExpression([filtered_mask])
                    mapped_sensor_grad = control.MapGradient({Kratos.YOUNG_MODULUS: masked_sensor_grad.GetContainerExpressions()[0]})
                    list_of_masked_sensor_gradients.append((filtered_mask, mapped_sensor_grad))
                    obj_grad += KratosOA.CollectiveExpression([mapped_sensor_grad])

                    sensor.AddContainerExpression("mask", mask)
                    sensor.AddContainerExpression("filtered_mask", filtered_mask)
                    sensor.AddContainerExpression("raw_gradient", raw_sensor_grad.GetContainerExpressions()[0])
                    sensor.AddContainerExpression("mapped_sensor_grad", mapped_sensor_grad)
                    sensor.AddContainerExpression("masked_sensor_grad", masked_sensor_grad.GetContainerExpressions()[0])

                    self.algorithm_data.GetBufferedData()[f"{sensor.GetName()}_mask"] = mask.Clone()
                    self.algorithm_data.GetBufferedData()[f"{sensor.GetName()}_filtered_mask"] = filtered_mask.Clone()
                    self.algorithm_data.GetBufferedData()[f"{sensor.GetName()}_raw_gradient"] = raw_sensor_grad.Clone()
                    self.algorithm_data.GetBufferedData()[f"{sensor.GetName()}_mapped_gradient"] = mapped_sensor_grad.Clone()
                    self.algorithm_data.GetBufferedData()[f"{sensor.GetName()}_masked_sensor_grad"] = masked_sensor_grad.Clone()

                self.ComputeSearchDirection(obj_grad)

                alpha = self._line_search_method.ComputeStep()

                self.ComputeControlUpdate(alpha, list_of_masked_sensor_gradients)

                self._FinalizeIteration()

                self.Output()

                self.UpdateControl()

                self.converged = self._convergence_criteria.IsConverged()

                self._optimization_problem.AdvanceStep()

        return self.converged