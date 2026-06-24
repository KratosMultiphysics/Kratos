import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.Process:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorSensitivityAggregationProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorSensitivityAggregationProcess(model, parameters["settings"], optimization_problem)

class SensorSensitivityAggregationProcess(Kratos.Process):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.Process.__init__(self)

        default_parameters = Kratos.Parameters("""{
            "sensor_group_name"         : "",
            "tensor_adaptor_names"      : [],
            "output_tensor_adaptor_name": "<TEST_ANALYSIS_NAME>_<TENSOR_ADAPTOR_NAME>"
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.optimization_problem = optimization_problem
        self.model = model

        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.tensor_adaptor_names = parameters["tensor_adaptor_names"].GetStringArray()
        self.output_tensor_adaptor_name = parameters["output_tensor_adaptor_name"].GetString()

    def ExecuteFinalizeSolutionStep(self):
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        sensor_model_part = self.model[self.sensor_group_name]

        list_of_sensors = GetSensors(sensor_group_data)
        current_sensor_name: str = sensor_model_part.ProcessInfo[KratosSI.SENSOR_NAME]
        test_analysis_name: str = sensor_model_part.ProcessInfo[KratosSI.TEST_ANALYSIS_NAME]

        found_sensor = False
        for sensor in list_of_sensors:
            if current_sensor_name == sensor.GetName():
                found_sensor = True
                for tensor_adaptor_name in self.tensor_adaptor_names:
                    tensor_adaptor: Kratos.TensorAdaptors.DoubleTensorAdaptor = sensor_group_data.GetUnBufferedData().GetValue(tensor_adaptor_name)
                    current_output_tensor_adaptor_name = self.output_tensor_adaptor_name.replace("<TEST_ANALYSIS_NAME>", test_analysis_name)
                    current_output_tensor_adaptor_name = current_output_tensor_adaptor_name.replace("<TENSOR_ADAPTOR_NAME>", tensor_adaptor_name)
                    sensor.AddTensorAdaptor(current_output_tensor_adaptor_name, tensor_adaptor)
                break

        if not found_sensor:
            raise RuntimeError(f"The sensor \"{current_sensor_name}\" not found.")


