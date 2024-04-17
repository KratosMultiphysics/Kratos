from pathlib import Path
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToJson

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return SensorPlacementOutputProcess(model, parameters["settings"], optimization_problem)

class SensorPlacementOutputProcess(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        default_settings = Kratos.Parameters("""{
            "model_part_name"        : "PLEASE_PROVIDE_A_MODEL_PART_NAME",
            "output_file_name"       : "",
            "sensor_status_threshold": 0.6,
            "output_every_step"      : false
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.output_file_name = settings["output_file_name"].GetString()
        self.sensor_status_threshold = settings["sensor_status_threshold"].GetDouble()
        self.output_every_step = settings["output_every_step"].GetBool()

        self.optimization_problem = optimization_problem

    def IsOutputStep(self) -> bool:
        return self.output_every_step

    def PrintOutput(self) -> None:
        output_file = Path(self.output_file_name.replace("<step>", str(self.optimization_problem.GetStep())))
        sensor_data = ComponentDataView("sensor", self.optimization_problem)
        PrintSensorListToJson(output_file, [sensor_data.GetUnBufferedData()["list_of_sensors"][i] for i, node in enumerate(self.model_part.Nodes) if node.GetValue(KratosDT.SENSOR_STATUS) > self.sensor_status_threshold])

    def ExecuteFinalize(self) -> None:
        output_file = Path(self.output_file_name.replace("<step>", "final"))
        sensor_data = ComponentDataView("sensor", self.optimization_problem)
        PrintSensorListToJson(output_file, [sensor_data.GetUnBufferedData()["list_of_sensors"][i] for i, node in enumerate(self.model_part.Nodes) if node.GetValue(KratosDT.SENSOR_STATUS) > self.sensor_status_threshold])
