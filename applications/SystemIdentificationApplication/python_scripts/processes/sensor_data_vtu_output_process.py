import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.SystemIdentificationApplication.processes.execution_point_process import ExecutionPointProcess

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.Process:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorDataVtuOutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorDataVtuOutputProcess(model, parameters["settings"], optimization_problem)

class SensorDataVtuOutputProcess(ExecutionPointProcess):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "execution_point": "",
            "output_path"    : ""
        }""")

        self.model = model
        self.optimization_problem = optimization_problem
        parameters.ValidateAndAssignDefaults(default_settings)

        super().__init__(parameters["execution_point"].GetString())

        self.output_path = parameters["output_path"].GetString()

        self.vtu_output_dict: 'dict[Kratos.ModelPart, Kratos.VtuOutput]' = {}

    def Execute(self) -> None:
        list_of_sensors = GetSensors(self.optimization_problem)

        for sensor in list_of_sensors:
            for _, vtu_output in self.vtu_output_dict.items():
                vtu_output.ClearCellContainerExpressions()
                vtu_output.ClearNodalContainerExpressions()

            for exp_name, exp in sensor.GetContainerExpressionsMap().items():
                if exp.GetModelPart() not in self.vtu_output_dict.keys():
                    self.vtu_output_dict[exp.GetModelPart()] = Kratos.VtuOutput(exp.GetModelPart())

                vtu_output = self.vtu_output_dict[exp.GetModelPart()]
                vtu_output.AddContainerExpression(exp_name, exp)

            for _, vtu_output in self.vtu_output_dict.items():
                vtu_output.PrintOutput(f"{self.output_path}/{vtu_output.GetModelPart().FullName()}_{sensor.GetName()}")
