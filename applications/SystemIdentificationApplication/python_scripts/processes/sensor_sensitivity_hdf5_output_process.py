import h5py
import KratosMultiphysics as Kratos

from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.OutputProcess:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorSensitivityHDF5OutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorSensitivityHDF5OutputProcess(parameters["settings"], optimization_problem)


class SensorSensitivityHDF5OutputProcess(Kratos.OutputProcess):
    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        default_parameters = Kratos.Parameters("""{
            "sensor_group_name"  : "",
            "hdf5_file_name"     : "",
            "hdf5_dataset_prefix": "/data"
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.optimization_problem = optimization_problem

        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.hdf5_file_name = parameters["hdf5_file_name"].GetString()
        self.hdf5_dataset_prefix = parameters["hdf5_dataset_prefix"].GetString()

    def IsOutputStep(self):
        return False

    def PrintOutput(self):
        pass

    def ExecuteFinalize(self):
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        list_of_sensors = GetSensors(sensor_group_data)

        with h5py.File(self.hdf5_file_name, "w") as h5_file:
            for sensor in list_of_sensors:
                for expression_name, expression in sensor.GetContainerExpressionsMap().items():
                    expression_data = expression.Evaluate()
                    current_dataset_name = f"{self.hdf5_dataset_prefix}/{sensor.GetName()}/{expression_name}"
                    dataset = h5_file.create_dataset(current_dataset_name, data=expression_data)
                    dataset.attrs["__container_type"] = expression.__class__.__name__
                    dataset.attrs["__model_part_name"] = expression.GetModelPart().FullName()

