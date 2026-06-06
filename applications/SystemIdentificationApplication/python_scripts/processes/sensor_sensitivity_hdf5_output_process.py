import h5py
import KratosMultiphysics as Kratos
from pathlib import Path
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.OutputProcess:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorSensitivityHDF5OutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorSensitivityHDF5OutputProcess(model, parameters["settings"], optimization_problem)


class SensorSensitivityHDF5OutputProcess(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        default_parameters = Kratos.Parameters("""{
            "sensor_group_name"  : "",
            "hdf5_file_name"     : "",
            "hdf5_dataset_prefix": "/data",
            "list_of_variables"  : []
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.optimization_problem = optimization_problem
        self.model = model

        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.hdf5_file_name = parameters["hdf5_file_name"].GetString()
        self.hdf5_dataset_prefix = parameters["hdf5_dataset_prefix"].GetString()
        self.list_of_variables = [Kratos.KratosGlobals.GetVariable(var_name) for var_name in parameters["list_of_variables"].GetStringArray()]

    def IsOutputStep(self):
        return False

    def PrintOutput(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        list_of_sensors = GetSensors(sensor_group_data)

        path = Path(self.hdf5_file_name.replace("<step>", str(self.optimization_problem.GetStep())))
        path.parent.mkdir(exist_ok=True, parents=True)

        with h5py.File(path, "w") as h5_file:
            for variable in self.list_of_variables:
                ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model[self.sensor_group_name].Nodes, variable)
                ta.CollectData()
                h5_file.create_dataset(f"{self.hdf5_dataset_prefix}/{variable.Name()}", data=ta.data)

            for sensor in list_of_sensors:
                for ta_name, ta in sensor.GetTensorAdaptorsMap().items():
                    current_dataset_name = f"{self.hdf5_dataset_prefix}/{sensor.GetName()}/{ta_name}"
                    dataset = h5_file.create_dataset(current_dataset_name, data=ta.data)
                    dataset.attrs["__container_type"] = ta.__class__.__name__

