import h5py
import numpy
import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

class SensorIO:
    def __init__(self, model_part: Kratos.ModelPart, h5_file: h5py.File, prefix: str) -> None:
        self.__model_part = model_part
        self.__h5_file = h5_file
        self.__prefix = prefix

    def Read(self, sensor: KratosDT.Sensors.Sensor) -> None:
        h5_sensor_data = self.__h5_file[f"{self.__prefix}/{sensor.GetName()}"]
        if "nodal" in h5_sensor_data.keys() and isinstance(h5_sensor_data["nodal"], h5py.Group):
            container_data_group = h5_sensor_data["nodal"]
            exp_type = Kratos.Expression.NodalExpression
            self.__ReadContainerData(sensor, container_data_group, exp_type)

        if "condition" in h5_sensor_data.keys() and isinstance(h5_sensor_data["condition"], h5py.Group):
            container_data_group = h5_sensor_data["condition"]
            exp_type = Kratos.Expression.ConditionExpression
            self.__ReadContainerData(sensor, container_data_group, exp_type)

        if "element" in h5_sensor_data.keys() and isinstance(h5_sensor_data["element"], h5py.Group):
            container_data_group = h5_sensor_data["element"]
            exp_type = Kratos.Expression.ElementExpression
            self.__ReadContainerData(sensor, container_data_group, exp_type)

    def Write(self, sensor: KratosDT.Sensors.Sensor) -> None:
        pass

    def __ReadContainerData(self, sensor: KratosDT.Sensors.Sensor, container_data_group: h5py.Group, exp_type: 'typing.Type[ContainerExpressionTypes]') -> None:
        exp = exp_type(self.__model_part)
        for k, v in container_data_group.items():
            if not k.endswith("_partition") and isinstance(v, h5py.Dataset):
                numpy_data = numpy.array(v[:])
                Kratos.Expression.CArrayExpressionIO.Read(exp, numpy_data)
                sensor.AddContainerExpression(k, exp.Clone())
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Read {sensor.GetName()} {k} variable data.")

class OpenSensorFile:
    """@brief A context responsible for managing the lifetime of sensor HDF5 files."""

    def __init__(self, model_part: Kratos.ModelPart, h5_file_name: str, prefix: str, access_mode: str):
        self.access_mode = access_mode
        self.h5_file_name = h5_file_name
        self.prefix = prefix
        self.model_part = model_part
        self.__h5_file: 'typing.Optional[SensorIO]' = None

    def __enter__(self) -> SensorIO:
        self.__h5_file = h5py.File(self.h5_file_name, self.access_mode)
        return SensorIO(self.model_part, self.__h5_file, self.prefix)

    def __exit__(self, exit_type, exit_value, exit_traceback) -> None:
        if self.__h5_file is not None:
            self.__h5_file.close()
            self.__h5_file = None

if __name__ == "__main__":
    import KratosMultiphysics.StructuralMechanicsApplication
    from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
    mdpa_file_name = "/mnt/suneth/data/PostDoc/6_Simulation_Data/24_heat_maps/24_cluster_redundancy_campaign_medium/simulations/model_file"
    model = Kratos.Model()
    model_part = model.CreateModelPart("test")
    Kratos.ModelPartIO(mdpa_file_name, Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(model_part)

    sensor_data_file = "/mnt/suneth/data/PostDoc/6_Simulation_Data/24_heat_maps/24_cluster_redundancy_campaign_medium/simulations/sensor_placement/sensor_data.json"
    with open(sensor_data_file, "r") as file_inp:
        list_of_sensors = GetSensors(model_part, Kratos.Parameters(file_inp.read())["list_of_sensors"].values())

    h5_file_name = "/mnt/suneth/data/PostDoc/6_Simulation_Data/24_heat_maps/24_cluster_redundancy_campaign_medium/simulations/sensor_placement/sensor_sensitivities.h5"
    with OpenSensorFile(model_part, h5_file_name, "/Structure.all_nodes_elements_model_part", "r") as sensor_io:
        for sensor in list_of_sensors:
            sensor_io.Read(sensor)
