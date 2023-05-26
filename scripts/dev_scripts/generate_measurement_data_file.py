import json
import dataclasses
from typing import List


@dataclasses.dataclass
class LoadDataContainer:
    type_of_load: str = "PointLoad",
    direction_normal: List[float] = [0, -1, 0],
    strength_in_N: float = 1,
    position_of_mesh_vertex: List[float] = [0, 0, 0],


@dataclasses.dataclass
class SensorDataContainer:
    type_of_sensor: str = "Displacement",
    measurement_direction_normal: List[float] = [0, 1, 0],
    position_of_mesh_node: List[float] = [0, 0, 0],
    mesh_node_id: float = 0,
    measured_value: float = 0,


@dataclasses.dataclass
class MeasurementDataContainer:
    description_of_load: LoadDataContainer = None,
    description_of_sensors: List[SensorDataContainer] = None,


load = LoadDataContainer(strength_in_N=2000,
                         position_of_mesh_vertex=[1.9679, 0.5, 0.],
                         direction_normal=[0, -1, 0])

sensors = [
    SensorDataContainer(measured_value=0.5,
                        position_of_mesh_node=[2.5, 0.5, 0.],
                        measurement_direction_normal=[0., -1., 0.]),
    SensorDataContainer(measured_value=0.5,
                        position_of_mesh_node=[1.4, 0.5, 0.],
                        measurement_direction_normal=[0., -1., 0.])
]

data_out = MeasurementDataContainer(description_of_load=load,
                                    description_of_sensors=sensors)


with open("MeasurementData.json", "w") as outfile:
    json.dump(dataclasses.asdict(data_out), outfile)
print(f"Written data: \n {data_out}")

file = open("MeasurementData.json")
data_in = json.load(file)
print(f"Reread data: \n {data_in}")
