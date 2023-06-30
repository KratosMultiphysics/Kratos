from dataclasses import dataclass, field
from typing import List


@dataclass
class sensor_infos:
    type_of_sensor: str = "DISPLACEMENT"
    measurement_direction_normal: List[float] = field(default_factory=lambda: [0.0, 1.0, 0.0])
    position_of_mesh_node: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    mesh_node_id: int = 0
    measured_value: float = 0.0


@dataclass
class sensors_infos:
    sensors: list[sensor_infos]


@dataclass
class load_info:
    type_of_load: str = "PointLoad"
    direction_normal: List[float] = field(default_factory=lambda: [0.0, 1.0, 0.0])
    strength_in_N: float = 1.0
    position_of_mesh_vertex: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])


@dataclass
class MeasurementFile:
    load_cases: list[tuple]
