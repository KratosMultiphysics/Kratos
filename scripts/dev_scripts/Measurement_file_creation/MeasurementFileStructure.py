import dataclasses as dc
from typing import List


@dc.dataclass
class LoadDataContainer:
    type_of_load: str = "PointLoad"
    position_of_mesh_vertex: List[float] = dc.field(default_factory=lambda: [0.0, 0.0, 0.0])
    direction_normal: List[float] = dc.field(default_factory=lambda: [0.0, 1.0, 0.0])
    strength_in_N: float = 1.0
    mesh_node_id: int = None


@dc.dataclass
class SensorDataContainer:
    type_of_sensor: str = "DISPLACEMENT"
    measurement_direction_normal: List[float] = dc.field(default_factory=lambda: [0.0, 1.0, 0.0])
    position_of_mesh_node: List[float] = dc.field(default_factory=lambda: [0.0, 0.0, 0.0])
    mesh_node_id: int = None
    measured_value: float = None


@dc.dataclass
class PerLoadCaseMeasurementDataContainer:
    load_info: LoadDataContainer = None
    sensors_infos: List[SensorDataContainer] = None


@dc.dataclass
class LoadCasesContainer:
    load_cases: List[PerLoadCaseMeasurementDataContainer] = None
