import dataclasses as dc
from typing import List
from typing import Union

"""_summary_
Dataclasses that define the Measurement.json file structure
LoadDataContainer is still under heavy construction and is probable to not stay as it is.
"""

@dc.dataclass
class LoadDataContainer:
    type_of_load: str = "PointLoad"
    position_of_mesh_vertex: List[float] = dc.field(default_factory=lambda: [0.0, 0.0, 0.0])
    direction_normal: List[float] = dc.field(default_factory=lambda: [0.0, 1.0, 0.0])
    strength_in_N: float = 1.0
    model_part_name: str = "load_model_part"


@dc.dataclass
class SensorDataContainer:
    type_of_sensor: str = "DISPLACEMENT"
    measurement_direction_normal: List[float] = dc.field(default_factory=lambda: [0.0, 1.0, 0.0])
    position_of_mesh_node: List[float] = dc.field(default_factory=lambda: [0.0, 0.0, 0.0])
    mesh_node_id: int = None
    measured_value: float = None

    def copy(self) -> "SensorDataContainer":
        return SensorDataContainer(type_of_sensor=self.type_of_sensor,
                                   measurement_direction_normal=[self.measurement_direction_normal[0], self.measurement_direction_normal[1], self.measurement_direction_normal[2]],
                                   position_of_mesh_node=[self.position_of_mesh_node[0], self.position_of_mesh_node[1], self.position_of_mesh_node[2]],
                                   mesh_node_id=self.mesh_node_id,
                                   measured_value=self.measured_value
                                   )


@dc.dataclass
class PerLoadCaseMeasurementDataContainer:
    load_info: LoadDataContainer = None
    sensors_infos: List[SensorDataContainer] = None


@dc.dataclass
class LoadCasesContainer:
    load_cases: List[PerLoadCaseMeasurementDataContainer] = None
