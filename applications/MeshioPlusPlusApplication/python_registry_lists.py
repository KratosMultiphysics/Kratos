from typing import List

python_modelers_to_be_registered: List[str] = [
    "modelers.meshio_input_modeler.MeshioInputModeler",
    "modelers.meshio_operation_modeler.MeshioOperationModeler",
    "modelers.meshio_interpolate_modeler.MeshioInterpolateModeler",
]

python_operations_to_be_registered: List[str] = []

python_processes_to_be_registered: List[str] = [
    "meshio_output_process.MeshioOutputProcess",
]

python_stages_to_be_registered: List[str] = []

python_orchestrators_to_be_registered: List[str] = []
