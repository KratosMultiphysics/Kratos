"""Mechanically proves the lazy-import contract: with torch and physicsnemo
made unimportable, the application and all its eagerly-importable modules
still import, and the lazily-failing entry points raise the documented
actionable ImportError."""

import subprocess
import sys

import KratosMultiphysics.KratosUnittest as KratosUnittest

_PROBE = r"""
import sys

class _Blocker:
    BLOCKED = ("torch", "physicsnemo", "torch_geometric", "torch_scatter", "pyvista", "onnxruntime", "pyacvd", "tritonclient", "gpytorch", "physicsnemo_curator", "cupy", "pxr", "tetgen", "tensordict")
    def find_module(self, name, path=None):
        if name.split(".")[0] in self.BLOCKED:
            return self
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in self.BLOCKED:
            raise ImportError(f"'{name}' is blocked for this test")
sys.meta_path.insert(0, _Blocker())
for name in list(sys.modules):
    if name.split(".")[0] in _Blocker.BLOCKED:
        del sys.modules[name]

# The application and every non-test module must import without ML packages.
import KratosMultiphysics
import KratosMultiphysics.PhysicsNeMoApplication
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.utilities import array_backend_utils
from KratosMultiphysics.PhysicsNeMoApplication.training import torch_dataset
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.processes.export import dataset_export_process
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import hybrid_initialization_process
from KratosMultiphysics.PhysicsNeMoApplication.processes import validation_metrics_process
from KratosMultiphysics.PhysicsNeMoApplication.processes.export import mesh_export_process
from KratosMultiphysics.PhysicsNeMoApplication.distributed import distributed_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import superresolution_process
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import grid_inference_process
from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
from KratosMultiphysics.PhysicsNeMoApplication.training import rollout_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import graph_inference_process
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import time_series_inference_process
from KratosMultiphysics.PhysicsNeMoApplication.physics import solver_residuals
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import sequence_inference_process
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import point_cloud_inference_process
from KratosMultiphysics.PhysicsNeMoApplication.processes.export import grid_dataset_export_process
from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import diffusion_inference_process
from KratosMultiphysics.PhysicsNeMoApplication.processes.export import cae_dataset_export_process
from KratosMultiphysics.PhysicsNeMoApplication.bridges import mapping_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges import rom_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import rom_surrogate_process
from KratosMultiphysics.PhysicsNeMoApplication.training import rom_temporal
from KratosMultiphysics.PhysicsNeMoApplication.bridges import cfd_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges import curator_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.export import curator_export_process
from KratosMultiphysics.PhysicsNeMoApplication.deployment import onnx_utils
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import onnx_inference_process
from KratosMultiphysics.PhysicsNeMoApplication.deployment import ood_guard_utils
from KratosMultiphysics.PhysicsNeMoApplication.deployment import uncertainty_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges import particle_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import particle_inference_process
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import domino_inference_process
from KratosMultiphysics.PhysicsNeMoApplication.physics import physics_informed
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import pinn_solve_process
from KratosMultiphysics.PhysicsNeMoApplication.physics import differentiable_residual
from KratosMultiphysics.PhysicsNeMoApplication.physics import diffusion_residual_operator
from KratosMultiphysics.PhysicsNeMoApplication.physics import sensitivity_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges import calculus_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import adaptive_remeshing
from KratosMultiphysics.PhysicsNeMoApplication.processes import adaptive_remesh_process
from KratosMultiphysics.PhysicsNeMoApplication.training import temporal_training
from KratosMultiphysics.PhysicsNeMoApplication.training import streaming_dataset
from KratosMultiphysics.PhysicsNeMoApplication.bridges import vfgn_bridge
from KratosMultiphysics.PhysicsNeMoApplication.distributed import graph_partition_utils
from KratosMultiphysics.PhysicsNeMoApplication.training import domino_finetune
from KratosMultiphysics.PhysicsNeMoApplication.deployment import triton_export
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import triton_inference_process
from KratosMultiphysics.PhysicsNeMoApplication.utilities import tensor_adaptor_dataset_utils
from KratosMultiphysics.PhysicsNeMoApplication.utilities import nvtx_utils
from KratosMultiphysics.PhysicsNeMoApplication.utilities import point_subsampling
from KratosMultiphysics.PhysicsNeMoApplication.deployment import geometry_guard_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges import globe_bridge
from KratosMultiphysics.PhysicsNeMoApplication.training import globe_training
from KratosMultiphysics.PhysicsNeMoApplication.training import aerojepa_pretraining
from KratosMultiphysics.PhysicsNeMoApplication.utilities import shallow_water_reference
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import tessellation, curved_tessellation, provenance, domain_mesh_builder, deformation, spatial, generate, operations, sampling
from KratosMultiphysics.PhysicsNeMoApplication.bridges import vtk_bridge
from KratosMultiphysics.PhysicsNeMoApplication.active_learning import (
    sample_io, kratos_label_strategy, query_strategies, metrology)
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.execution_backends import (
    base_backend, in_process_backend, subprocess_backend)

# The lazy entry points must fail with the documented actionable message.
for fn, expected in (
        (lambda: torch_bridge.KratosTensorToTorch(None), "pip install torch"),
        (lambda: torch_dataset._TryImportTorch(), "pip install torch"),
        (lambda: kratos_label_strategy._TryImportPhysicsNemo(), "pip install nvidia-physicsnemo"),
        (lambda: domain_mesh_builder._TryImportPhysicsNemo(), "pip install nvidia-physicsnemo"),
        (lambda: validation_metrics_process._TryImportPhysicsNemoMetrics(), "pip install nvidia-physicsnemo"),
        (lambda: distributed_utils._TryImportDistributedManager(), "pip install nvidia-physicsnemo"),
        (lambda: graph_bridge._TryImportPyG(), "pip install torch_geometric"),
        (lambda: training_utils._TryImportTorch(), "pip install torch"),
        (lambda: training_utils._TryImportStaticCapture(), "pip install -U nvidia-physicsnemo"),
        (lambda: training_utils._TryImportCheckpointing(), "pip install -U nvidia-physicsnemo"),
        (lambda: training_utils._TryImportMuon(), "pip install -U nvidia-physicsnemo"),
        (lambda: training_utils._TryImportLaunchLogger(), "pip install -U nvidia-physicsnemo"),
        (lambda: validation_metrics_process._TryImportDistributionMetrics(), "pip install nvidia-physicsnemo"),
        (lambda: validation_metrics_process._TryImportPowerSpectrum(), "pip install nvidia-physicsnemo"),
        (lambda: graph_partition_utils._TryImportAutogradCollectives(), "pip install -U nvidia-physicsnemo"),
        (lambda: rollout_utils._TryImportTorch(), "pip install torch"),
        (lambda: query_strategies._TryImportPhysicsNemo(), "pip install nvidia-physicsnemo"),
        (lambda: metrology._TryImportPhysicsNemo(), "pip install nvidia-physicsnemo"),
        (lambda: diffusion_utils._TryImportTorch(), "pip install torch"),
        (lambda: diffusion_utils._TryImportPhysicsNemoDiffusion(), "pip install nvidia-physicsnemo"),
        (lambda: diffusion_utils._TryImportPhysicsNemoDiffusionProtocol(), "pip install -U nvidia-physicsnemo"),
        (lambda: diffusion_utils._TryImportPhysicsNemoMultiDiffusion(), "pip install -U nvidia-physicsnemo"),
        (lambda: diffusion_utils._TryImportTensorDict(), "pip install nvidia-physicsnemo"),
        (lambda: diffusion_residual_operator._TryImportTorch(), "pip install torch"),
        (lambda: grid_bridge._TryImportTorch(), "pip install torch"),
        (lambda: rom_temporal._TryImportTorch(), "pip install torch"),
        (lambda: rom_temporal._TryImportSequenceModel(), "pip install nvidia-physicsnemo"),
        (lambda: cfd_bridge._TryImportPhysicsNemoCfd(), "github.com/NVIDIA/physicsnemo-cfd"),
        (lambda: cfd_bridge._TryImportCfdMetrics(), "github.com/NVIDIA/physicsnemo-cfd"),
        (lambda: cfd_bridge._TryImportPyVista(), "pip install pyvista"),
        (lambda: cfd_bridge._TryImportDominoScaling(), "github.com/NVIDIA/physicsnemo-cfd"),
        (lambda: cfd_bridge._TryImportCfdEvaluationWrappers(), "github.com/NVIDIA/physicsnemo-cfd"),
        (lambda: curator_bridge._TryImportCurator(), "github.com/NVIDIA/physicsnemo-curator"),
        (lambda: curator_bridge._TryImportCuratorMeshSinks(), "github.com/NVIDIA/physicsnemo-curator"),
        (lambda: onnx_utils._TryImportOnnxExport(), "pip install nvidia-physicsnemo"),
        (lambda: onnx_utils._TryImportOnnxRuntime(), "pip install onnxruntime"),
        (lambda: ood_guard_utils._TryImportOODGuard(), "pip install nvidia-physicsnemo"),
        (lambda: ood_guard_utils._TryImportTorch(), "pip install torch"),
        (lambda: geometry_guard_utils._TryImportTorch(), "pip install torch"),
        (lambda: geometry_guard_utils._TryImportGeometryGuardrail(), "pip install -U nvidia-physicsnemo"),
        (lambda: point_subsampling._TryImportTorch(), "pip install torch"),
        (lambda: point_subsampling._TryImportFarthestPointSampling(), "pip install -U nvidia-physicsnemo"),
        (lambda: globe_bridge._TryImportTorch(), "pip install torch"),
        (lambda: globe_bridge._TryImportGlobe(), "pip install -U nvidia-physicsnemo"),
        (lambda: globe_training._TryImportTorch(), "pip install torch"),
        (lambda: aerojepa_pretraining._TryImportTorch(), "pip install torch"),
        (lambda: aerojepa_pretraining._TryImportAeroJepa(), "pip install -U nvidia-physicsnemo"),
        (lambda: uncertainty_utils._TryImportTorch(), "pip install torch"),
        (lambda: uncertainty_utils._TryImportConcreteDropout(), "pip install nvidia-physicsnemo"),
        (lambda: uncertainty_utils._TryImportGpytorch(), "pip install gpytorch"),
        (lambda: uncertainty_utils._TryImportFieldGpHead(), "pip install -U nvidia-physicsnemo"),
        (lambda: particle_bridge._TryImportNeighborSearch(), "pip install nvidia-physicsnemo"),
        (lambda: array_backend_utils._TryImportCuPy(), "pip install cupy"),
        (lambda: physics_informed._TryImportTorch(), "pip install torch"),
        (lambda: physics_informed._TryImportPhysicsNemoSym(), "pip install nvidia-physicsnemo"),
        (lambda: differentiable_residual._TryImportTorch(), "pip install torch"),
        (lambda: sensitivity_utils._TryImportTorch(), "pip install torch"),
        (lambda: calculus_bridge._TryImportTorch(), "pip install torch"),
        (lambda: deformation._TryImportTorch(), "pip install torch"),
        (lambda: deformation._TryImportDeformers(), "pip install -U nvidia-physicsnemo"),
        (lambda: deformation._TryImportEnergies(), "pip install -U nvidia-physicsnemo"),
        (lambda: deformation._TryImportMeshDeformers(), "pip install -U nvidia-physicsnemo"),
        (lambda: operations._TryImportTorch(), "pip install torch"),
        (lambda: operations._TryImportMeshOperations(), "pip install -U nvidia-physicsnemo"),
        (lambda: sampling._TryImportTorch(), "pip install torch"),
        (lambda: sampling._TryImportSampling(), "pip install -U nvidia-physicsnemo"),
        (lambda: domain_mesh_builder._TryImportZarrIo(), "pip install -U nvidia-physicsnemo"),
        (lambda: vtk_bridge._TryImportPyVista(), "pip install pyvista"),
        (lambda: vtk_bridge._TryImportVtkReader(), "pip install -U nvidia-physicsnemo"),
        (lambda: vtk_bridge._TryImportFromPyVista(), "pip install -U nvidia-physicsnemo"),
        (lambda: vtk_bridge._TryImportTorch(), "pip install torch"),
        (lambda: spatial._TryImportTorch(), "pip install torch"),
        (lambda: spatial._TryImportSignedDistanceField(), "pip install -U nvidia-physicsnemo"),
        (lambda: generate._TryImportTorch(), "pip install torch"),
        (lambda: generate._TryImportGenerators(), "pip install -U nvidia-physicsnemo"),
        (lambda: generate._TryImportFillInterior(), "pip install -U nvidia-physicsnemo"),
        (lambda: generate.SdfPrimitives(), "pip install -U nvidia-physicsnemo"),
        (lambda: grid_bridge._TryImportGridVectorOperators(), "pip install torch nvidia-physicsnemo"),
        (lambda: calculus_bridge._TryImportPhysicsNemoCalculus(), "pip install nvidia-physicsnemo"),
        (lambda: adaptive_remeshing._TryImportTorch(), "pip install torch"),
        (lambda: temporal_training._TryImportTorch(), "pip install torch"),
        (lambda: streaming_dataset._TryImportTorch(), "pip install torch"),
        (lambda: vfgn_bridge._TryImportTorch(), "pip install torch"),
        (lambda: graph_partition_utils._TryImportTorch(), "pip install torch"),
        (lambda: vfgn_bridge._TryImportVfgn(), "pip install torch_scatter"),
        (lambda: streaming_dataset._TryImportIterableDatasetBase(), "pip install -U nvidia-physicsnemo"),
        (lambda: triton_export._TryImportTorch(), "pip install torch"),
        (lambda: domino_finetune._TryImportTorch(), "pip install torch"),
        (lambda: domino_finetune._TryImportPeft(), "pip install -U nvidia-physicsnemo"),
        (lambda: triton_inference_process._TryImportTritonClient("http"), "pip install tritonclient"),
        (lambda: adaptive_remeshing._TryImportPartitionCells(), "pip install nvidia-physicsnemo"),
        (lambda: grid_bridge._TryImportGridDerivatives(), "pip install torch nvidia-physicsnemo"),
):
    try:
        fn()
    except ImportError as err:
        assert expected in str(err), f"missing actionable hint in: {err}"
    else:
        raise AssertionError("lazy entry point did not raise ImportError")

# The array backend must degrade silently: every converted path calls this
# per step, so on a machine without cupy it has to answer False, not raise.
assert array_backend_utils.IsCuPyAvailable() is False, "IsCuPyAvailable must be False without cupy"
assert array_backend_utils.ResolveArrayModule("cupy", size_hint=10**9)[0] is __import__("numpy"), \
    "a cupy request must fall back to numpy when cupy is unavailable"

print("IMPORT_CONTRACT_OK")
"""


class TestImportContract(KratosUnittest.TestCase):
    def test_ImportsWithoutMLPackages(self):
        result = subprocess.run(
            [sys.executable, "-c", _PROBE],
            capture_output=True, text=True, timeout=240)
        self.assertEqual(result.returncode, 0, msg=f"stderr:\n{result.stderr[-3000:]}")
        self.assertIn("IMPORT_CONTRACT_OK", result.stdout)


if __name__ == '__main__':
    KratosUnittest.main()
