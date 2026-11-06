---
title: Versions and compatibility
keywords: physicsnemo version 2.2 extras dependencies compatibility upgrade
tags: [Versions_And_Compatibility.md]
sidebar: physicsnemo_application
summary: The physicsnemo release this application tracks, the optional extras and what each pulls in, what changed between 2.1 and 2.2 that the bridge absorbed, and how to upgrade without being fooled by green tests.
---

# Versions and compatibility

This application is exercised against **physicsnemo 2.2.0**. Upstream's 2.2.1 is a fix-only release (Python 3.14 packaging and annotation introspection) and the 2.3.0 changelog is still empty, so 2.2 is the current surface. Where a page says an API "does not exist" or "raises", it means in 2.2.

## What is installed where this was written

The versions the test suite, the notebooks and the Examples cases were last run against. None of them is a *requirement* - the application declares no ML dependency at all (see [Overview](../General/Overview.html)) - but if something behaves differently for you, this is the baseline to compare with.

| Package | Version here | Role |
|---|---|---|
| `nvidia-physicsnemo` | 2.2.0 | everything ML |
| `torch` | 2.13.0 (CUDA 13.0 build) | everything ML |
| `tensordict` | 0.14.0 | the mesh tensorclass and datapipe samples; 2.2 requires `tensordict[zarr]>=0.14` |
| `warp-lang` | 1.16.0 | Warp kernels: ball query, SDF, remeshing, the FFD deformer; 2.2 requires `>=1.14` |
| `zarr` | 3.3.0 | mesh and datapipe Zarr I/O |
| `torch-geometric`, `torch_scatter` | 2.8.0, 2.1.2 | the graph processes (`PYG_AVAILABLE` upstream checks both) |
| `onnx`, `onnxscript`, `onnxruntime` | 1.22.0, 0.7.1, 1.29.0 | ONNX export (torch's exporter needs `onnxscript`) and CPU inference; the GPU build was verified once at 1.29.0 in a throwaway environment |
| `gpytorch` | 1.15.2 | the GP uncertainty head (`uq-extras`) |
| `pyvista` | 0.48.4 | mesh rendering, `cfd_bridge`, the notebooks' figures |
| `nvidia-physicsnemo-cfd` | 0.0.3a0 (source clone) | the CFD metric registry and hybrid initialization; not on PyPI |
| `tritonclient` | 2.71.0 | the Triton client (needs its protocol extra) |
| `cupy-cuda13x` | 14.1.1 | the opt-in array backend (upstream's `cu13` extra pins `<14`; 14.1.1 works here) |
| `usd-core` | 26.8 | OpenUSD export |
| `tetgen` | 0.8.4 | the opt-in exact-boundary tetrahedral fill (AGPL) |

Absent here, and the reason a few tests self-skip: `torch_sparse` and `dgl` (GraphCast), `onnxruntime-gpu` (kept out of the reference environment because it overwrites the CPU build file for file).

## The extras, and what each pulls in

`pip install "nvidia-physicsnemo[<extra>]"`. Read from the installed package's metadata:

| Extra | Pulls in | You want it for |
|---|---|---|
| `cu12` / `cu13` | `cuml`, `cupy`, `nvidia-dali`, `pylibraft` for that CUDA major, plus torch/torchvision | GPU-side data loading and clustering kernels; pick the one matching your driver |
| `mesh-extras` | `matplotlib`, `pyvista`, `vtk` | `physicsnemo.mesh` visualization and pyvista conversion |
| `datapipes-extras` | `dask`, `netcdf4`, `tfrecord`, `xarray`, `zarr` | the climate and Zarr readers |
| `gnns` | `torch-geometric`, `torch-scatter`, `torch-sparse`, `torch-cluster`, plus `pyvista`, `vtk`, `stl`, `scipy`, `mlflow`, `wandb` | every graph model; GraphCast needs `torch-sparse` |
| `model-extras`, `nn-extras`, `utils-extras` | `scipy`, `stl`, `vtk`, `mlflow`, `wandb`, `line-profiler` | logging and profiling helpers, STL I/O |
| `sym` | `sympy` | `physicsnemo.sym` (the module itself is bundled; only SymPy is extra) |
| `uq-extras` | `gpytorch` | the variational GP heads |
| `natten-cu12` / `natten-cu13` | `natten` | neighborhood attention layers |
| `transformer-engine-cu12` / `-cu13` | `transformer-engine[core,pytorch]` | `use_te=True` on Transolver-family models (fp8 attention) |

The core dependencies (always installed) include `torch`, `tensordict[zarr]`, `warp-lang`, `h5py`, `onnx`, `hydra-core`/`omegaconf`, `s3fs`/`fsspec`, `timm`, `nvtx`, `treelib` and `einops`.

## What changed from 2.1 to 2.2, as seen from this bridge

Every one of these either broke something here or would have. The lesson that recurs: **a layout change that keeps shapes identical passes every test whose fixture is symmetric.**

| Change in 2.2 | Effect on this application |
|---|---|
| Mesh-calculus gradients are derivative-first `(N, D, C)` from *every* backend (2.1's LSQ was channel-major) | `bridges.calculus_bridge` normalizes all backends to its own stable `(N, C, D)` contract and refuses physicsnemo older than 2.2 rather than guessing. All 469 tests of the time passed with the wrong layout because the canary gradient was a symmetric matrix; the canary is now asymmetric and non-square |
| `ShardTensor` moved from `physicsnemo.distributed` to `physicsnemo.domain_parallel` | `distributed.domain_parallel_utils` imports from the new location |
| `remesh` is Warp-backed (`pyacvd` dropped) and its count targets output *vertices*, not cells; it raises for anything but a surface in 3-D | `bridges.mesh_bridge.adaptive_remeshing` documents the new semantics; `remesh` still drops all point and cell data |
| `Mesh.save`/`load` gained a Zarr backend (`mesh.io.to_zarr`/`from_zarr`) | `bridges.mesh_bridge.domain_mesh_builder`'s `SaveMeshZarr`/`LoadMeshZarr`, and the default backend of `CuratorExportProcess` - AI-ready Zarr without the git-only curator |
| GeoTransolver and FLARE promoted out of `experimental` (with `ShardTensor` support and activation checkpointing) | `processes.inference.point_cloud_inference_process` reaches both; the OOD guard was decoupled from GeoTransolver at the same time |
| `nn.functional.signed_distance_field` returns a 3-tuple `(sdf, hit_points, hit_faces)` | `bridges.mesh_bridge.spatial` unpacks it |
| `fill_interior` gained exact-boundary 2-D filling; its `n = 3` raises `NotImplementedError` | the tetrahedral fill in `bridges.mesh_bridge.generate` uses its own winding-number carve, with `tetgen` as the opt-in exact-boundary backend |
| `export_to_onnx_stream` no longer runs the model twice, but still exposes no `dynamic_axes` | `deployment.triton_export` keeps exporting through `torch.onnx.export` so the entity axis stays dynamic |
| `datapipes/protocols.py` rewritten (`_PrefetchResult` became `HostPayload`) | nothing here subclasses it; `training.streaming_dataset` inherits `IterableDatasetBase` and torch's `IterableDataset` together, which is unchanged |
| `integrate` gained `nan_policy`; `integrate_cell_data`/`integrate_point_data` deprecated in favour of `integrate(...)` | `bridges.calculus_bridge` uses `integrate` |
| The legacy diffusion modules (`samplers.legacy_deterministic_sampler`, `metrics.legacy_losses`, `preconditioners.legacy`) now warn they "will be deprecated in a future release" | migrated: the bridge runs on the protocol API (`WrapDiffusionModel`, `MSEDSMLoss` over an `EDMNoiseScheduler`, `samplers.sample`), and `"api": "legacy"` still reaches these modules because 2.2 ships no protocol equivalent of CorrDiff's `RegressionLoss`/`ResidualLoss` |
| `shrink_and_perturb_`, the mesh deformers (`sobolev_deform`, `shrinkwrap`, RBF, FFD), the fixed-topology energies, the grid divergence/curl/Laplacian functionals, `farthest_point_sampling`, FSDP2 checkpoint support | all of them shipped here: `TrainModel`'s warm restart, the deformation layer (`sobolev_deform`/`shrinkwrap` through `DeformMesh`, which take a `Mesh` rather than points), the fixed-topology energies, `grid_bridge`'s vector operators, `utilities.point_subsampling`'s token budgets, and `SaveTrainedModel`'s DTensor gather |
| `poisson_sample_indices_fixed` removed | unused here |

## How to check what you have

```python
import physicsnemo, torch, tensordict
print(physicsnemo.__version__, torch.__version__, tensordict.__version__)
```

The lazy-import helpers in every module produce an actionable message naming the missing package, so the fastest check of a fresh environment is simply running the small test suite and reading the skip reasons:

```bash
python3 applications/PhysicsNeMoApplication/tests/test_PhysicsNeMoApplication.py -l small
```

## Upgrading physicsnemo: compare skip counts, not pass/fail

Two things in this suite make a green run after an upgrade *insufficient* evidence:

1. **Optional-dependency gates skip green.** The graph, ONNX, Triton, GP-head, pyvista and cfd tests each self-skip when their import fails. If an upgrade renames a symbol the gate probes (the `PYG_AVAILABLE` flag, or the logger name the OOD guard capture keys on), whole classes turn into skips and the run stays green. Record the skip count before and after; it must not grow.
2. **Symmetric fixtures hide layout flips.** A unit cube hides every length scale, and a field with a symmetric gradient hides a transpose. The suite now carries asymmetric canaries and a non-unit extent for exactly this reason, but a new bridge you write should get one too.

The upstream [release notes](https://docs.nvidia.com/physicsnemo/latest/release-notes/index.html) and the repository `CHANGELOG.md` list every change per release; the table above is the subset that touched this bridge.

Back to [PhysicsNeMo Basics](Overview.html), or on to [Where things live](../General/Module_Map.html).
