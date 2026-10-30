---
title: Troubleshooting and traps
keywords: troubleshooting traps symptoms errors normalization mdlus pmsh strides mpi deadlock cuda onnx sigma_data dropout
tags: [Troubleshooting_And_Traps.md]
sidebar: physicsnemo_application
summary: The failures this application's history has actually produced, each as symptom, cause and fix, grouped by where they bite - most of them silent, none of them hypothetical.
---

# Troubleshooting and traps

Every entry below was hit for real while building or using the application. The dangerous ones do not raise: they return plausible numbers. Read the *symptom* column first.

## Predictions look plausible and are wrong

| Symptom | Cause | Fix |
|---|---|---|
| a deployed model's field is off by a constant factor or offset, no error anywhere | the model was trained on normalized targets and the raw output was written onto a physical variable; a pretrained DoMINO emits dimensionless fields (a raw 0.1386 was really -609 Pa) | put `"output_normalization"` in the model card; for DoMINO set `scaling_factors_file`, `normalization` and `redimensionalize`. See [Inference](../Inference/Inference.html) |
| a particle rollout drifts (measured 18 %) though training converged | features were standardized at training and fed raw at deployment | the symmetric `"input_normalization"` card key; `torch_dataset.MakeNormalizationCardEntries` writes both |
| a test with a `Doubler` stand-in passes while the real model fails | a linear stand-in cannot detect a missing affine shift | use affine stand-ins and a non-unit fixture extent |
| everything works on the unit cube and nothing generalizes | on `[0, 1]^3` every length scale is one, so `normalize_coordinates` is the identity | test on a stretched box |
| a PINN converges to a tiny loss with half the amplitude | on a planar cloud the out-of-plane second derivative is unconstrained and the network cancels the source with curvature in `z` | `dim=2` on the builtin PDEs; the process normalizes coordinates *inside* the model so autodiff stays physical |
| a diffusion ensemble's mean is worse than its coarse input | fields of standard deviation 1e-3 fed raw to an EDM loss with `sigma_data = 0.5` - the noise wins | scale fields to `sigma_data` and un-scale everything reported |
| MC dropout returns zero uncertainty | dropout layers are silent no-ops in eval mode | the `"mc_dropout"` method re-enables them; if you sample by hand, switch them to train mode |
| an ensemble's spread stays small while the error grows tenfold | all members are wrong *together* outside the training range | this is expected; use the OOD guard for the input side, and calibration metrics to size the bars |

## Files and formats

| Symptom | Cause | Fix |
|---|---|---|
| `Module.save` refuses the file name | physicsnemo requires the `.mdlus` extension | `SaveTrainedModel` enforces it with a clear error |
| `MeshReader` finds no samples in a folder full of meshes | the reader globs `**/*.pmsh`; `Mesh.save(prefix)` adds no extension | name exports `*.pmsh` (the exporter does) |
| a mesh series trains in the wrong order | the reader sorts lexicographically, `mesh_10` before `mesh_2` | `"zero_pad_steps"` on `mesh_export_process` |
| a DoMINO datapipe misreads the global parameters | `Kratos.Parameters.keys()` is alphabetical, not insertion-ordered | `"global_params_order"` is mandatory when order matters |
| an FNO will not export to ONNX | `aten::fft_rfftn` is unsupported by both exporters | distill into an MLP or convolutional model for serving; the thermal lifecycle example does |
| `torch.onnx.export` fails on import | torch's exporter needs `onnxscript` | `pip install onnxscript` |
| `onnxruntime-gpu` installs but the CUDA provider is missing, or vice versa | `onnxruntime` and `onnxruntime-gpu` overwrite each other file for file | one per environment |
| `Kratos.Vector` holds garbage after being built from a numpy slice | `Kratos.Vector` ignores numpy strides | `numpy.ascontiguousarray` first |

## Meshes and geometry

| Symptom | Cause | Fix |
|---|---|---|
| gaps in a tessellated hexahedral mesh | the legacy `"fan"` mode splits shared faces inconsistently on unstructured numbering | the default `"smallest_id_diagonal"` |
| a gradient or divergence computed on a volume mesh is wrong | upstream's DEC gradient and divergence are silently wrong on codimension-0 meshes | `bridges.calculus_bridge` refuses them; use the LSQ backend |
| a signed distance is positive inside the body | upstream's boundary extraction winds triangles inconsistently (signed volume 0) | `spatial.BoundarySurface` re-orients; `marching_cubes` output is consistently wound |
| `remesh` returns a mesh with no fields | it drops all point and cell data by design | re-attach or re-interpolate |
| `mesh_implicit_domain` returns a worse mesh, or trips its coverage guard | it differentiates the level set to project boundary vertices and breaks under `torch.no_grad()` | the wrapper forces `enable_grad`; do not wrap the call yourself |
| a cube trips the implicit mesher's coverage guard | sharp features at moderate cell size | use `feature_points`, or a smoother primitive |
| an RBF deformation raises a singular system | thin-plate splines are singular for coplanar control points | `polynomial=False`, or non-coplanar controls |
| `connectivity_param` gives an empty particle graph | its default is far below the element size | match it to the particle spacing |
| a bistride graph model dies with "mat1 and mat2 must have the same dtype" | `graph.pos` feeds the edge MLPs and was float64 while the model is float32 | cast positions to the model dtype |

## Solvers, residuals, adjoints

| Symptom | Cause | Fix |
|---|---|---|
| a residual-driven remesh only coarsens | the residual of a *converged* solve is machine zero; only the surrogate's imperfect state carries a refinement signal | score the surrogate's field, not the solver's |
| sensitivities wrong by six orders of magnitude, condition number 1e18 | `analysis.Finalize()` released the DOF fixities; the adjoint was assembled on an unconstrained system | re-fix the supports before assembling |
| a shape sensitivity is zero | small-displacement elements read only `X0`; perturbing `X` changes nothing (loads read `X`) | perturb both, restore by writing the saved value back rather than `+h, -2h, +h` |
| an objective disagrees with Kratos's response by exactly the node count | Kratos's `point_temperature` *averages* over the part; `weighted_sum` sums | match the normalization |
| a physics loss on a grid scores an FD-exact quadratic as 3.33 instead of 0 | upstream stencils are wrong on the outermost shell of non-periodic fields | `boundary_trim` on the grid loss modes |
| a PINN result changes between `"auto"` and `"cpu"` | PINN training is device-sensitive | pin `"cpu"` for reproducibility |
| a first Newton iteration is *slower* after hybrid initialization | a rough seed closer to the solution loses to a cold start; Newton rewards smooth seeds | smooth the prediction, or use `output_interval` to seed less often |

## MPI

| Symptom | Cause | Fix |
|---|---|---|
| the MPI suite hangs after an unrelated change | a collective (`GatherFieldToRank0`, `GlobalNumberOfNodes()`, a rank-local precondition skipping one) inside a rank-0 guard; it "worked" until import timing changed | every collective runs on every rank; reduce preconditions with `SumAll` first |
| the co-simulation surrogate deadlocks under MPI | both wrappers used the rank-zero communicator; rank 0 picked the serial mapper, others the MPI one | keep one side distributed |
| `Metis` hangs the peer rank | it cannot partition an element-free mesh | `partition_mdpa: false` for point clouds |
| rendezvous fails between two test classes | `MASTER_PORT` left by an earlier class | set it explicitly; do not destroy the process group in `tearDownClass` |
| gathering FSDP2 DTensors segfaults in gloo | `fully_shard()` built a CUDA mesh because a GPU was visible | pass an explicit CPU `DeviceMesh` |
| a checkpoint saved from an FSDP2 model will not load | each rank wrote its own shard as a DTensor | `SaveTrainedModel` gathers to full tensors on rank 0 |

## Environment

| Symptom | Cause | Fix |
|---|---|---|
| an edited module still runs the old code | `bin/Release/KratosMultiphysics/...` holds copies made at install time | copy the file there, or rebuild |
| an application imports fine yesterday and segfaults today | ABI-stale install after a reconfigure with a shorter application list | see [Installation and environment](Installation_And_Environment.html) |
| `libmkl_rt.so.2: cannot open shared object file` inside a structural adjoint | MKL directories missing from `LD_LIBRARY_PATH` | add them |
| `pip install` refuses | PEP 668 externally managed interpreter | `--break-system-packages`, or a venv with a `.pth` to the user site |
| a first CUDA forward of a graph model fails in `scatter_add` | the PyG graph was never moved to the device | fixed in the process; if you build graphs by hand, `.to(device)` them |
| `nvidia-physicsnemo-cfd` cannot be pip-installed | it is not on PyPI | `pip install git+https://github.com/NVIDIA/physicsnemo-cfd` |
| a physicsnemo upgrade turns green with more skips | a renamed upstream symbol turned test classes into skips | compare skip counts across the upgrade |

If a symptom is not here and you find its cause, add the row: this page is the application's memory.
