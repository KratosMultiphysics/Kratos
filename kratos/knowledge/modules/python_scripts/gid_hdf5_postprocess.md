---
type: Python Module
title: gid_hdf5_postprocess
resource: python_scripts/gid_hdf5_postprocess.py#L1-L314
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-fluiddynamicsapplication
    resolved_by: tree-sitter
    confidence: exact
  - target: external/numpy-as-np
    resolved_by: tree-sitter
    confidence: exact
  - target: external/collections
    resolved_by: tree-sitter
    confidence: exact
  - target: external/h5py
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [Factory](../../functions/python_scripts/gid_hdf5_postprocess/Factory.md)
- [GiDHDF5OutputProcess](../../classes/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess.md)
- [__init__](../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__.md)
- [IsOutputStep](../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/IsOutputStep.md)
- [PrintOutput](../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/PrintOutput.md)
- [_GetListOfVariables](../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/_GetListOfVariables.md)
- [__categorize_geometries_by_geometry_type](../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__categorize_geometries_by_geometry_type.md)
- [ExecuteInitialize](../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/ExecuteInitialize.md)
- [__WriteMesh](../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__WriteMesh.md)
- [WriteResultsAttribute](../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/WriteResultsAttribute.md)
- [_WriteNodalSolutionStepResults](../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/_WriteNodalSolutionStepResults.md)

# Imports

- `KratosMultiphysics`
- `KratosMultiphysics.FluidDynamicsApplication`
- `numpy as np`
- `collections`
- `h5py`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)