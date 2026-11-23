---
type: C++ Module
title: add_containers_to_python
resource: python/add_containers_to_python.cpp#L1-L759
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/pybind11-numpy-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/pybind11-operators-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-data-value-container-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-variables-list-data-value-container-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-flags-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-variable-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-nd-data-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-flags-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-constitutive-law-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-convection-diffusion-settings-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-radiation-settings-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-quaternion-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/python-numpy-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-cfd-variables-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [FlagsOr](../../functions/python/add_containers_to_python/FlagsOr.md)
- [FlagsAnd](../../functions/python/add_containers_to_python/FlagsAnd.md)
- [FlagsSet1](../../functions/python/add_containers_to_python/FlagsSet1.md)
- [FlagsSet2](../../functions/python/add_containers_to_python/FlagsSet2.md)
- [Array1DModifier](../../classes/python/add_containers_to_python/Array1DModifier.md)
- [Resize](../../functions/python/add_containers_to_python/Array1DModifier/Resize.md)
- [MoveSlice](../../functions/python/add_containers_to_python/Array1DModifier/MoveSlice.md)
- [VariableIndexingUtility](../../functions/python/add_containers_to_python/VariableIndexingUtility.md)
- [DataValueContainerIndexingUtility](../../functions/python/add_containers_to_python/DataValueContainerIndexingUtility.md)
- [AddNDData](../../functions/python/add_containers_to_python/AddNDData.md)
- [AddContainersToPython](../../functions/python/add_containers_to_python/AddContainersToPython.md)

# Imports

- `<pybind11/numpy.h>`
- `<pybind11/operators.h>`
- `"containers/data_value_container.h"`
- `"containers/variables_list_data_value_container.h"`
- `"containers/flags.h"`
- `"containers/variable.h"`
- `"containers/nd_data.h"`
- `"includes/define_python.h"`
- `"includes/kratos_flags.h"`
- `"includes/constitutive_law.h"`
- `"includes/convection_diffusion_settings.h"`
- `"includes/radiation_settings.h"`
- `"utilities/quaternion.h"`
- `"python/numpy_utils.h"`
- `"includes/cfd_variables.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)