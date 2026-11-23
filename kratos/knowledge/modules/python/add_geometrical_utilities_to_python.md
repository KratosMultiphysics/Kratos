---
type: C++ Module
title: add_geometrical_utilities_to_python
resource: python/add_geometrical_utilities_to_python.cpp#L1-L493
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/pybind11-numpy-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/python-add-geometrical-utilities-to-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-normal-calculation-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-body-normal-calculation-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-body-distance-calculation-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-signed-distance-calculation-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-brute-force-point-locator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-binbased-fast-point-locator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-binbased-fast-point-locator-conditions-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-binbased-nodes-in-element-locator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-embedded-skin-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-geometry-tester-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-cutting-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-geometrical-transformation-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-iso-printer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-interval-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-convect-particles-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-mls-shape-functions-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-tessellation-utilities-delaunator-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-rbf-shape-functions-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-assign-master-slave-constraints-to-neighbours-utility-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [VectorizedFindHelper](../../functions/python/add_geometrical_utilities_to_python/VectorizedFindHelper.md)
- [InterpolateMeshVariableToSkinDouble](../../functions/python/add_geometrical_utilities_to_python/InterpolateMeshVariableToSkinDouble.md)
- [InterpolateMeshVariableToSkinArray](../../functions/python/add_geometrical_utilities_to_python/InterpolateMeshVariableToSkinArray.md)
- [InterpolateDiscontinuousMeshVariableToSkinDouble](../../functions/python/add_geometrical_utilities_to_python/InterpolateDiscontinuousMeshVariableToSkinDouble.md)
- [InterpolateDiscontinuousMeshVariableToSkinArray](../../functions/python/add_geometrical_utilities_to_python/InterpolateDiscontinuousMeshVariableToSkinArray.md)
- [AuxiliaryCalculateMLSShapeFunctions](../../functions/python/add_geometrical_utilities_to_python/AuxiliaryCalculateMLSShapeFunctions.md)
- [AuxiliaryCalculateMLSShapeFunctionsAndGradients](../../functions/python/add_geometrical_utilities_to_python/AuxiliaryCalculateMLSShapeFunctionsAndGradients.md)
- [AddGeometricalUtilitiesToPython](../../functions/python/add_geometrical_utilities_to_python/AddGeometricalUtilitiesToPython.md)

# Imports

- `<pybind11/numpy.h>`
- `"includes/define_python.h"`
- `"processes/process.h"`
- `"python/add_geometrical_utilities_to_python.h"`
- `"utilities/normal_calculation_utils.h"`
- `"utilities/body_normal_calculation_utils.h"`
- `"utilities/body_distance_calculation_utils.h"`
- `"utilities/signed_distance_calculation_utils.h"`
- `"utilities/brute_force_point_locator.h"`
- `"utilities/binbased_fast_point_locator.h"`
- `"utilities/binbased_fast_point_locator_conditions.h"`
- `"utilities/binbased_nodes_in_element_locator.h"`
- `"utilities/embedded_skin_utility.h"`
- `"utilities/geometry_tester.h"`
- `"utilities/cutting_utility.h"`
- `"utilities/geometrical_transformation_utilities.h"`
- `"utilities/iso_printer.h"`
- `"utilities/interval_utility.h"`
- `"utilities/convect_particles_utilities.h"`
- `"utilities/mls_shape_functions_utility.h"`
- `"utilities/tessellation_utilities/delaunator_utilities.h"`
- `"utilities/rbf_shape_functions_utility.h"`
- `"utilities/assign_master_slave_constraints_to_neighbours_utility.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)