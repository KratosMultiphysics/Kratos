---
type: C++ Module
title: add_processes_to_python
resource: python/add_processes_to_python.cpp#L1-L720
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/includes-define-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-parameters-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-output-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/python-add-processes-to-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-calculate-embedded-nodal-variable-from-skin-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-edge-based-gradient-recovery-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-fast-transfer-between-model-parts-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-find-nodal-h-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-find-nodal-neighbours-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-find-conditions-neighbours-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-find-global-nodal-neighbours-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-find-global-nodal-neighbours-for-entities-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-find-global-nodal-entity-neighbours-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-find-intersected-geometrical-objects-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-flux-corrected-transport-convection-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-calculate-nodal-area-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-entity-erase-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-eliminate-isolated-nodes-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-calculate-distance-to-path-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-calculate-signed-distance-to-3d-skin-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-calculate-embedded-signed-distance-to-3d-skin-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-calculate-signed-distance-to-3d-condition-skin-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-structured-mesh-generator-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-tetrahedral-mesh-orientation-check-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-variational-distance-calculation-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-levelset-convection-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-apply-constant-scalarvalue-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-apply-constant-vectorvalue-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-check-skin-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-replace-elements-and-condition-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-compute-nodal-gradient-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-compute-nodal-normal-divergence-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-assign-scalar-variable-to-entities-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-assign-scalar-input-to-entities-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-assign-scalar-field-to-entities-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-reorder-and-optimize-modelpart-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-calculate-distance-to-skin-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-calculate-discontinuous-distance-to-skin-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-apply-ray-casting-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-apply-ray-casting-interface-recognition-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-simple-mortar-mapper-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-simple-mortar-mapper-wrapper-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-skin-detection-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-sub-model-part-skin-detection-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-apply-periodic-boundary-condition-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-integration-values-extrapolation-to-nodes-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-time-averaging-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-from-json-check-result-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-set-initial-state-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-split-internal-interfaces-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-parallel-distance-calculation-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-generic-find-elements-neighbours-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-check-same-modelpart-using-skin-distance-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-calculate-nodal-distance-to-skin-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spaces-ublas-space-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-linear-solver-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [CalculateDiscontinuousEmbeddedVariableFromSkinDouble](../../functions/python/add_processes_to_python/CalculateDiscontinuousEmbeddedVariableFromSkinDouble.md)
- [CalculateDiscontinuousEmbeddedVariableFromSkinArray](../../functions/python/add_processes_to_python/CalculateDiscontinuousEmbeddedVariableFromSkinArray.md)
- [CalculateEmbeddedVariableFromSkinDouble](../../functions/python/add_processes_to_python/CalculateEmbeddedVariableFromSkinDouble.md)
- [CalculateEmbeddedVariableFromSkinArray](../../functions/python/add_processes_to_python/CalculateEmbeddedVariableFromSkinArray.md)
- [DefineSimpleMortarMapperProcess](../../functions/python/add_processes_to_python/DefineSimpleMortarMapperProcess.md)
- [AddProcessesToPython](../../functions/python/add_processes_to_python/AddProcessesToPython.md)

# Imports

- `"includes/define_python.h"`
- `"includes/kratos_parameters.h"`
- `"processes/process.h"`
- `"processes/output_process.h"`
- `"python/add_processes_to_python.h"`
- `"processes/calculate_embedded_nodal_variable_from_skin_process.h"`
- `"processes/edge_based_gradient_recovery_process.h"`
- `"processes/fast_transfer_between_model_parts_process.h"`
- `"processes/find_nodal_h_process.h"`
- `"processes/find_nodal_neighbours_process.h"`
- `"processes/find_conditions_neighbours_process.h"`
- `"processes/find_global_nodal_neighbours_process.h"`
- `"processes/find_global_nodal_neighbours_for_entities_process.h"`
- `"processes/find_global_nodal_entity_neighbours_process.h"`
- `"processes/find_intersected_geometrical_objects_process.h"`
- `"processes/flux_corrected_transport_convection_process.h"`
- `"processes/calculate_nodal_area_process.h"`
- `"processes/entity_erase_process.h"`
- `"processes/eliminate_isolated_nodes_process.h"`
- `"processes/calculate_distance_to_path_process.h"`
- `"processes/calculate_signed_distance_to_3d_skin_process.h"`
- `"processes/calculate_embedded_signed_distance_to_3d_skin_process.h"`
- `"processes/calculate_signed_distance_to_3d_condition_skin_process.h"`
- `"processes/structured_mesh_generator_process.h"`
- `"processes/tetrahedral_mesh_orientation_check.h"`
- `"processes/variational_distance_calculation_process.h"`
- `"processes/levelset_convection_process.h"`
- `"processes/apply_constant_scalarvalue_process.h"`
- `"processes/apply_constant_vectorvalue_process.h"`
- `"processes/check_skin_process.h"`
- `"processes/replace_elements_and_condition_process.h"`
- `"processes/compute_nodal_gradient_process.h"`
- `"processes/compute_nodal_normal_divergence_process.h"`
- `"processes/assign_scalar_variable_to_entities_process.h"`
- `"processes/assign_scalar_input_to_entities_process.h"`
- `"processes/assign_scalar_field_to_entities_process.h"`
- `"processes/reorder_and_optimize_modelpart_process.h"`
- `"processes/calculate_distance_to_skin_process.h"`
- `"processes/calculate_discontinuous_distance_to_skin_process.h"`
- `"processes/apply_ray_casting_process.h"`
- `"processes/apply_ray_casting_interface_recognition_process.h"`
- `"processes/simple_mortar_mapper_process.h"`
- `"processes/simple_mortar_mapper_wrapper_process.h"`
- `"processes/skin_detection_process.h"`
- `"processes/sub_model_part_skin_detection_process.h"`
- `"processes/apply_periodic_boundary_condition_process.h"`
- `"processes/integration_values_extrapolation_to_nodes_process.h"`
- `"processes/time_averaging_process.h"`
- `"processes/from_json_check_result_process.h"`
- `"processes/set_initial_state_process.h"`
- `"processes/split_internal_interfaces_process.h"`
- `"processes/parallel_distance_calculation_process.h"`
- `"processes/generic_find_elements_neighbours_process.h"`
- `"processes/check_same_modelpart_using_skin_distance_process.h"`
- `"processes/calculate_nodal_distance_to_skin_process.h"`
- `"spaces/ublas_space.h"`
- `"linear_solvers/linear_solver.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)