---
type: C++ Module
title: calculate_embedded_nodal_variable_from_skin_process
resource: processes/calculate_embedded_nodal_variable_from_skin_process.h#L1-L837
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/containers-model-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-key-hash-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-flags-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/factories-linear-solver-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/elements-embedded-nodal-variable-calculation-element-simplex-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-find-intersected-geometrical-objects-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-builder-and-solvers-residualbased-block-builder-and-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-schemes-residualbased-incrementalupdate-static-scheme-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-strategies-residualbased-linear-strategy-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-intersection-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-variable-utils-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [EmbeddedNodalVariableFromSkinTypeHelperClass](../../classes/processes/calculate_embedded_nodal_variable_from_skin_process/EmbeddedNodalVariableFromSkinTypeHelperClass.md)
- [AddUnknownVariable](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/EmbeddedNodalVariableFromSkinTypeHelperClass<double>/AddUnknownVariable.md)
- [AddUnknownVariable](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/EmbeddedNodalVariableFromSkinTypeHelperClass<array_1d<double,3>>/AddUnknownVariable.md)
- [AddUnknownVariableDofs](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/EmbeddedNodalVariableFromSkinTypeHelperClass<double>/AddUnknownVariableDofs.md)
- [AddUnknownVariableDofs](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/EmbeddedNodalVariableFromSkinTypeHelperClass<array_1d<double, 3>>/AddUnknownVariableDofs.md)
- [CalculateEmbeddedNodalVariableFromSkinProcess](../../classes/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess.md)
- [CalculateEmbeddedNodalVariableFromSkinProcess](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/CalculateEmbeddedNodalVariableFromSkinProcess.md)
- [CalculateEmbeddedNodalVariableFromSkinProcess](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/CalculateEmbeddedNodalVariableFromSkinProcess-2.md)
- [Execute](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/Execute.md)
- [Clear](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/Clear.md)
- [StaticGetDefaultParameters](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/StaticGetDefaultParameters.md)
- [GetDefaultParameters](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/GetDefaultParameters.md)
- [Info](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/Info.md)
- [PrintInfo](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/PrintInfo.md)
- [PrintData](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/PrintData.md)
- [GenerateIntersectedEdgesElementsModelPart](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/GenerateIntersectedEdgesElementsModelPart.md)
- [SetObtainedEmbeddedNodalValues](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/SetObtainedEmbeddedNodalValues.md)
- [AddIntersectedElementsVariables](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddIntersectedElementsVariables.md)
- [AddIntersectedElementsModelPartDOFs](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddIntersectedElementsModelPartDOFs.md)
- [AddIntersectedElementsModelPartElements](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddIntersectedElementsModelPartElements.md)
- [SetLinearStrategy](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/SetLinearStrategy.md)
- [CalculateEmbeddedNodalVariableFromSkinProcess](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/CalculateEmbeddedNodalVariableFromSkinProcess-3.md)
- [CalculateIntersections](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/CalculateIntersections.md)
- [ClearIntersections](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/ClearIntersections.md)
- [ComputeEdgeIntersection](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/ComputeEdgeIntersection.md)
- [AddEdgeNodes](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddEdgeNodes.md)
- [pSetEdgeElementGeometry](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/pSetEdgeElementGeometry.md)
- [SetEdgePair](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/SetEdgePair.md)
- [CalculateEmbeddedNodalVariableFromSkinProcess](../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/CalculateEmbeddedNodalVariableFromSkinProcess-4.md)

# Imports

- `"containers/model.h"`
- `"includes/define.h"`
- `"includes/key_hash.h"`
- `"includes/kratos_flags.h"`
- `"factories/linear_solver_factory.h"`
- `"elements/embedded_nodal_variable_calculation_element_simplex.h"`
- `"processes/process.h"`
- `"processes/find_intersected_geometrical_objects_process.h"`
- `"solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"`
- `"solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"`
- `"solving_strategies/strategies/residualbased_linear_strategy.h"`
- `"utilities/intersection_utilities.h"`
- `"utilities/variable_utils.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)