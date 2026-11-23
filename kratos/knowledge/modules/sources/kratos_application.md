---
type: C++ Module
title: kratos_application
resource: sources/kratos_application.cpp#L1-L456
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-variables-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-flags-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kernel-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-node-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-element-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-condition-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-constitutive-law-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-geometrical-object-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-master-slave-constraint-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/factories-standard-linear-solver-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/factories-standard-preconditioner-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/modeler-coloring-voxel-mesher-coloring-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/modeler-key-plane-generation-key-plane-generation-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/modeler-entity-generation-voxel-mesher-entity-generation-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/modeler-operation-voxel-mesher-operation-factory-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [KratosApplication](../../functions/sources/kratos_application/KratosApplication/KratosApplication.md)
- [RegisterKratosCore](../../functions/sources/kratos_application/KratosApplication/RegisterKratosCore.md)
- [DeregisterVariables](../../functions/sources/kratos_application/KratosApplication/DeregisterVariables.md)
- [DeregisterComponent](../../functions/sources/kratos_application/KratosApplication/DeregisterComponent.md)
- [DeregisterCommonComponents](../../functions/sources/kratos_application/KratosApplication/DeregisterCommonComponents.md)
- [DeregisterApplication](../../functions/sources/kratos_application/KratosApplication/DeregisterApplication.md)
- [DeregisterMappers](../../functions/sources/kratos_application/KratosApplication/DeregisterMappers.md)

# Imports

- `"includes/define.h"`
- `"includes/variables.h"`
- `"includes/kratos_flags.h"`
- `"includes/kernel.h"`
- `"includes/node.h"`
- `"includes/element.h"`
- `"includes/condition.h"`
- `"includes/constitutive_law.h"`
- `"includes/geometrical_object.h"`
- `"includes/master_slave_constraint.h"`
- `"factories/standard_linear_solver_factory.h"`
- `"factories/standard_preconditioner_factory.h"`
- `"modeler/coloring/voxel_mesher_coloring_factory.h"`
- `"modeler/key_plane_generation/key_plane_generation_factory.h"`
- `"modeler/entity_generation/voxel_mesher_entity_generation_factory.h"`
- `"modeler/operation/voxel_mesher_operation_factory.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)