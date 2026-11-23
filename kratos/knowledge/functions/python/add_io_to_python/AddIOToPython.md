---
type: C++ Function
title: AddIOToPython
resource: python/add_io_to_python.cpp#L49-L262
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/gid_io/WriteSphereMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteCircleMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteClusterMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteLocalAxesOnNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteNodalFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteLocalAxesOnNodesNonHistorical
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/PrintFlagsOnGaussPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddIOToPython(pybind11::module& m)`

# Calls

- [WriteSphereMesh](../../../functions/includes/gid_io/WriteSphereMesh.md)
- [WriteCircleMesh](../../../functions/includes/gid_io/WriteCircleMesh.md)
- [WriteClusterMesh](../../../functions/includes/gid_io/WriteClusterMesh.md)
- [WriteLocalAxesOnNodes](../../../functions/includes/gid_io/WriteLocalAxesOnNodes.md)
- [WriteNodalFlags](../../../functions/includes/gid_io/WriteNodalFlags.md)
- [WriteLocalAxesOnNodesNonHistorical](../../../functions/includes/gid_io/WriteLocalAxesOnNodesNonHistorical.md)
- [PrintFlagsOnGaussPoints](../../../functions/includes/gid_io/PrintFlagsOnGaussPoints.md)