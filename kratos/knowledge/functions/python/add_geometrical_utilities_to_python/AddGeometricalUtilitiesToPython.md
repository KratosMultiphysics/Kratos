---
type: C++ Function
title: AddGeometricalUtilitiesToPython
resource: python/add_geometrical_utilities_to_python.cpp#L211-L491
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/binbased_fast_point_locator/FindPointOnMeshSimplified
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/AssignMasterSlaveConstraintsToNodes
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/kratos_python/PYBIND11_MODULE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddGeometricalUtilitiesToPython(pybind11::module &m)`

# Calls

- [FindPointOnMeshSimplified](../../../functions/utilities/binbased_fast_point_locator/FindPointOnMeshSimplified.md)
- [AssignMasterSlaveConstraintsToNodes](../../../functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/AssignMasterSlaveConstraintsToNodes.md)

# Called by

- [PYBIND11_MODULE](../../../functions/python/kratos_python/PYBIND11_MODULE.md)