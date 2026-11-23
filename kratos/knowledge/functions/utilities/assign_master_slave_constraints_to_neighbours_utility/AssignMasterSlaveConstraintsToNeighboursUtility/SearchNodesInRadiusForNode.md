---
type: C++ Method
title: SearchNodesInRadiusForNode
resource: utilities/assign_master_slave_constraints_to_neighbours_utility.cpp#L53-L69
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/AssignMasterSlaveConstraintsToNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignMasterSlaveConstraintsToNeighboursUtility::SearchNodesInRadiusForNode( NodeType::Pointer pSlaveNode, double const Radius, ResultNodesContainerType& rCloudOfNodes, ResultNodesContainerType& rLocalResults)`

# Called by

- [AssignMasterSlaveConstraintsToNodes](../../../../functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/AssignMasterSlaveConstraintsToNodes.md)