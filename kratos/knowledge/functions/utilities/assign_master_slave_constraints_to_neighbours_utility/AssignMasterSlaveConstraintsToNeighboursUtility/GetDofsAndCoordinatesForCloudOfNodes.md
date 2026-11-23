---
type: C++ Method
title: GetDofsAndCoordinatesForCloudOfNodes
resource: utilities/assign_master_slave_constraints_to_neighbours_utility.cpp#L113-L155
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/HasDofFor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/AssignMasterSlaveConstraintsToNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignMasterSlaveConstraintsToNeighboursUtility::GetDofsAndCoordinatesForCloudOfNodes( const ResultNodesContainerType& CloudOfNodesArray, const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList, std::vector<DofPointerVectorType>& rCloudOfDofs, Matrix& rCloudOfNodesCoordinates )`

# Calls

- [HasDofFor](../../../../functions/includes/node/HasDofFor.md)

# Called by

- [AssignMasterSlaveConstraintsToNodes](../../../../functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/AssignMasterSlaveConstraintsToNodes.md)