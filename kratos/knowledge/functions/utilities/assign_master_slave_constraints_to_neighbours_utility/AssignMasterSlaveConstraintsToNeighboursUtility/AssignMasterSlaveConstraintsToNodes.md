---
type: C++ Method
title: AssignMasterSlaveConstraintsToNodes
resource: utilities/assign_master_slave_constraints_to_neighbours_utility.cpp#L157-L237
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/SearchNodesInRadiusForNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/GetDofsAndCoordinatesForSlaveNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/GetDofsAndCoordinatesForCloudOfNodes
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/add_geometrical_utilities_to_python/AddGeometricalUtilitiesToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_master_slave_constraints_to_neighbours_process/AssignMasterSlaveConstraintsToNeighboursProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_master_slave_constraints_to_neighbours_process/AssignMasterSlaveConstraintsToNeighboursProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignMasterSlaveConstraintsToNeighboursUtility::AssignMasterSlaveConstraintsToNodes( NodesContainerType pSlaveNodes, double const Radius, ModelPart& rComputingModelPart, const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList, double const MinNumOfNeighNodes )`

# Calls

- [SearchNodesInRadiusForNode](../../../../functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/SearchNodesInRadiusForNode.md)
- [GetDofsAndCoordinatesForSlaveNode](../../../../functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/GetDofsAndCoordinatesForSlaveNode.md)
- [GetDofsAndCoordinatesForCloudOfNodes](../../../../functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/GetDofsAndCoordinatesForCloudOfNodes.md)

# Called by

- [AddGeometricalUtilitiesToPython](../../../../functions/python/add_geometrical_utilities_to_python/AddGeometricalUtilitiesToPython.md)
- [ExecuteInitialize](../../../../functions/python_scripts/assign_master_slave_constraints_to_neighbours_process/AssignMasterSlaveConstraintsToNeighboursProcess/ExecuteInitialize.md)
- [ExecuteInitializeSolutionStep](../../../../functions/python_scripts/assign_master_slave_constraints_to_neighbours_process/AssignMasterSlaveConstraintsToNeighboursProcess/ExecuteInitializeSolutionStep.md)