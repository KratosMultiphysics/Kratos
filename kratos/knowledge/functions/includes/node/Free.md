---
type: C++ Function
title: Free
resource: includes/node.h#L877-L892
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FreeDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/openmp_utils/IsInParallel
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_node_to_python/NodeFree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/StoreFixityData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_dofs/TestDofs/testDofExport
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_graph_utilities/TestGraphUtilities/test_graph_and_connected_components_non_consecutive
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_processes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_fix_processes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void Free(const TVariableType& rDofVariable)`

# Calls

- [GetVariable](../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [FreeDof](../../../functions/includes/dof/Dof/FreeDof.md)
- [IsInParallel](../../../functions/utilities/openmp_utils/IsInParallel.md)

# Called by

- [Execute](../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute.md)
- [NodeFree](../../../functions/python/add_node_to_python/NodeFree.md)
- [StoreFixityData](../../../functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/StoreFixityData.md)
- [testDofExport](../../../functions/tests/test_dofs/TestDofs/testDofExport.md)
- [test_graph_and_connected_components_non_consecutive](../../../functions/tests/test_graph_utilities/TestGraphUtilities/test_graph_and_connected_components_non_consecutive.md)
- [test_assign_processes](../../../functions/tests/test_processes/TestProcesses/test_assign_processes.md)
- [test_fix_processes](../../../functions/tests/test_processes/TestProcesses/test_fix_processes.md)