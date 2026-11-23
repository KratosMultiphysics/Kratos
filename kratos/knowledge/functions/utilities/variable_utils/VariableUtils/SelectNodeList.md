---
type: C++ Method
title: SelectNodeList
resource: utilities/variable_utils.cpp#L29-L48
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python/add_variable_utils_to_python/AddVariableUtilsToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_select_node_list
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPart::NodesContainerType VariableUtils::SelectNodeList( const DoubleVarType& Variable, const double Value, const NodesContainerType& rOriginNodes )`

# Called by

- [AddVariableUtilsToPython](../../../../functions/python/add_variable_utils_to_python/AddVariableUtilsToPython.md)
- [test_select_node_list](../../../../functions/tests/test_variable_utils/TestVariableUtils/test_select_node_list.md)