---
type: Python Function
title: SolutionLoopPointOutputProcesses
resource: tests/test_processes.py#L2676-L2717
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/process_factory/KratosProcessFactory/ConstructListOfProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/SetNodalValuesForPointOutputProcesses
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_processes/TestProcesses/test_point_output_process_node
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_point_output_process_element
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_point_output_process_condition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_point_output_process_interval
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_point_output_process_restart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_point_output_process_restart_with_restart_time_no_found
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_point_output_process_failed_restart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_multiple_point_output_process
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_line_output_process
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_csv_points_output_process
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def SolutionLoopPointOutputProcesses(model_part, settings, end_time, delta_time):`

# Calls

- [GetModel](../../../functions/python_scripts/project/Project/GetModel.md)
- [ConstructListOfProcesses](../../../functions/python_scripts/process_factory/KratosProcessFactory/ConstructListOfProcesses.md)
- [SetNodalValuesForPointOutputProcesses](../../../functions/tests/test_processes/SetNodalValuesForPointOutputProcesses.md)

# Called by

- [test_point_output_process_node](../../../functions/tests/test_processes/TestProcesses/test_point_output_process_node.md)
- [test_point_output_process_element](../../../functions/tests/test_processes/TestProcesses/test_point_output_process_element.md)
- [test_point_output_process_condition](../../../functions/tests/test_processes/TestProcesses/test_point_output_process_condition.md)
- [test_point_output_process_interval](../../../functions/tests/test_processes/TestProcesses/test_point_output_process_interval.md)
- [test_point_output_process_restart](../../../functions/tests/test_processes/TestProcesses/test_point_output_process_restart.md)
- [test_point_output_process_restart_with_restart_time_no_found](../../../functions/tests/test_processes/TestProcesses/test_point_output_process_restart_with_restart_time_no_found.md)
- [test_point_output_process_failed_restart](../../../functions/tests/test_processes/TestProcesses/test_point_output_process_failed_restart.md)
- [test_multiple_point_output_process](../../../functions/tests/test_processes/TestProcesses/test_multiple_point_output_process.md)
- [test_line_output_process](../../../functions/tests/test_processes/TestProcesses/test_line_output_process.md)
- [test_csv_points_output_process](../../../functions/tests/test_processes/TestProcesses/test_csv_points_output_process.md)