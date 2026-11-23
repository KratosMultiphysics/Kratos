---
type: Python Method
title: ConstructListOfProcesses
resource: python_scripts/process_factory.py#L12-L76
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_CreateProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_processes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_rotated_system
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_scalar_input_to_nodes_historical
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_scalar_input_to_nodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_scalar_input_to_elements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_scalar_input_to_nodes_historical_json
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_scalar_input_to_nodes_json
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_scalar_input_to_elements_json
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_scalar_value_to_conditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_scalar_value_to_elements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_scalar_value_to_constraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_scalar_field_to_conditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_scalar_field_scalar_variable_to_conditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_scalar_field_component_to_conditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_acceleration_to_nodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_vector_variable_to_conditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_vector_variable_to_elements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_vector_variable_to_conditions_analytic_direction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_flag_process
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_fix_processes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_sub_model_part_entities_boolean_operation_process
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_average_master_slave_constraints_process
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/SolutionLoopPointOutputProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalHistoricalAvgProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalNonHistoricalAvgProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingElementalNonHistoricalAvgProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalHistoricalRMSProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalNonHistoricalRMSProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingElementalNonHistoricalRMSProcess
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ConstructListOfProcesses( self, process_list ):`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [_CreateProcesses](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_CreateProcesses.md)
- [test_assign_processes](../../../../functions/tests/test_processes/TestProcesses/test_assign_processes.md)
- [test_rotated_system](../../../../functions/tests/test_processes/TestProcesses/test_rotated_system.md)
- [test_assign_scalar_input_to_nodes_historical](../../../../functions/tests/test_processes/TestProcesses/test_assign_scalar_input_to_nodes_historical.md)
- [test_assign_scalar_input_to_nodes](../../../../functions/tests/test_processes/TestProcesses/test_assign_scalar_input_to_nodes.md)
- [test_assign_scalar_input_to_elements](../../../../functions/tests/test_processes/TestProcesses/test_assign_scalar_input_to_elements.md)
- [test_assign_scalar_input_to_nodes_historical_json](../../../../functions/tests/test_processes/TestProcesses/test_assign_scalar_input_to_nodes_historical_json.md)
- [test_assign_scalar_input_to_nodes_json](../../../../functions/tests/test_processes/TestProcesses/test_assign_scalar_input_to_nodes_json.md)
- [test_assign_scalar_input_to_elements_json](../../../../functions/tests/test_processes/TestProcesses/test_assign_scalar_input_to_elements_json.md)
- [test_assign_scalar_value_to_conditions](../../../../functions/tests/test_processes/TestProcesses/test_assign_scalar_value_to_conditions.md)
- [test_assign_scalar_value_to_elements](../../../../functions/tests/test_processes/TestProcesses/test_assign_scalar_value_to_elements.md)
- [test_assign_scalar_value_to_constraints](../../../../functions/tests/test_processes/TestProcesses/test_assign_scalar_value_to_constraints.md)
- [test_assign_scalar_field_to_conditions](../../../../functions/tests/test_processes/TestProcesses/test_assign_scalar_field_to_conditions.md)
- [test_assign_scalar_field_scalar_variable_to_conditions](../../../../functions/tests/test_processes/TestProcesses/test_assign_scalar_field_scalar_variable_to_conditions.md)
- [test_assign_scalar_field_component_to_conditions](../../../../functions/tests/test_processes/TestProcesses/test_assign_scalar_field_component_to_conditions.md)
- [test_assign_acceleration_to_nodes](../../../../functions/tests/test_processes/TestProcesses/test_assign_acceleration_to_nodes.md)
- [test_assign_vector_variable_to_conditions](../../../../functions/tests/test_processes/TestProcesses/test_assign_vector_variable_to_conditions.md)
- [test_assign_vector_variable_to_elements](../../../../functions/tests/test_processes/TestProcesses/test_assign_vector_variable_to_elements.md)
- [test_assign_vector_variable_to_conditions_analytic_direction](../../../../functions/tests/test_processes/TestProcesses/test_assign_vector_variable_to_conditions_analytic_direction.md)
- [test_assign_flag_process](../../../../functions/tests/test_processes/TestProcesses/test_assign_flag_process.md)
- [test_fix_processes](../../../../functions/tests/test_processes/TestProcesses/test_fix_processes.md)
- [test_sub_model_part_entities_boolean_operation_process](../../../../functions/tests/test_processes/TestProcesses/test_sub_model_part_entities_boolean_operation_process.md)
- [test_assign_average_master_slave_constraints_process](../../../../functions/tests/test_processes/TestProcesses/test_assign_average_master_slave_constraints_process.md)
- [SolutionLoopPointOutputProcesses](../../../../functions/tests/test_processes/SolutionLoopPointOutputProcesses.md)
- [testTimeAveragingNodalHistoricalAvgProcess](../../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalHistoricalAvgProcess.md)
- [testTimeAveragingNodalNonHistoricalAvgProcess](../../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalNonHistoricalAvgProcess.md)
- [testTimeAveragingElementalNonHistoricalAvgProcess](../../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingElementalNonHistoricalAvgProcess.md)
- [testTimeAveragingNodalHistoricalRMSProcess](../../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalHistoricalRMSProcess.md)
- [testTimeAveragingNodalNonHistoricalRMSProcess](../../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingNodalNonHistoricalRMSProcess.md)
- [testTimeAveragingElementalNonHistoricalRMSProcess](../../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/testTimeAveragingElementalNonHistoricalRMSProcess.md)