---
type: C++ Function
title: exists
resource: sources/kratos_filesystem.cpp#L27-L30
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateClassFrom
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/application_generator/utils/io/CheckNameAvail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/benchmark/generate_plot_google_benchmark/load_json
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__remove_list_files
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__remove_post_results_files
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/LoadRestart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/ci_utilities/check_valid_environment_configuration_exists
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sequential_orchestrator/TestSequentialOrchestrator/test_sequential_orchestrator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool exists(const std::string& rPath)`

# Called by

- [_generateClassFrom](../../../functions/python_scripts/application_generator/applicationGenerator/ApplicationGenerator/_generateClassFrom.md)
- [CheckNameAvail](../../../functions/python_scripts/application_generator/utils/io/CheckNameAvail.md)
- [load_json](../../../functions/python_scripts/benchmark/generate_plot_google_benchmark/load_json.md)
- [__remove_list_files](../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__remove_list_files.md)
- [__remove_post_results_files](../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__remove_post_results_files.md)
- [LoadRestart](../../../functions/python_scripts/restart_utility/RestartUtility/LoadRestart.md)
- [check_valid_environment_configuration_exists](../../../functions/python_scripts/testing/ci_utilities/check_valid_environment_configuration_exists.md)
- [test_sequential_orchestrator](../../../functions/tests/test_sequential_orchestrator/TestSequentialOrchestrator/test_sequential_orchestrator.md)