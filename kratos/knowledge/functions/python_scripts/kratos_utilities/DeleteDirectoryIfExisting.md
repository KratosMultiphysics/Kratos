---
type: Python Function
title: DeleteDirectoryIfExisting
resource: python_scripts/kratos_utilities.py#L17-L24
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/ensight_output_process/EnSightOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteDirectoryIfExistingAndEmpty
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/tikz_output_process/TikZOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtk_output_process/VtkOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtu_output_process/VtuOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_file_name_data_collector/TestFileNameDataCollector/test_GetSortedFileNamesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_generic_find_elemental_neighbours_process/TestGenericFindElementalNeighboursProcess/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_point_output_process_node
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sequential_orchestrator/TestSequentialOrchestrator/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tikz_output_process/TestTikZOutputProcess/tearDown
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/tearDown
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def DeleteDirectoryIfExisting(directory_name):`

# Called by

- [tearDown](../../../functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/tearDown.md)
- [tearDown](../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/tearDown.md)
- [__init__](../../../functions/python_scripts/ensight_output_process/EnSightOutputProcess/__init__.md)
- [DeleteDirectoryIfExistingAndEmpty](../../../functions/python_scripts/kratos_utilities/DeleteDirectoryIfExistingAndEmpty.md)
- [__init__](../../../functions/python_scripts/tikz_output_process/TikZOutputProcess/__init__.md)
- [__init__](../../../functions/python_scripts/vtk_output_process/VtkOutputProcess/__init__.md)
- [__init__](../../../functions/python_scripts/vtu_output_process/VtuOutputProcess/__init__.md)
- [tearDown](../../../functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/tearDown.md)
- [tearDown](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/tearDown.md)
- [test_GetSortedFileNamesList](../../../functions/tests/test_file_name_data_collector/TestFileNameDataCollector/test_GetSortedFileNamesList.md)
- [tearDown](../../../functions/tests/test_generic_find_elemental_neighbours_process/TestGenericFindElementalNeighboursProcess/tearDown.md)
- [test_point_output_process_node](../../../functions/tests/test_processes/TestProcesses/test_point_output_process_node.md)
- [tearDown](../../../functions/tests/test_restart/TestRestart/tearDown.md)
- [tearDown](../../../functions/tests/test_sequential_orchestrator/TestSequentialOrchestrator/tearDown.md)
- [tearDown](../../../functions/tests/test_tikz_output_process/TestTikZOutputProcess/tearDown.md)
- [tearDown](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/tearDown.md)