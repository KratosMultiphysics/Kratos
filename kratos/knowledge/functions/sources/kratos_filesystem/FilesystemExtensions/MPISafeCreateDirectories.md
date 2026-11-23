---
type: C++ Method
title: MPISafeCreateDirectories
resource: sources/kratos_filesystem.cpp#L124-L132
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/GetOutputFileName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/GetOutputFileName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/_ValidateFileName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/Save
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/CreateOutputFolder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__ValidateAndAssignOutputFolderPath
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideInputToPartitions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_filesystem/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FilesystemExtensions::MPISafeCreateDirectories(const std::filesystem::path& rPath)`

# Called by

- [GetOutputFileName](../../../../functions/input_output/ensight_output/EnSightOutput/GetOutputFileName.md)
- [GetOutputFileName](../../../../functions/input_output/vtk_output/VtkOutput/GetOutputFileName.md)
- [_ValidateFileName](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/_ValidateFileName.md)
- [Save](../../../../functions/python_scripts/project/Project/Save.md)
- [CreateOutputFolder](../../../../functions/python_scripts/restart_utility/RestartUtility/CreateOutputFolder.md)
- [__ValidateAndAssignOutputFolderPath](../../../../functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__ValidateAndAssignOutputFolderPath.md)
- [DivideInputToPartitions](../../../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitions.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_filesystem/KRATOS_TEST_CASE_IN_SUITE.md)