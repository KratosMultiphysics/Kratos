---
type: C++ Function
title: GidIO
resource: includes/gid_io.h#L144-L174
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/gid_io/InitializeResultFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/SetUpMeshContainers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/SetUpGaussPointContainers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output/GiDOutput/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/_InitializeGiDIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io/TestGidIO/test_DoubleFreeError
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/initializeOutputFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GidIO( const std::string& rDatafilename, const GiD_PostMode Mode, const MultiFileFlag UseMultipleFilesFlag, const WriteDeformedMeshFlag WriteDeformedFlag, const WriteConditionsFlag WriteConditions, const bool InitializeGaussPointContainers=true ) : mResultFileName(rDatafilename), mMeshFileName(rDatafilename), mWriteDeformed(WriteDeformedFlag), mWriteConditions(WriteConditions), mUseMultiFile(UseMultipleFilesFlag), mMode(Mode)`

# Calls

- [InitializeResultFile](../../../functions/includes/gid_io/InitializeResultFile.md)
- [SetUpMeshContainers](../../../functions/includes/gid_io/SetUpMeshContainers.md)
- [SetUpGaussPointContainers](../../../functions/includes/gid_io/SetUpGaussPointContainers.md)
- [GetData](../../../functions/sources/gid_io/GidIOBase/GetData.md)

# Called by

- [__init__](../../../functions/python_scripts/gid_output/GiDOutput/__init__.md)
- [_InitializeGiDIO](../../../functions/python_scripts/gid_output_process/GiDOutputProcess/_InitializeGiDIO.md)
- [test_DoubleFreeError](../../../functions/tests/test_gid_io/TestGidIO/test_DoubleFreeError.md)
- [initializeOutputFile](../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/initializeOutputFile.md)