---
type: C++ Method
title: ModelPartIO
resource: sources/model_part_io.cpp#L33-L67
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/SetOutputFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPartIO::ModelPartIO(std::filesystem::path const& Filename, const Flags Options) : mNumberOfLines(1) , mBaseFilename(Filename) , mOptions(Options)`

# Calls

- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [SetOutputFile](../../../../functions/utilities/timer/Timer/SetOutputFile.md)