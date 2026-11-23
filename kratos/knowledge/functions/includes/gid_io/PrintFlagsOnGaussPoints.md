---
type: C++ Function
title: PrintFlagsOnGaussPoints
resource: includes/gid_io.h#L1342-L1356
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_gauss_point_container/GidGaussPointsContainer/PrintFlagsResults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/add_io_to_python/AddIOToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__write_elemental_conditional_flags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/writeResults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PrintFlagsOnGaussPoints( const Kratos::Flags& rFlag, const std::string& rFlagName, const ModelPart& rModelPart, const double SolutionTag )`

# Calls

- [Start](../../../functions/utilities/timer/Timer/Start.md)
- [PrintFlagsResults](../../../functions/includes/gid_gauss_point_container/GidGaussPointsContainer/PrintFlagsResults.md)
- [Stop](../../../functions/utilities/timer/Timer/Stop.md)

# Called by

- [AddIOToPython](../../../functions/python/add_io_to_python/AddIOToPython.md)
- [__write_elemental_conditional_flags](../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_elemental_conditional_flags.md)
- [writeResults](../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/writeResults.md)