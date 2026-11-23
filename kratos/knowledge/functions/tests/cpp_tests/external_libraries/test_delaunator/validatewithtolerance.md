---
type: C++ Function
title: validatewithtolerance
resource: tests/cpp_tests/external_libraries/test_delaunator.cpp#L36-L98
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/SetOutputFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/CloseOutputFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/external_libraries/test_delaunator/validate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void validatewithtolerance( const std::vector<double>& rCoordinates, const double Tolerance, const double HullArea = 0.0, const std::string& Debug = "" )`

# Calls

- [CreateNewProperties](../../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [SetOutputFile](../../../../../functions/utilities/timer/Timer/SetOutputFile.md)
- [Start](../../../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../../../functions/utilities/timer/Timer/Stop.md)
- [CloseOutputFile](../../../../../functions/utilities/timer/Timer/CloseOutputFile.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [validate](../../../../../functions/tests/cpp_tests/external_libraries/test_delaunator/validate.md)