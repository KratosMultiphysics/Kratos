---
type: C++ Function
title: GetCadJsonOutput
resource: input_output/cad_json_output.h#L87-L97
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/cad_io_modeler/CadIoModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void GetCadJsonOutput( ModelPart& rModelPart, std::string& rGeometryFile, SizeType EchoLevel = 0)`

# Calls

- [PrettyPrintJsonString](../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)

# Called by

- [SetupModelPart](../../../functions/modeler/cad_io_modeler/CadIoModeler/SetupModelPart.md)