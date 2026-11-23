---
type: C++ Method
title: Modeler
resource: modeler/modeler.h#L66-L73
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`explicit Modeler( Parameters ModelerParameters = Parameters()) : mParameters(ModelerParameters) , mEchoLevel( ModelerParameters.Has("echo_level") ? ModelerParameters["echo_level"].GetInt() : 0)`

# Calls

- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)