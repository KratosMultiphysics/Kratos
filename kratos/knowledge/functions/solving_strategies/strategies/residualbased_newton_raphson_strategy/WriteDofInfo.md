---
type: C++ Function
title: WriteDofInfo
resource: solving_strategies/strategies/residualbased_newton_raphson_strategy.h#L1452-L1464
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/EchoInfo
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void WriteDofInfo(std::string FileName, const TSystemVectorType& rDX)`

# Calls

- [GetVariable](../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)

# Called by

- [EchoInfo](../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/EchoInfo.md)