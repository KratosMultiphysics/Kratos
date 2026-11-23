---
type: C++ Method
title: ConstructMatrixStructure
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h#L990-L1106
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void ConstructMatrixStructure( typename TSchemeType::Pointer pScheme, TSystemMatrixType& rA, ModelPart& rModelPart )`

# Calls

- [Start](../../../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../../../functions/utilities/timer/Timer/Stop.md)