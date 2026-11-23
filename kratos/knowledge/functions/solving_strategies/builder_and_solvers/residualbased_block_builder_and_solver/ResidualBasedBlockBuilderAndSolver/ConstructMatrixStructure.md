---
type: C++ Method
title: ConstructMatrixStructure
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h#L1434-L1541
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

`virtual void ConstructMatrixStructure( typename TSchemeType::Pointer pScheme, TSystemMatrixType& A, ModelPart& rModelPart)`

# Calls

- [Start](../../../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../../../functions/utilities/timer/Timer/Stop.md)