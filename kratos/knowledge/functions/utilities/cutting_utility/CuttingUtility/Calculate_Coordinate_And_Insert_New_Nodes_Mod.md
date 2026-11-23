---
type: C++ Method
title: Calculate_Coordinate_And_Insert_New_Nodes_Mod
resource: utilities/cutting_utility.h#L579-L654
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/X0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Y0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Z0
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/cutting_utility/CuttingUtility/GenerateCut
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Calculate_Coordinate_And_Insert_New_Nodes_Mod(ModelPart& this_model_part, ModelPart& new_model_part, const DenseVector<array_1d<int, 2 > >& Position_Node, const DenseVector<int> &List_New_Nodes, const array_1d<double, 3 >& versor, const array_1d<double, 3 >& Xp, double tolerance)//,`

# Calls

- [X0](../../../../functions/includes/node/X0.md)
- [Y0](../../../../functions/includes/node/Y0.md)
- [Z0](../../../../functions/includes/node/Z0.md)

# Called by

- [GenerateCut](../../../../functions/utilities/cutting_utility/CuttingUtility/GenerateCut.md)