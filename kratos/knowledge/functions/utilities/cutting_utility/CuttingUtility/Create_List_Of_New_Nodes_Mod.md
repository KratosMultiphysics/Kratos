---
type: C++ Method
title: Create_List_Of_New_Nodes_Mod
resource: utilities/cutting_utility.h#L510-L575
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/cutting_utility/CuttingUtility/GenerateCut
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Create_List_Of_New_Nodes_Mod(ModelPart& this_model_part, ModelPart& new_model_part, boost::numeric::ublas::compressed_matrix<int>& Coord, DenseVector<int> &List_New_Nodes, DenseVector<array_1d<int, 2 > >& Position_Node) //plane number = 1 -- inf (but the first one should be always one!`

# Called by

- [GenerateCut](../../../../functions/utilities/cutting_utility/CuttingUtility/GenerateCut.md)