---
type: C++ Method
title: GenerateCut
resource: utilities/cutting_utility.h#L144-L182
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/cutting_utility/CuttingUtility/CSR_Row_Matrix_Mod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/FirstLoop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/Create_List_Of_New_Nodes_Mod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/Calculate_Coordinate_And_Insert_New_Nodes_Mod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output/GiDOutput/define_cuts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__define_output_plane
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateCut(ModelPart& mr_model_part, ModelPart& mr_new_model_part, const array_1d<double, 3 >& versor, const array_1d<double, 3 >& Xp, int plane_number ,double tolerance_factor)`

# Calls

- [CSR_Row_Matrix_Mod](../../../../functions/utilities/cutting_utility/CuttingUtility/CSR_Row_Matrix_Mod.md)
- [FirstLoop](../../../../functions/utilities/cutting_utility/CuttingUtility/FirstLoop.md)
- [Create_List_Of_New_Nodes_Mod](../../../../functions/utilities/cutting_utility/CuttingUtility/Create_List_Of_New_Nodes_Mod.md)
- [Calculate_Coordinate_And_Insert_New_Nodes_Mod](../../../../functions/utilities/cutting_utility/CuttingUtility/Calculate_Coordinate_And_Insert_New_Nodes_Mod.md)

# Called by

- [define_cuts](../../../../functions/python_scripts/gid_output/GiDOutput/define_cuts.md)
- [__define_output_plane](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__define_output_plane.md)