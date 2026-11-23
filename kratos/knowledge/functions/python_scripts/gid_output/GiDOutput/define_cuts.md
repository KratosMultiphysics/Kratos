---
type: Python Method
title: define_cuts
resource: python_scripts/gid_output.py#L110-L140
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/cutting_utility/CuttingUtility/CuttingUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/FindSmallestEdge
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/AddVariablesToCutModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/GenerateCut
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def define_cuts(self, model_part, cut_list):`

# Calls

- [CuttingUtility](../../../../functions/utilities/cutting_utility/CuttingUtility/CuttingUtility.md)
- [FindSmallestEdge](../../../../functions/utilities/cutting_utility/CuttingUtility/FindSmallestEdge.md)
- [AddVariablesToCutModelPart](../../../../functions/utilities/cutting_utility/CuttingUtility/AddVariablesToCutModelPart.md)
- [GenerateCut](../../../../functions/utilities/cutting_utility/CuttingUtility/GenerateCut.md)