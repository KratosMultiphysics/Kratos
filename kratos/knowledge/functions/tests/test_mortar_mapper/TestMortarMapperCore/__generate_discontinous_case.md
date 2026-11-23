---
type: Python Method
title: __generate_discontinous_case
resource: tests/test_mortar_mapper.py#L219-L256
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveNodesFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __generate_discontinous_case(self, inverted):`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [RemoveNodesFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveNodesFromAllLevels.md)
- [RemoveConditionsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels.md)