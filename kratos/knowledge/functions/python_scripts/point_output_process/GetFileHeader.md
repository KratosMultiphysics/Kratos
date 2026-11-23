---
type: Python Function
title: GetFileHeader
resource: python_scripts/point_output_process.py#L206-L229
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/point_output_process/IsArrayVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/point_output_process/PointOutputProcess/__SearchPoint
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetFileHeader(search_config, entity_type, entity_id, point, output_variables):`

# Calls

- [IsArrayVariable](../../../functions/python_scripts/point_output_process/IsArrayVariable.md)

# Called by

- [__SearchPoint](../../../functions/python_scripts/point_output_process/PointOutputProcess/__SearchPoint.md)