---
type: Python Function
title: GetTemplateString
resource: python_scripts/tikz_output_process.py#L206-L258
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/tikz_output_process/TikZOutputProcess/ExecuteBeforeSolutionLoop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetTemplateString(model_part, cmyk_colors_fill, include_nodes_label = True, include_elements_label = True, include_axis = True, include_grid = True, landscape_mode = True, include_caption = False, sans_fonts = False):`

# Called by

- [ExecuteBeforeSolutionLoop](../../../functions/python_scripts/tikz_output_process/TikZOutputProcess/ExecuteBeforeSolutionLoop.md)