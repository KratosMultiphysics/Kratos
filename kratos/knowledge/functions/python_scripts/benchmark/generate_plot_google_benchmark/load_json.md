---
type: Python Function
title: load_json
resource: python_scripts/benchmark/generate_plot_google_benchmark.py#L21-L39
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/exists
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/benchmark/generate_plot_google_benchmark/build_dataframe
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def load_json(filename: str) -> dict:`

# Calls

- [exists](../../../../functions/sources/kratos_filesystem/exists.md)

# Called by

- [build_dataframe](../../../../functions/python_scripts/benchmark/generate_plot_google_benchmark/build_dataframe.md)