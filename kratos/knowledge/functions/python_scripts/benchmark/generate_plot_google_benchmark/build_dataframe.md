---
type: Python Function
title: build_dataframe
resource: python_scripts/benchmark/generate_plot_google_benchmark.py#L54-L84
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/benchmark/generate_plot_google_benchmark/load_json
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/benchmark/generate_plot_google_benchmark/extract_context
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/benchmark/generate_plot_google_benchmark/main
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def build_dataframe(filenames: list[str]) -> tuple[pd.DataFrame, dict]:`

# Calls

- [load_json](../../../../functions/python_scripts/benchmark/generate_plot_google_benchmark/load_json.md)
- [extract_context](../../../../functions/python_scripts/benchmark/generate_plot_google_benchmark/extract_context.md)

# Called by

- [main](../../../../functions/python_scripts/benchmark/generate_plot_google_benchmark/main.md)