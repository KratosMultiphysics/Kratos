---
type: Python Function
title: plot_benchmarks
resource: python_scripts/benchmark/generate_plot_google_benchmark.py#L138-L342
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/benchmark/generate_plot_google_benchmark/_truncate
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/benchmark/generate_plot_google_benchmark/main
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def plot_benchmarks(df: pd.DataFrame, filenames: list[str], contexts: dict, divisor: float, unit: str, output: str | None = None):`

# Calls

- [_truncate](../../../../functions/python_scripts/benchmark/generate_plot_google_benchmark/_truncate.md)

# Called by

- [main](../../../../functions/python_scripts/benchmark/generate_plot_google_benchmark/main.md)