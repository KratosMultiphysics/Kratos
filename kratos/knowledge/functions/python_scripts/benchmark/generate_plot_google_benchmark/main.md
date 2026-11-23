---
type: Python Function
title: main
resource: python_scripts/benchmark/generate_plot_google_benchmark.py#L344-L356
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/benchmark/generate_plot_google_benchmark/build_dataframe
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/benchmark/generate_plot_google_benchmark/smart_time_unit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/benchmark/generate_plot_google_benchmark/print_summary_table
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/benchmark/generate_plot_google_benchmark/plot_benchmarks
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def main(filenames: list[str], output: str | None = None):`

# Calls

- [build_dataframe](../../../../functions/python_scripts/benchmark/generate_plot_google_benchmark/build_dataframe.md)
- [smart_time_unit](../../../../functions/python_scripts/benchmark/generate_plot_google_benchmark/smart_time_unit.md)
- [print_summary_table](../../../../functions/python_scripts/benchmark/generate_plot_google_benchmark/print_summary_table.md)
- [plot_benchmarks](../../../../functions/python_scripts/benchmark/generate_plot_google_benchmark/plot_benchmarks.md)