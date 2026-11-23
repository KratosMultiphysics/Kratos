---
type: Python Method
title: _finalize_results
resource: python_scripts/gid_output.py#L165-L166
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/gid_io/FinalizeResults
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output/GiDOutput/write_results
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _finalize_results(self):`

# Calls

- [FinalizeResults](../../../../functions/includes/gid_io/FinalizeResults.md)

# Called by

- [write_results](../../../../functions/python_scripts/gid_output/GiDOutput/write_results.md)