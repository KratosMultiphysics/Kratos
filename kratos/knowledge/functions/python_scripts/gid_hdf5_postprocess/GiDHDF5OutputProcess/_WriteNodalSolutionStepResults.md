---
type: Python Method
title: _WriteNodalSolutionStepResults
resource: python_scripts/gid_hdf5_postprocess.py#L284-L311
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/WriteResultsAttribute
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _WriteNodalSolutionStepResults(self, step, time, list_of_variables_with_history):`

# Calls

- [WriteResultsAttribute](../../../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/WriteResultsAttribute.md)

# Called by

- [PrintOutput](../../../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/PrintOutput.md)