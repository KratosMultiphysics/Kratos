---
type: Python Method
title: _GetListOfVariables
resource: python_scripts/gid_hdf5_postprocess.py#L144-L149
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _GetListOfVariables(self,var_names):`

# Calls

- [GetVariable](../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [__init__](../../../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__.md)