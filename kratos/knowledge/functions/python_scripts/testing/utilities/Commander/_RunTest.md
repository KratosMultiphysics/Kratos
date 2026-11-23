---
type: Python Method
title: _RunTest
resource: python_scripts/testing/utilities.py#L132-L170
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/testing/utilities/PrintTestHeader
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/PrintTestFooter
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/testing/utilities/Commander/RunPythonTests
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/RunCppTests
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/RunMPIPythonTests
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/Commander/RunMPICppTests
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _RunTest(self, test_suit_name, command, timer, working_dir=os.getcwd()): # Print test header`

# Calls

- [PrintTestHeader](../../../../../functions/python_scripts/testing/utilities/PrintTestHeader.md)
- [PrintTestFooter](../../../../../functions/python_scripts/testing/utilities/PrintTestFooter.md)

# Called by

- [RunPythonTests](../../../../../functions/python_scripts/testing/utilities/Commander/RunPythonTests.md)
- [RunCppTests](../../../../../functions/python_scripts/testing/utilities/Commander/RunCppTests.md)
- [RunMPIPythonTests](../../../../../functions/python_scripts/testing/utilities/Commander/RunMPIPythonTests.md)
- [RunMPICppTests](../../../../../functions/python_scripts/testing/utilities/Commander/RunMPICppTests.md)