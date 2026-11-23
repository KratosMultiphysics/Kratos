---
type: Python Method
title: testFlagAndEqual
resource: tests/test_flags.py#L99-L128
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/AsFalse
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsDefined
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def testFlagAndEqual(self): # both true | both false | opposite sets | first true | first false | second true | second false`

# Calls

- [AsFalse](../../../../functions/containers/flags/AsFalse.md)
- [IsDefined](../../../../functions/containers/flags/IsDefined.md)