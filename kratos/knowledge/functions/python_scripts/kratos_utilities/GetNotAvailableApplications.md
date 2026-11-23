---
type: Python Function
title: GetNotAvailableApplications
resource: python_scripts/kratos_utilities.py#L107-L111
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/GetListOfAvailableApplications
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/KratosUnittest/TestCase/skipTestIfApplicationsNotAvailable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/skipIfApplicationsNotAvailable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetNotAvailableApplications(*application_names):`

# Calls

- [GetListOfAvailableApplications](../../../functions/python_scripts/kratos_utilities/GetListOfAvailableApplications.md)

# Called by

- [skipTestIfApplicationsNotAvailable](../../../functions/python_scripts/KratosUnittest/TestCase/skipTestIfApplicationsNotAvailable.md)
- [skipIfApplicationsNotAvailable](../../../functions/python_scripts/KratosUnittest/skipIfApplicationsNotAvailable.md)