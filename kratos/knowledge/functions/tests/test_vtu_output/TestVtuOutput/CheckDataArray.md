---
type: Python Method
title: CheckDataArray
resource: tests/test_vtu_output.py#L469-L474
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/check_pvtu
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckVtuFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CheckDataArray(test_class: kratos_unittest.TestCase, xml_element: ET.Element, number_of_components: int, name: str, data_type: str, data_format: 'typing.Optional[str]' = None):`

# Called by

- [check_pvtu](../../../../functions/mpi/future/tests/test_mpi_vtu_output/check_pvtu.md)
- [CheckVtuFile](../../../../functions/tests/test_vtu_output/TestVtuOutput/CheckVtuFile.md)