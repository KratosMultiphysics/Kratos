---
type: C++ Method
title: Initialize
resource: tests/test_utilities/test_element.cpp#L374-L390
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/constitutive_law/ConstitutiveLaw/InitializeMaterial
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TestElement::Initialize(const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetDefaultIntegrationMethod](../../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [InitializeMaterial](../../../../../functions/sources/constitutive_law/ConstitutiveLaw/InitializeMaterial.md)