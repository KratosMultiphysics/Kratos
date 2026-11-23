---
type: C++ Method
title: InitializeMaterial
resource: sources/constitutive_law.cpp#L552-L556
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_utilities/test_element/TestElement/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ConstitutiveLaw::InitializeMaterial(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const Vector& rShapeFunctionsValues)`

# Called by

- [Initialize](../../../../functions/tests/test_utilities/test_element/TestElement/Initialize.md)