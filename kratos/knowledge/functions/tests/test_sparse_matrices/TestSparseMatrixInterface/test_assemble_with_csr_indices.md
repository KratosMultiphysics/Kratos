---
type: Python Method
title: test_assemble_with_csr_indices
resource: tests/test_sparse_matrices.py#L219-L275
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_assemble_with_csr_indices(self): # Create model part`

# Calls

- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [assertVectorAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)