---
type: C++ Method
title: TetrahedraSplitMode
resource: utilities/split_tetrahedra.h#L91-L136
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_split_tetrahedra/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/discont_utils/DiscontinuousShapeFunctionsUtilities/CalculateDiscontinuousShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_tetrahedra_3d_4/DivideTetrahedra3D4<TPointType>/GenerateDivision
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/enrichment_utilities/CalculateTetrahedraEnrichedShapeFuncions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/enrichment_utilities/EnrichmentUtilities/CalculateEnrichedShapeFuncions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/enrichment_utilities_duplicate_dofs/CalculateTetrahedraEnrichedShapeFuncions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_tetrahedra_utilities/CalculateSplitTetrahedraShapeFuncions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void TetrahedraSplitMode( int aux_ids[11], int edge_ids[6])`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_split_tetrahedra/KRATOS_TEST_CASE_IN_SUITE.md)
- [CalculateDiscontinuousShapeFunctions](../../../../functions/utilities/discont_utils/DiscontinuousShapeFunctionsUtilities/CalculateDiscontinuousShapeFunctions.md)
- [GenerateDivision](../../../../functions/utilities/divide_tetrahedra_3d_4/DivideTetrahedra3D4<TPointType>/GenerateDivision.md)
- [CalculateTetrahedraEnrichedShapeFuncions](../../../../functions/utilities/enrichment_utilities/CalculateTetrahedraEnrichedShapeFuncions.md)
- [CalculateEnrichedShapeFuncions](../../../../functions/utilities/enrichment_utilities/EnrichmentUtilities/CalculateEnrichedShapeFuncions.md)
- [CalculateTetrahedraEnrichedShapeFuncions](../../../../functions/utilities/enrichment_utilities_duplicate_dofs/CalculateTetrahedraEnrichedShapeFuncions.md)
- [CalculateSplitTetrahedraShapeFuncions](../../../../functions/utilities/split_tetrahedra_utilities/CalculateSplitTetrahedraShapeFuncions.md)