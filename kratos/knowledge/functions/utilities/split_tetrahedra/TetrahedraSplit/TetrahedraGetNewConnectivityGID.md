---
type: C++ Method
title: TetrahedraGetNewConnectivityGID
resource: utilities/split_tetrahedra.h#L150-L164
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
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

`static inline void TetrahedraGetNewConnectivityGID( const int tetra_index, const int t[56], const int aux_ids[11], int* id0, int* id1, int* id2, int* id3)`

# Called by

- [CalculateDiscontinuousShapeFunctions](../../../../functions/utilities/discont_utils/DiscontinuousShapeFunctionsUtilities/CalculateDiscontinuousShapeFunctions.md)
- [GenerateDivision](../../../../functions/utilities/divide_tetrahedra_3d_4/DivideTetrahedra3D4<TPointType>/GenerateDivision.md)
- [CalculateTetrahedraEnrichedShapeFuncions](../../../../functions/utilities/enrichment_utilities/CalculateTetrahedraEnrichedShapeFuncions.md)
- [CalculateEnrichedShapeFuncions](../../../../functions/utilities/enrichment_utilities/EnrichmentUtilities/CalculateEnrichedShapeFuncions.md)
- [CalculateTetrahedraEnrichedShapeFuncions](../../../../functions/utilities/enrichment_utilities_duplicate_dofs/CalculateTetrahedraEnrichedShapeFuncions.md)
- [CalculateSplitTetrahedraShapeFuncions](../../../../functions/utilities/split_tetrahedra_utilities/CalculateSplitTetrahedraShapeFuncions.md)