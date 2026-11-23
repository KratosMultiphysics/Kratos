---
type: C++ Module
title: exact_mortar_segmentation_utility
resource: utilities/exact_mortar_segmentation_utility.h#L1-L725
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/iomanip
    resolved_by: tree-sitter
    confidence: exact
  - target: external/cmath
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-mortar-classes-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-line-2d-2-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-3d-3-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-3d-4-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/integration-triangle-gauss-legendre-integration-points-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-math-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-mortar-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [ExactMortarIntegrationUtility](../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility.md)
- [SetIntegrationOrder](../../functions/utilities/exact_mortar_segmentation_utility/SetIntegrationOrder.md)
- [SetDistanceThreshold](../../functions/utilities/exact_mortar_segmentation_utility/SetDistanceThreshold.md)
- [SetEchoLevel](../../functions/utilities/exact_mortar_segmentation_utility/SetEchoLevel.md)
- [SetZeroToleranceFactor](../../functions/utilities/exact_mortar_segmentation_utility/SetZeroToleranceFactor.md)
- [SetConsiderDelaunator](../../functions/utilities/exact_mortar_segmentation_utility/SetConsiderDelaunator.md)
- [CheckAllInside](../../functions/utilities/exact_mortar_segmentation_utility/CheckAllInside.md)
- [Clipping2D](../../functions/utilities/exact_mortar_segmentation_utility/Clipping2D.md)
- [GetNormalVector2D](../../functions/utilities/exact_mortar_segmentation_utility/GetNormalVector2D.md)
- [AnglePoints](../../functions/utilities/exact_mortar_segmentation_utility/AnglePoints.md)
- [CheckPoints](../../functions/utilities/exact_mortar_segmentation_utility/CheckPoints.md)
- [FastTriagleCheck2D](../../functions/utilities/exact_mortar_segmentation_utility/FastTriagleCheck2D.md)
- [PushBackPoints](../../functions/utilities/exact_mortar_segmentation_utility/PushBackPoints.md)
- [PushBackPoints](../../functions/utilities/exact_mortar_segmentation_utility/PushBackPoints-2.md)
- [CheckInside](../../functions/utilities/exact_mortar_segmentation_utility/CheckInside.md)

# Imports

- `<iostream>`
- `<iomanip>`
- `<cmath>`
- `"includes/mortar_classes.h"`
- `"geometries/line_2d_2.h"`
- `"geometries/triangle_3d_3.h"`
- `"geometries/quadrilateral_3d_4.h"`
- `"integration/triangle_gauss_legendre_integration_points.h"`
- `"utilities/math_utils.h"`
- `"utilities/mortar_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)