---
type: C++ Module
title: embedded_skin_utility
resource: utilities/embedded_skin_utility.cpp#L1-L419
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/variable-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/embedded-skin-utility-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/modified-shape-functions-triangle-2d-3-ausas-modified-shape-functions-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/modified-shape-functions-tetrahedra-3d-4-ausas-modified-shape-functions-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/modified-shape-functions-triangle-2d-3-modified-shape-functions-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/modified-shape-functions-tetrahedra-3d-4-modified-shape-functions-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-divide-triangle-2d-3-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-divide-tetrahedra-3d-4-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GenerateSkin](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GenerateSkin.md)
- [InterpolateMeshVariableToSkin](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/InterpolateMeshVariableToSkin.md)
- [InterpolateMeshVariableToSkin](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/InterpolateMeshVariableToSkin-2.md)
- [InterpolateDiscontinuousMeshVariableToSkin](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/InterpolateDiscontinuousMeshVariableToSkin.md)
- [InterpolateDiscontinuousMeshVariableToSkin](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/InterpolateDiscontinuousMeshVariableToSkin-2.md)
- [Clear](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/Clear.md)
- [ComputeElementSkin](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/ComputeElementSkin.md)
- [RenumberAndAddSkinEntities](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/RenumberAndAddSkinEntities.md)
- [pCreateNewConditionGeometry](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/pCreateNewConditionGeometry.md)
- [pCreateNewCondition](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/pCreateNewCondition.md)
- [GetConditionType](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<2>/GetConditionType.md)
- [GetConditionType](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<3>/GetConditionType.md)
- [SetSkinEntitiesProperties](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetSkinEntitiesProperties.md)
- [ElementIsSplit](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/ElementIsSplit.md)
- [SetDistancesVector](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetDistancesVector.md)
- [SetDivideGeometryUtility](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetDivideGeometryUtility.md)
- [pCreateModifiedShapeFunctions](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/pCreateModifiedShapeFunctions.md)
- [GetModifiedShapeFunctionsValues](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GetModifiedShapeFunctionsValues.md)
- [GetModifiedShapeFunctionsValuesOnEdge](../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GetModifiedShapeFunctionsValuesOnEdge.md)

# Imports

- `"variable_utils.h"`
- `"embedded_skin_utility.h"`
- `"modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"`
- `"modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"`
- `"modified_shape_functions/triangle_2d_3_modified_shape_functions.h"`
- `"modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"`
- `"utilities/divide_triangle_2d_3.h"`
- `"utilities/divide_tetrahedra_3d_4.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)