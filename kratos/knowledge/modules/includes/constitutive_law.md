---
type: C++ Module
title: constitutive_law
resource: includes/constitutive_law.h#L1-L1605
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-serializer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-variables-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-node-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-properties-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-geometry-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-math-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-process-info-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-ublas-interface-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-parameters-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-data-value-container-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-flags-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-initial-state-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [Features](../../classes/includes/constitutive_law/Features.md)
- [Features](../../functions/includes/constitutive_law/Features/Features.md)
- [SetOptions](../../functions/includes/constitutive_law/Features/SetOptions.md)
- [SetStrainSize](../../functions/includes/constitutive_law/Features/SetStrainSize.md)
- [SetSpaceDimension](../../functions/includes/constitutive_law/Features/SetSpaceDimension.md)
- [SetStrainMeasure](../../functions/includes/constitutive_law/Features/SetStrainMeasure.md)
- [SetStrainMeasures](../../functions/includes/constitutive_law/Features/SetStrainMeasures.md)
- [Parameters](../../classes/includes/constitutive_law/Parameters.md)
- [Parameters](../../functions/includes/constitutive_law/Parameters/Parameters.md)
- [Parameters](../../functions/includes/constitutive_law/Parameters/Parameters-2.md)
- [Parameters](../../functions/includes/constitutive_law/Parameters/Parameters-3.md)
- [CheckAllParameters](../../functions/includes/constitutive_law/Parameters/CheckAllParameters.md)
- [CheckShapeFunctions](../../functions/includes/constitutive_law/Parameters/CheckShapeFunctions.md)
- [CheckInfoMaterialGeometry](../../functions/includes/constitutive_law/Parameters/CheckInfoMaterialGeometry.md)
- [CheckMechanicalVariables](../../functions/includes/constitutive_law/Parameters/CheckMechanicalVariables.md)
- [Set](../../functions/includes/constitutive_law/Parameters/Set.md)
- [Reset](../../functions/includes/constitutive_law/Parameters/Reset.md)
- [SetOptions](../../functions/includes/constitutive_law/Parameters/SetOptions.md)
- [SetDeterminantF](../../functions/includes/constitutive_law/Parameters/SetDeterminantF.md)
- [SetShapeFunctionsValues](../../functions/includes/constitutive_law/Parameters/SetShapeFunctionsValues.md)
- [SetShapeFunctionsDerivatives](../../functions/includes/constitutive_law/Parameters/SetShapeFunctionsDerivatives.md)
- [SetDeformationGradientF](../../functions/includes/constitutive_law/Parameters/SetDeformationGradientF.md)
- [SetStrainVector](../../functions/includes/constitutive_law/Parameters/SetStrainVector.md)
- [SetStressVector](../../functions/includes/constitutive_law/Parameters/SetStressVector.md)
- [SetConstitutiveMatrix](../../functions/includes/constitutive_law/Parameters/SetConstitutiveMatrix.md)
- [SetProcessInfo](../../functions/includes/constitutive_law/Parameters/SetProcessInfo.md)
- [SetMaterialProperties](../../functions/includes/constitutive_law/Parameters/SetMaterialProperties.md)
- [SetElementGeometry](../../functions/includes/constitutive_law/Parameters/SetElementGeometry.md)
- [IsSetDeterminantF](../../functions/includes/constitutive_law/Parameters/IsSetDeterminantF.md)
- [IsSetShapeFunctionsValues](../../functions/includes/constitutive_law/Parameters/IsSetShapeFunctionsValues.md)
- [IsSetShapeFunctionsDerivatives](../../functions/includes/constitutive_law/Parameters/IsSetShapeFunctionsDerivatives.md)
- [IsSetDeformationGradientF](../../functions/includes/constitutive_law/Parameters/IsSetDeformationGradientF.md)
- [IsSetStrainVector](../../functions/includes/constitutive_law/Parameters/IsSetStrainVector.md)
- [IsSetStressVector](../../functions/includes/constitutive_law/Parameters/IsSetStressVector.md)
- [IsSetConstitutiveMatrix](../../functions/includes/constitutive_law/Parameters/IsSetConstitutiveMatrix.md)
- [IsSetProcessInfo](../../functions/includes/constitutive_law/Parameters/IsSetProcessInfo.md)
- [IsSetMaterialProperties](../../functions/includes/constitutive_law/Parameters/IsSetMaterialProperties.md)
- [IsSetElementGeometry](../../functions/includes/constitutive_law/Parameters/IsSetElementGeometry.md)
- [SetInitialState](../../functions/includes/constitutive_law/SetInitialState.md)
- [pGetInitialState](../../functions/includes/constitutive_law/pGetInitialState.md)
- [HasInitialState](../../functions/includes/constitutive_law/HasInitialState.md)
- [AddInitialStressVectorContribution](../../functions/includes/constitutive_law/AddInitialStressVectorContribution.md)
- [AddInitialStrainVectorContribution](../../functions/includes/constitutive_law/AddInitialStrainVectorContribution.md)
- [AddInitialDeformationGradientMatrixContribution](../../functions/includes/constitutive_law/AddInitialDeformationGradientMatrixContribution.md)
- [RequiresInitializeMaterialResponse](../../functions/includes/constitutive_law/RequiresInitializeMaterialResponse.md)
- [RequiresFinalizeMaterialResponse](../../functions/includes/constitutive_law/RequiresFinalizeMaterialResponse.md)
- [HasSameType](../../functions/includes/constitutive_law/HasSameType.md)
- [HasSameType](../../functions/includes/constitutive_law/HasSameType-2.md)
- [Info](../../functions/includes/constitutive_law/Info.md)
- [PrintInfo](../../functions/includes/constitutive_law/PrintInfo.md)
- [PrintData](../../functions/includes/constitutive_law/PrintData.md)
- [save](../../functions/includes/constitutive_law/save.md)
- [load](../../functions/includes/constitutive_law/load.md)

# Imports

- `"includes/define.h"`
- `"includes/serializer.h"`
- `"includes/variables.h"`
- `"includes/node.h"`
- `"includes/properties.h"`
- `"geometries/geometry.h"`
- `"utilities/math_utils.h"`
- `"includes/process_info.h"`
- `"includes/ublas_interface.h"`
- `"includes/kratos_parameters.h"`
- `"containers/data_value_container.h"`
- `"containers/flags.h"`
- `"includes/initial_state.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)