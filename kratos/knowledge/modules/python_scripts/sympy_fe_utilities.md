---
type: Python Module
title: sympy_fe_utilities
resource: python_scripts/sympy_fe_utilities.py#L1-L825
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/re
    resolved_by: tree-sitter
    confidence: exact
  - target: external/sympy
    resolved_by: tree-sitter
    confidence: exact
  - target: external/sympy-codegen-rewriting
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [DefineMatrix](../../functions/python_scripts/sympy_fe_utilities/DefineMatrix.md)
- [DefineSymmetricMatrix](../../functions/python_scripts/sympy_fe_utilities/DefineSymmetricMatrix.md)
- [DefineSymmetricFourthOrderTensor](../../functions/python_scripts/sympy_fe_utilities/DefineSymmetricFourthOrderTensor.md)
- [DefineVector](../../functions/python_scripts/sympy_fe_utilities/DefineVector.md)
- [DefineShapeFunctions](../../functions/python_scripts/sympy_fe_utilities/DefineShapeFunctions.md)
- [StrainToMatrix](../../functions/python_scripts/sympy_fe_utilities/StrainToMatrix.md)
- [StrainToVoigt](../../functions/python_scripts/sympy_fe_utilities/StrainToVoigt.md)
- [VoigtToMatrix](../../functions/python_scripts/sympy_fe_utilities/VoigtToMatrix.md)
- [MatrixToVoigt](../../functions/python_scripts/sympy_fe_utilities/MatrixToVoigt.md)
- [ConvertTensorToVoigtMatrix](../../functions/python_scripts/sympy_fe_utilities/ConvertTensorToVoigtMatrix.md)
- [ConvertVoigtMatrixToTensor](../../functions/python_scripts/sympy_fe_utilities/ConvertVoigtMatrixToTensor.md)
- [DoubleContraction](../../functions/python_scripts/sympy_fe_utilities/DoubleContraction.md)
- [MatrixB](../../functions/python_scripts/sympy_fe_utilities/MatrixB.md)
- [grad_sym_voigtform](../../functions/python_scripts/sympy_fe_utilities/grad_sym_voigtform.md)
- [DfjDxi](../../functions/python_scripts/sympy_fe_utilities/DfjDxi.md)
- [DfiDxj](../../functions/python_scripts/sympy_fe_utilities/DfiDxj.md)
- [div](../../functions/python_scripts/sympy_fe_utilities/div.md)
- [SubstituteMatrixValue](../../functions/python_scripts/sympy_fe_utilities/SubstituteMatrixValue.md)
- [SubstituteScalarValue](../../functions/python_scripts/sympy_fe_utilities/SubstituteScalarValue.md)
- [SubstituteMatrixValueByVoigtSymbols](../../functions/python_scripts/sympy_fe_utilities/SubstituteMatrixValueByVoigtSymbols.md)
- [Compute_RHS](../../functions/python_scripts/sympy_fe_utilities/Compute_RHS.md)
- [Compute_LHS](../../functions/python_scripts/sympy_fe_utilities/Compute_LHS.md)
- [Compute_RHS_and_LHS](../../functions/python_scripts/sympy_fe_utilities/Compute_RHS_and_LHS.md)
- [_Indentation](../../functions/python_scripts/sympy_fe_utilities/_Indentation.md)
- [_CodeGen](../../functions/python_scripts/sympy_fe_utilities/_CodeGen.md)
- [_VariableDeclaration](../../functions/python_scripts/sympy_fe_utilities/_VariableDeclaration.md)
- [_Suffix](../../functions/python_scripts/sympy_fe_utilities/_Suffix.md)
- [_ReplaceIndices](../../functions/python_scripts/sympy_fe_utilities/_ReplaceIndices.md)
- [_GetVoigtToTensorConversionIndices](../../functions/python_scripts/sympy_fe_utilities/_GetVoigtToTensorConversionIndices.md)
- [_GetTensorToVoigtConversionIndices](../../functions/python_scripts/sympy_fe_utilities/_GetTensorToVoigtConversionIndices.md)
- [OutputScalar](../../functions/python_scripts/sympy_fe_utilities/OutputScalar.md)
- [OutputVector](../../functions/python_scripts/sympy_fe_utilities/OutputVector.md)
- [OutputMatrix](../../functions/python_scripts/sympy_fe_utilities/OutputMatrix.md)
- [OutputSymbolicVariable](../../functions/python_scripts/sympy_fe_utilities/OutputSymbolicVariable.md)
- [OutputSymbolicVariableDeclaration](../../functions/python_scripts/sympy_fe_utilities/OutputSymbolicVariableDeclaration.md)
- [_AuxiliaryOutputCollectionFactors](../../functions/python_scripts/sympy_fe_utilities/_AuxiliaryOutputCollectionFactors.md)
- [OutputMatrix_CollectingFactors](../../functions/python_scripts/sympy_fe_utilities/OutputMatrix_CollectingFactors.md)
- [OutputVector_CollectingFactors](../../functions/python_scripts/sympy_fe_utilities/OutputVector_CollectingFactors.md)
- [OutputScalar_CollectingFactors](../../functions/python_scripts/sympy_fe_utilities/OutputScalar_CollectingFactors.md)

# Imports

- `re`
- `sympy`
- `sympy.codegen.rewriting`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)