---
type: C++ Function
title: CalculateRotationOperatorPureShapeSensitivities
resource: utilities/coordinate_transformation_utilities.h#L225-L254
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/SolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/VectorNormDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/UnitVectorDerivative
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void CalculateRotationOperatorPureShapeSensitivities( TLocalMatrixType& rRotationMatrixShapeDerivative, const std::size_t DerivativeNodeIndex, const std::size_t DerivativeDirectionIndex, const GeometryType::PointType& rThisPoint) const`

# Calls

- [SolutionStepsDataHas](../../../functions/includes/node/SolutionStepsDataHas.md)
- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)
- [VectorNormDerivative](../../../functions/utilities/coordinate_transformation_utilities/VectorNormDerivative.md)
- [UnitVectorDerivative](../../../functions/utilities/coordinate_transformation_utilities/UnitVectorDerivative.md)