---
type: C++ Method
title: AssemblyMortarOperators
resource: processes/simple_mortar_mapper_process.cpp#L380-L446
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/mortar_utilities/LengthCheck
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/GaussPointUnitNormal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/FastProjectDirection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mortar_classes/MortarOperator/CalculateMortarOperators
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/simple_mortar_mapper_process/PerformMortarOperations
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::AssemblyMortarOperators( const std::vector<array_1d<PointType,TDim>>& rGeometricalObjectsPointSlave, const GeometryType& rSlaveGeometry, const GeometryType& rMasterGeometry, const array_1d<double, 3>& rMasterNormal, MortarKinematicVariablesType& rThisKinematicVariables, MortarOperatorType& rThisMortarOperators, const IntegrationMethod& rThisIntegrationMethod, const BoundedMatrixType Ae )`

# Calls

- [LengthCheck](../../../../functions/utilities/mortar_utilities/LengthCheck.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)
- [GaussPointUnitNormal](../../../../functions/utilities/mortar_utilities/GaussPointUnitNormal.md)
- [FastProjectDirection](../../../../functions/utilities/geometrical_projection_utilities/FastProjectDirection.md)
- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [CalculateMortarOperators](../../../../functions/includes/mortar_classes/MortarOperator/CalculateMortarOperators.md)

# Called by

- [PerformMortarOperations](../../../../functions/processes/simple_mortar_mapper_process/PerformMortarOperations.md)