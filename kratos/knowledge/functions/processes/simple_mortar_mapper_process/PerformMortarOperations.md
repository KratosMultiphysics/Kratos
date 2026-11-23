---
type: C++ Function
title: PerformMortarOperations
resource: processes/simple_mortar_mapper_process.h#L519-L673
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ComputeResidualMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/LumpMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_frobenius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssembleRHSAndLHS
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PerformMortarOperations( MatrixType& rA, std::vector<VectorType>& rb, IntMap& rInverseConectivityDatabase, typename TClassType::Pointer pIndexesPairs, GeometricalObject& rGeometricalObject, ExactMortarIntegrationUtilityType& rIntegrationUtility, MortarKinematicVariablesType& rThisKineticVariables, MortarOperatorType& rThisMortarOperators, const IndexType Iteration )`

# Calls

- [GetDouble](../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [AssemblyMortarOperators](../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators.md)
- [ComputeResidualMatrix](../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ComputeResidualMatrix.md)
- [LumpMatrix](../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/LumpMatrix.md)
- [norm_frobenius](../../../functions/includes/amatrix_interface/norm_frobenius.md)
- [AssembleRHSAndLHS](../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssembleRHSAndLHS.md)