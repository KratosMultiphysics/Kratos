---
type: C++ Function
title: visit
resource: python/variable_indexing_python.h#L82-L94
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/controllers/output_controller/OutputController/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/controllers/output_controller/OutputController/Evaluate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/controllers/output_controller/OutputController/Update
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/controllers/output_controller/OutputController/GetCurrentControlValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/controllers/output_controller/OutputController/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetEntityName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/AddConnectivityData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/for
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/AddFields
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/GetOutputContainerList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/AddTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/PerformSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteMultigridLoop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/Solve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/MakeLhsTopology
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/Assemble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ApplyConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ApplyDirichletConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Clear
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ProjectGrid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Impl/WriteIntermediateState
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/ConnectivityIdsTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/CollectData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/FixityTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/CollectData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/StoreData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/flags_tensor_adaptor/FlagsTensorAdaptor/FlagsTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/flags_tensor_adaptor/FlagsTensorAdaptor/CollectData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/flags_tensor_adaptor/FlagsTensorAdaptor/StoreData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/gauss_point_variable_tensor_adaptor/GaussPointVariableTensorAdaptor/GaussPointVariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/gauss_point_variable_tensor_adaptor/GaussPointVariableTensorAdaptor/CollectData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/gauss_point_variable_tensor_adaptor/GaussPointVariableTensorAdaptor/Info
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/GeometriesTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/HistoricalVariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/CollectData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/StoreData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/Info
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/Info
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/VariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/CollectData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/StoreData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/Info
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_appended_data_element_wrapper/XmlAppendedDataElementWrapper/Get
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/xml_utilities/xml_data_element_wrapper/XmlDataElementWrapper/Get
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void visit(TClassType& ThisClass) const`

# Called by

- [Check](../../../functions/controllers/output_controller/OutputController/Check.md)
- [Evaluate](../../../functions/controllers/output_controller/OutputController/Evaluate.md)
- [Update](../../../functions/controllers/output_controller/OutputController/Update.md)
- [GetCurrentControlValue](../../../functions/controllers/output_controller/OutputController/GetCurrentControlValue.md)
- [PrintData](../../../functions/controllers/output_controller/OutputController/PrintData.md)
- [GetEntityName](../../../functions/input_output/vtu_output/GetEntityName.md)
- [AddConnectivityData](../../../functions/input_output/vtu_output/AddConnectivityData.md)
- [for](../../../functions/input_output/vtu_output/for.md)
- [AddFields](../../../functions/input_output/vtu_output/AddFields.md)
- [GetOutputContainerList](../../../functions/input_output/vtu_output/VtuOutput/GetOutputContainerList.md)
- [AddTensorAdaptor](../../../functions/input_output/vtu_output/VtuOutput/AddTensorAdaptor.md)
- [WriteUnstructuredGridData](../../../functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData.md)
- [WriteIntegrationPointData](../../../functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData.md)
- [PrintData](../../../functions/input_output/vtu_output/VtuOutput/PrintData.md)
- [InitializeSolutionStep](../../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/InitializeSolutionStep.md)
- [PerformSolutionStep](../../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/PerformSolutionStep.md)
- [SolveSubStep](../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep.md)
- [ExecuteMultigridLoop](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteMultigridLoop.md)
- [Solve](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/Solve.md)
- [MakeLhsTopology](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/MakeLhsTopology.md)
- [Assemble](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/Assemble.md)
- [ApplyConstraints](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ApplyConstraints.md)
- [ApplyDirichletConditions](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ApplyDirichletConditions.md)
- [Clear](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Clear.md)
- [ProjectGrid](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ProjectGrid.md)
- [WriteIntermediateState](../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Impl/WriteIntermediateState.md)
- [ConnectivityIdsTensorAdaptor](../../../functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/ConnectivityIdsTensorAdaptor.md)
- [Check](../../../functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/Check.md)
- [CollectData](../../../functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/CollectData.md)
- [FixityTensorAdaptor](../../../functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/FixityTensorAdaptor.md)
- [Check](../../../functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/Check.md)
- [CollectData](../../../functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/CollectData.md)
- [StoreData](../../../functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/StoreData.md)
- [FlagsTensorAdaptor](../../../functions/tensor_adaptors/flags_tensor_adaptor/FlagsTensorAdaptor/FlagsTensorAdaptor.md)
- [CollectData](../../../functions/tensor_adaptors/flags_tensor_adaptor/FlagsTensorAdaptor/CollectData.md)
- [StoreData](../../../functions/tensor_adaptors/flags_tensor_adaptor/FlagsTensorAdaptor/StoreData.md)
- [GaussPointVariableTensorAdaptor](../../../functions/tensor_adaptors/gauss_point_variable_tensor_adaptor/GaussPointVariableTensorAdaptor/GaussPointVariableTensorAdaptor.md)
- [CollectData](../../../functions/tensor_adaptors/gauss_point_variable_tensor_adaptor/GaussPointVariableTensorAdaptor/CollectData.md)
- [Info](../../../functions/tensor_adaptors/gauss_point_variable_tensor_adaptor/GaussPointVariableTensorAdaptor/Info.md)
- [GeometriesTensorAdaptor](../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/GeometriesTensorAdaptor.md)
- [CollectData](../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectData.md)
- [HistoricalVariableTensorAdaptor](../../../functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/HistoricalVariableTensorAdaptor.md)
- [Check](../../../functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/Check.md)
- [CollectData](../../../functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/CollectData.md)
- [StoreData](../../../functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/StoreData.md)
- [Info](../../../functions/tensor_adaptors/historical_variable_tensor_adaptor/HistoricalVariableTensorAdaptor/Info.md)
- [Info](../../../functions/tensor_adaptors/tensor_adaptor/TensorAdaptor<TDataType>/Info.md)
- [VariableTensorAdaptor](../../../functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/VariableTensorAdaptor.md)
- [Check](../../../functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/Check.md)
- [CollectData](../../../functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/CollectData.md)
- [StoreData](../../../functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/StoreData.md)
- [Info](../../../functions/tensor_adaptors/variable_tensor_adaptor/VariableTensorAdaptor/Info.md)
- [Get](../../../functions/utilities/xml_utilities/xml_appended_data_element_wrapper/XmlAppendedDataElementWrapper/Get.md)
- [Get](../../../functions/utilities/xml_utilities/xml_data_element_wrapper/XmlDataElementWrapper/Get.md)