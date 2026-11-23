---
type: C++ Function
title: RetrieveGlobalIndexedPointersMap
resource: utilities/global_pointer_utilities.h#L95-L248
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/data_communicator/Recv
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Send
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Broadcast
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/global_pointer_utilities/ComputeGpMap
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/global_pointer_utilities/ExtractById
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/utilities/debug_utilities/CheckNonHistoricalVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/debug_utilities/CheckHistoricalVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CreatePointerCommunicator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_global_nodal_neighbours_for_entities_process/FindNodalNeighboursForEntitiesProcess<TContainerType>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/sensitivity_builder_scheme/SensitivityBuilderScheme/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointers
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::unordered_map< int, GlobalPointer<typename TContainerType::value_type> > RetrieveGlobalIndexedPointersMap( const TContainerType& rContainer, const std::vector<int>& rIdList, const DataCommunicator& rDataCommunicator )`

# Calls

- [Recv](../../../functions/includes/data_communicator/Recv.md)
- [Send](../../../functions/includes/data_communicator/Send.md)
- [Broadcast](../../../functions/includes/data_communicator/Broadcast.md)
- [ComputeGpMap](../../../functions/utilities/global_pointer_utilities/ComputeGpMap.md)
- [ExtractById](../../../functions/utilities/global_pointer_utilities/ExtractById.md)

# Called by

- [CheckNonHistoricalVariable](../../../functions/mpi/utilities/debug_utilities/CheckNonHistoricalVariable.md)
- [CheckHistoricalVariable](../../../functions/mpi/utilities/debug_utilities/CheckHistoricalVariable.md)
- [CreatePointerCommunicator](../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CreatePointerCommunicator.md)
- [Execute](../../../functions/processes/find_global_nodal_neighbours_for_entities_process/FindNodalNeighboursForEntitiesProcess<TContainerType>/Execute.md)
- [Initialize](../../../functions/solving_strategies/schemes/sensitivity_builder_scheme/SensitivityBuilderScheme/Initialize.md)
- [RetrieveGlobalIndexedPointers](../../../functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointers.md)