---
type: C++ Method
title: ComputeCommunicationScheduling
resource: utilities/communication_coloring_utilities.cpp#L56-L131
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/communication_coloring_utilities/MPIColoringUtilities/HasEdge
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrepareNonLocalCommunications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_sparse_graph/DistributedSparseGraph/Finalize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_vector_exporter/DistributedVectorExporter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_vector_importer/DistributedVectorImporter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/utilities/test_mpi_coloring_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<TContainerType>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_global_nodal_neighbours_for_entities_process/FindNodalNeighboursForEntitiesProcess<TContainerType>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/pointer_communicator/GlobalPointerCommunicator/ComputeCommunicationPlan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/pointer_map_communicator/ComputeCommunicationPlan
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<int> MPIColoringUtilities::ComputeCommunicationScheduling( const std::vector<int>& rLocalDestinationIds, const DataCommunicator& rComm )`

# Calls

- [HasEdge](../../../../functions/utilities/communication_coloring_utilities/MPIColoringUtilities/HasEdge.md)

# Called by

- [PrepareNonLocalCommunications](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrepareNonLocalCommunications.md)
- [Finalize](../../../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/Finalize.md)
- [DistributedVectorExporter](../../../../functions/containers/distributed_vector_exporter/DistributedVectorExporter.md)
- [DistributedVectorImporter](../../../../functions/containers/distributed_vector_importer/DistributedVectorImporter.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/utilities/test_mpi_coloring_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [Execute](../../../../functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<TContainerType>/Execute.md)
- [Execute](../../../../functions/processes/find_global_nodal_neighbours_for_entities_process/FindNodalNeighboursForEntitiesProcess<TContainerType>/Execute.md)
- [ComputeCommunicationPlan](../../../../functions/utilities/pointer_communicator/GlobalPointerCommunicator/ComputeCommunicationPlan.md)
- [ComputeCommunicationPlan](../../../../functions/utilities/pointer_map_communicator/ComputeCommunicationPlan.md)