//     __ __      __  __________  ___                ___            __  _           
//    / //_/___ _/ / / /  _/ __ \/   |  ____  ____  / (_)________ _/ /_(_)___  ____ 
//   / ,< / __ `/ /_/ // // /_/ / /| | / __ \/ __ \/ / / ___/ __ `/ __/ / __ \/ __ \
//  / /| / /_/ / __  // // ____/ ___ |/ /_/ / /_/ / / / /__/ /_/ / /_/ / / /_/ / / /
// /_/ |_\__,_/_/ /_/___/_/   /_/  |_/ .___/ .___/_/_/\___/\__,_/\__/_/\____/_/ /_/ 
//                                  /_/   /_/                                       
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <algorithm>
#include <unordered_set>

// External includes

// Project includes
#include "processes/graph_coloring_process.h"
#include "custom_processes/kahip_divide_heterogeneous_input_process.h"
#include "custom_utilities/kahip_csr_converter.h"

// Utilities shared with Metis (graph / domains bookkeeping is identical)
// We duplicate the few helpers we need rather than pulling in MetisApplication
// as a dependency.

namespace Kratos {

namespace {

/// Build per-partition domain graph from connectivity and partition arrays.
/// Marks (node_rank, element_rank) as connected when they differ.
void CalculateDomainsGraph(
    IO::GraphType& rDomainsGraph,
    std::size_t NumberOfEntities,
    const IO::ConnectivitiesContainerType& rConnectivities,
    const std::vector<kahip_idx>& rNodePartition,
    const std::vector<kahip_idx>& rEntityPartition)
{
    for (std::size_t i = 0; i < NumberOfEntities; ++i) {
        for (std::size_t node_id : rConnectivities[i]) {
            const std::size_t node_rank   = rNodePartition[node_id - 1];
            const std::size_t entity_rank = rEntityPartition[i];
            if (node_rank != entity_rank) {
                rDomainsGraph(node_rank, entity_rank) = 1;
                rDomainsGraph(entity_rank, node_rank) = 1;
            }
        }
    }
}

/// Collect the set of partitions each node must appear in.
/// A node must be visible to every partition that has an element using it.
void DividingNodes(
    IO::PartitionIndicesContainerType& rNodesAllPartitions,
    const IO::ConnectivitiesContainerType& rGeometryConnectivities,
    const IO::ConnectivitiesContainerType& rElementConnectivities,
    const IO::ConnectivitiesContainerType& rConditionConnectivities,
    const IO::ConnectivitiesContainerType& rConstraintConnectivities,
    const std::vector<kahip_idx>& rNodePartition,
    const std::vector<kahip_idx>& rGeometryPartition,
    const std::vector<kahip_idx>& rElementPartition,
    const std::vector<kahip_idx>& rConditionPartition,
    const std::vector<kahip_idx>& rConstraintPartition)
{
    const std::size_t num_nodes = rNodePartition.size();
    rNodesAllPartitions.resize(num_nodes);

    auto process_connectivity = [&](
        const IO::ConnectivitiesContainerType& connectivities,
        const std::vector<kahip_idx>& entity_partition)
    {
        for (std::size_t i = 0; i < entity_partition.size(); ++i) {
            const int ep = entity_partition[i];
            for (std::size_t node_id : connectivities[i]) {
                const int np = rNodePartition[node_id - 1];
                if (ep != np) {
                    rNodesAllPartitions[node_id - 1].push_back(ep);
                }
            }
        }
    };

    process_connectivity(rGeometryConnectivities,  rGeometryPartition);
    process_connectivity(rElementConnectivities,    rElementPartition);
    process_connectivity(rConditionConnectivities,  rConditionPartition);
    process_connectivity(rConstraintConnectivities, rConstraintPartition);

    // Add the node's own partition
    for (std::size_t i = 0; i < num_nodes; ++i) {
        rNodesAllPartitions[i].push_back(rNodePartition[i]);
    }
}

void DividingEntities(
    IO::PartitionIndicesContainerType& rEntitiesAllPartitions,
    const std::vector<kahip_idx>& rEntityPartition)
{
    const std::size_t n = rEntityPartition.size();
    rEntitiesAllPartitions.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        rEntitiesAllPartitions[i].push_back(rEntityPartition[i]);
    }
}

} // anonymous namespace

// ─── Constructors ────────────────────────────────────────────────────────────

KaHIPDivideHeterogeneousInputProcess::KaHIPDivideHeterogeneousInputProcess(
    IO& rIO,
    SizeType NumberOfPartitions,
    Parameters rSettings,
    bool SynchronizeConditions)
    : mrIO(rIO),
      mNumberOfPartitions(NumberOfPartitions),
      mSynchronizeConditions(SynchronizeConditions),
      mVerbosity(0),
      mNumNodes(0),
      mPartitioner([&]() {
          // Extract top-level settings that belong to the process, pass the rest
          // to the partitioner so it can call ValidateAndAssignDefaults cleanly.
          Parameters partitioner_settings(R"({})");
          const std::vector<std::string> partitioner_keys = {
              "preconfiguration", "imbalance", "seed", "suppress_output", "num_trials"};
          for (const auto& key : partitioner_keys) {
              if (rSettings.Has(key)) {
                  partitioner_settings.AddValue(key, rSettings[key]);
              }
          }
          if (rSettings.Has("verbosity")) {
              mVerbosity = rSettings["verbosity"].GetInt();
          }
          if (rSettings.Has("synchronize_conditions")) {
              mSynchronizeConditions = rSettings["synchronize_conditions"].GetBool();
          }
          return KaHIPPartitioner(partitioner_settings);
      }())
{
}

KaHIPDivideHeterogeneousInputProcess::KaHIPDivideHeterogeneousInputProcess(
    IO& rIO,
    SizeType NumberOfPartitions,
    int Dimension,
    int Verbosity,
    bool SynchronizeConditions)
    : mrIO(rIO),
      mNumberOfPartitions(NumberOfPartitions),
      mSynchronizeConditions(SynchronizeConditions),
      mVerbosity(Verbosity),
      mNumNodes(0),
      mPartitioner(KaHIPPartitioner::GetDefaultParameters())
{
}

// ─── Execute ─────────────────────────────────────────────────────────────────

void KaHIPDivideHeterogeneousInputProcess::Execute()
{
    PartitioningInfo part_info;
    ExecutePartitioning(part_info);

    mrIO.DivideInputToPartitions(mNumberOfPartitions, part_info);
}

// ─── GetNodesPartitions ───────────────────────────────────────────────────────

void KaHIPDivideHeterogeneousInputProcess::GetNodesPartitions(
    std::vector<idxtype>& rNodePartition,
    SizeType& rNumNodes)
{
    // Read the nodal connectivity graph from the IO object.
    // ReadNodalGraph() returns a ConnectivitiesContainerType where entry i holds the
    // 1-indexed IDs of all neighbours of node i+1.
    IO::ConnectivitiesContainerType kratos_node_connectivities;
    rNumNodes = mrIO.ReadNodalGraph(kratos_node_connectivities);

    const SizeType num_nodes_in_mesh = mrIO.ReadNodesNumber();
    KRATOS_ERROR_IF(rNumNodes != num_nodes_in_mesh)
        << "KaHIPDivideHeterogeneousInputProcess: inconsistent node count — "
        << "connected nodes = " << rNumNodes
        << ", mesh nodes = " << num_nodes_in_mesh << std::endl;

    // Convert to KaHIP CSR format (0-indexed, both directions of each edge present)
    std::vector<kahip_idx> xadj, adjncy;
    KaHIPCSRConverter::ConvertToCSRFormat(kratos_node_connectivities, xadj, adjncy);

    KRATOS_INFO_IF("KaHIPDivideHeterogeneousInputProcess", mVerbosity > 0)
        << "Partitioning nodal graph: " << rNumNodes << " nodes, "
        << adjncy.size() / 2 << " edges, "
        << mNumberOfPartitions << " partitions." << std::endl;

    // Empty weight vectors → kaffpa uses uniform weights
    std::vector<int>       vwgt;
    std::vector<kahip_idx> adjcwgt;

    int n = static_cast<int>(rNumNodes);
    const std::vector<int> partition = mPartitioner.PartitionGraph(
        n, xadj, adjncy, vwgt, adjcwgt, static_cast<int>(mNumberOfPartitions));

    // Cache the node count so PartitionElementsSynchronous and
    // PartitionConstraintsSynchronous can size mNodeConnectivities correctly.
    mNumNodes = static_cast<int>(rNumNodes);

    rNodePartition.resize(rNumNodes);
    for (SizeType i = 0; i < rNumNodes; ++i) {
        rNodePartition[i] = static_cast<idxtype>(partition[i]);
    }

    PrintDebugData("Node Partition", rNodePartition);
}

// ─── ExecutePartitioning ─────────────────────────────────────────────────────

void KaHIPDivideHeterogeneousInputProcess::ExecutePartitioning(PartitioningInfo& rPartitioningInfo)
{
    SizeType number_of_nodes;
    std::vector<idxtype> node_partition;
    GetNodesPartitions(node_partition, number_of_nodes);

    // ── Geometries ──────────────────────────────────────────────────────────
    IO::ConnectivitiesContainerType geometry_connectivities;
    const SizeType number_of_geometries = mrIO.ReadGeometriesConnectivities(geometry_connectivities);
    KRATOS_ERROR_IF(number_of_geometries != geometry_connectivities.size())
        << "KaHIPDivideHeterogeneousInputProcess: read " << number_of_geometries
        << " geometries but connectivity list has " << geometry_connectivities.size()
        << " entries (non-correlative IDs?)" << std::endl;

    std::vector<idxtype> geometry_partition;
    if (mSynchronizeConditions)
        PartitionGeometriesSynchronous(node_partition, geometry_connectivities, geometry_partition);
    else
        PartitionMesh(node_partition, geometry_connectivities, geometry_partition);

    // ── Elements ────────────────────────────────────────────────────────────
    IO::ConnectivitiesContainerType element_connectivities;
    const SizeType number_of_elements = mrIO.ReadElementsConnectivities(element_connectivities);
    KRATOS_ERROR_IF(number_of_elements != element_connectivities.size())
        << "KaHIPDivideHeterogeneousInputProcess: read " << number_of_elements
        << " elements but connectivity list has " << element_connectivities.size()
        << " entries (non-correlative IDs?)" << std::endl;

    std::vector<idxtype> element_partition;
    if (mSynchronizeConditions)
        PartitionElementsSynchronous(node_partition, element_connectivities, element_partition);
    else
        PartitionMesh(node_partition, element_connectivities, element_partition);

    // ── Conditions ──────────────────────────────────────────────────────────
    IO::ConnectivitiesContainerType condition_connectivities;
    const SizeType number_of_conditions = mrIO.ReadConditionsConnectivities(condition_connectivities);
    KRATOS_ERROR_IF(number_of_conditions != condition_connectivities.size())
        << "KaHIPDivideHeterogeneousInputProcess: read " << number_of_conditions
        << " conditions but connectivity list has " << condition_connectivities.size()
        << " entries (non-correlative IDs?)" << std::endl;

    std::vector<idxtype> condition_partition;
    if (mSynchronizeConditions)
        PartitionConditionsSynchronous(
            node_partition, element_partition,
            condition_connectivities, element_connectivities,
            condition_partition);
    else
        PartitionMesh(node_partition, condition_connectivities, condition_partition);

    // ── Master-Slave Constraints ─────────────────────────────────────────────
    IO::ConnectivitiesContainerType constraint_connectivities;
    const SizeType number_of_constraints = mrIO.ReadConstraintsConnectivities(constraint_connectivities);
    KRATOS_ERROR_IF(number_of_constraints != constraint_connectivities.size())
        << "KaHIPDivideHeterogeneousInputProcess: read " << number_of_constraints
        << " constraints but connectivity list has " << constraint_connectivities.size()
        << " entries (non-correlative IDs?)" << std::endl;

    std::vector<idxtype> constraint_partition;
    if (mSynchronizeConditions)
        PartitionConstraintsSynchronous(
            node_partition, constraint_connectivities, constraint_partition);
    else
        PartitionMesh(node_partition, constraint_connectivities, constraint_partition);

    // ── Relocate hanging nodes ───────────────────────────────────────────────
    // Nodes that belong to a partition where no local element/condition/geometry
    // references them are re-assigned to avoid isolated nodes.
    RedistributeHangingNodes(
        node_partition,
        geometry_partition,  geometry_connectivities,
        element_partition,   element_connectivities,
        condition_partition, condition_connectivities,
        constraint_partition, constraint_connectivities);

    // ── Graph colouring (for MPI communication schedule) ─────────────────────
    GraphType domain_graph = zero_matrix<int>(mNumberOfPartitions);
    CalculateDomainsGraph(domain_graph, number_of_geometries,  geometry_connectivities,  node_partition, geometry_partition);
    CalculateDomainsGraph(domain_graph, number_of_elements,    element_connectivities,   node_partition, element_partition);
    CalculateDomainsGraph(domain_graph, number_of_conditions,  condition_connectivities, node_partition, condition_partition);
    CalculateDomainsGraph(domain_graph, number_of_constraints, constraint_connectivities, node_partition, constraint_partition);

    int number_of_colors;
    GraphColoringProcess(mNumberOfPartitions, domain_graph, rPartitioningInfo.Graph, number_of_colors).Execute();

    KRATOS_INFO_IF("KaHIPDivideHeterogeneousInputProcess", mVerbosity > 0)
        << "Number of colors: " << number_of_colors << std::endl;

    KRATOS_INFO_IF("KaHIPDivideHeterogeneousInputProcess", mVerbosity > 2)
        << "Graph: " << rPartitioningInfo.Graph << std::endl;

    // ── Build per-partition node/entity lists ────────────────────────────────
    DividingNodes(
        rPartitioningInfo.NodesAllPartitions,
        geometry_connectivities, element_connectivities,
        condition_connectivities, constraint_connectivities,
        node_partition, geometry_partition, element_partition,
        condition_partition, constraint_partition);

    DividingEntities(rPartitioningInfo.GeometriesAllPartitions,  geometry_partition);
    DividingEntities(rPartitioningInfo.ElementsAllPartitions,    element_partition);
    DividingEntities(rPartitioningInfo.ConditionsAllPartitions,  condition_partition);
    DividingEntities(rPartitioningInfo.ConstraintsAllPartitions, constraint_partition);

    if (mVerbosity > 1) {
        auto& r_nodes_all = rPartitioningInfo.NodesAllPartitions;
        KRATOS_INFO("KaHIPDivideHeterogeneousInputProcess")
            << "Final list of nodes known by each partition:" << std::endl;
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            std::cout << "Node #" << (i + 1) << " -> ";
            for (std::size_t p : r_nodes_all[i]) std::cout << p << ",";
            std::cout << std::endl;
        }
    }

    rPartitioningInfo.NodesPartitions.assign(node_partition.begin(), node_partition.end());
    rPartitioningInfo.GeometriesPartitions.assign(geometry_partition.begin(), geometry_partition.end());
    rPartitioningInfo.ElementsPartitions.assign(element_partition.begin(), element_partition.end());
    rPartitioningInfo.ConditionsPartitions.assign(condition_partition.begin(), condition_partition.end());
    rPartitioningInfo.ConstraintsPartitions.assign(constraint_partition.begin(), constraint_partition.end());
}

// ─── PartitionMesh ───────────────────────────────────────────────────────────

void KaHIPDivideHeterogeneousInputProcess::PartitionMesh(
    const std::vector<idxtype>& rNodePartition,
    const IO::ConnectivitiesContainerType& rConnectivities,
    std::vector<idxtype>& rEntityPartition)
{
    const SizeType num_entities = rConnectivities.size();
    std::vector<int> partition_weights(mNumberOfPartitions, 0);

    rEntityPartition.resize(num_entities, -1);

    // First pass: entities whose nodes all belong to the same partition go there.
    auto it_entity = rConnectivities.begin();
    for (auto it_part = rEntityPartition.begin(); it_part != rEntityPartition.end(); ++it_part, ++it_entity) {
        const int my_partition = rNodePartition[(*it_entity)[0] - 1];
        SizeType same_partition_nodes = 1;

        for (auto it_node = it_entity->begin() + 1; it_node != it_entity->end(); ++it_node) {
            if (rNodePartition[*it_node - 1] == my_partition)
                ++same_partition_nodes;
            else
                break;
        }

        if (same_partition_nodes == it_entity->size()) {
            *it_part = my_partition;
            partition_weights[my_partition]++;
        }
    }

    // Second pass: boundary entities are assigned to the partition that owns most of their nodes.
    const int max_weight = static_cast<int>(1.03 * num_entities / mNumberOfPartitions);
    it_entity = rConnectivities.begin();
    for (auto it_part = rEntityPartition.begin(); it_part != rEntityPartition.end(); ++it_part, ++it_entity) {
        if (*it_part != -1) continue;

        const SizeType nodes_in_entity = it_entity->size();
        std::vector<int> neighbour_partitions(nodes_in_entity, -1);
        std::vector<int> neighbour_weights(nodes_in_entity, 0);
        SizeType found_neighbours = 0;

        for (std::size_t node_id : *it_entity) {
            const int node_partition = rNodePartition[node_id - 1];
            SizeType i = 0;
            for (; i < found_neighbours; ++i) {
                if (node_partition == neighbour_partitions[i]) {
                    neighbour_weights[i]++;
                    break;
                }
            }
            if (i == found_neighbours) {
                neighbour_weights[found_neighbours]   = 1;
                neighbour_partitions[found_neighbours] = node_partition;
                ++found_neighbours;
            }
        }

        const int majority = neighbour_partitions[FindMax(found_neighbours, neighbour_weights)];
        if (partition_weights[majority] < max_weight) {
            *it_part = majority;
            partition_weights[majority]++;
        } else {
            // Majority partition is full; try the others
            bool assigned = false;
            for (SizeType i = 0; i < found_neighbours; ++i) {
                const int dest = neighbour_partitions[i];
                if (partition_weights[dest] < max_weight) {
                    *it_part = dest;
                    partition_weights[dest]++;
                    assigned = true;
                    break;
                }
            }
            if (!assigned) {
                // All candidate partitions are full — force majority
                *it_part = majority;
                partition_weights[majority]++;
            }
        }
    }

    PrintDebugData("Mesh Partition", rEntityPartition);
}

// ─── PartitionGeometriesSynchronous ──────────────────────────────────────────

void KaHIPDivideHeterogeneousInputProcess::PartitionGeometriesSynchronous(
    const std::vector<idxtype>& rNodePartition,
    const IO::ConnectivitiesContainerType& rGeomConnectivities,
    std::vector<idxtype>& rGeomPartition)
{
    const SizeType num_geom = rGeomConnectivities.size();
    rGeomPartition.resize(num_geom, -1);

    // First pass: all nodes in same partition
    auto it_geom = rGeomConnectivities.begin();
    for (auto it_part = rGeomPartition.begin(); it_part != rGeomPartition.end(); ++it_part, ++it_geom) {
        const int my_partition = rNodePartition[(*it_geom)[0] - 1];
        SizeType same = 1;
        for (auto it_node = it_geom->begin() + 1; it_node != it_geom->end(); ++it_node) {
            if (rNodePartition[*it_node - 1] == my_partition) ++same; else break;
        }
        if (same == it_geom->size()) *it_part = my_partition;
    }

    // Second pass: boundary geometries → majority partition (no balance constraint for geometries)
    it_geom = rGeomConnectivities.begin();
    for (auto it_part = rGeomPartition.begin(); it_part != rGeomPartition.end(); ++it_part, ++it_geom) {
        if (*it_part != -1) continue;

        const SizeType nodes_in_geom = it_geom->size();
        std::vector<int> neighbour_partitions(nodes_in_geom, -1);
        std::vector<int> neighbour_weights(nodes_in_geom, 0);
        SizeType found = 0;

        for (std::size_t node_id : *it_geom) {
            const int np = rNodePartition[node_id - 1];
            SizeType i = 0;
            for (; i < found; ++i) {
                if (np == neighbour_partitions[i]) { neighbour_weights[i]++; break; }
            }
            if (i == found) {
                neighbour_weights[found]   = 1;
                neighbour_partitions[found] = np;
                ++found;
            }
        }
        *it_part = neighbour_partitions[FindMax(found, neighbour_weights)];
    }

    PrintDebugData("Geometry Partition", rGeomPartition);
}

// ─── PartitionElementsSynchronous ────────────────────────────────────────────

void KaHIPDivideHeterogeneousInputProcess::PartitionElementsSynchronous(
    const std::vector<idxtype>& rNodePartition,
    const IO::ConnectivitiesContainerType& rElemConnectivities,
    std::vector<idxtype>& rElemPartition)
{
    const SizeType num_elements = rElemConnectivities.size();

    // Build node→element connectivity for later condition synchronisation.
    mNodeConnectivities.assign(mNumNodes, std::unordered_set<std::size_t>());
    for (std::size_t i = 0; i < num_elements; ++i) {
        for (std::size_t node_id : rElemConnectivities[i]) {
            mNodeConnectivities[node_id - 1].insert(i);
        }
    }

    rElemPartition.resize(num_elements, -1);

    // First pass: elements fully interior to one partition
    auto it_elem = rElemConnectivities.begin();
    for (auto it_part = rElemPartition.begin(); it_part != rElemPartition.end(); ++it_part, ++it_elem) {
        const int my_partition = rNodePartition[(*it_elem)[0] - 1];
        SizeType same = 1;
        for (auto it_node = it_elem->begin() + 1; it_node != it_elem->end(); ++it_node) {
            if (rNodePartition[*it_node - 1] == my_partition) ++same; else break;
        }
        if (same == it_elem->size()) *it_part = my_partition;
    }

    // Second pass: boundary elements → majority partition
    it_elem = rElemConnectivities.begin();
    for (auto it_part = rElemPartition.begin(); it_part != rElemPartition.end(); ++it_part, ++it_elem) {
        if (*it_part != -1) continue;

        const SizeType nodes_in_elem = it_elem->size();
        std::vector<int> neighbour_partitions(nodes_in_elem, -1);
        std::vector<int> neighbour_weights(nodes_in_elem, 0);
        SizeType found = 0;

        for (std::size_t node_id : *it_elem) {
            const int np = rNodePartition[node_id - 1];
            SizeType i = 0;
            for (; i < found; ++i) {
                if (np == neighbour_partitions[i]) { neighbour_weights[i]++; break; }
            }
            if (i == found) {
                neighbour_weights[found]   = 1;
                neighbour_partitions[found] = np;
                ++found;
            }
        }
        *it_part = neighbour_partitions[FindMax(found, neighbour_weights)];
    }

    PrintDebugData("Element Partition", rElemPartition);
}

// ─── PartitionConditionsSynchronous ──────────────────────────────────────────

void KaHIPDivideHeterogeneousInputProcess::PartitionConditionsSynchronous(
    const std::vector<idxtype>& rNodePartition,
    const std::vector<idxtype>& rElemPartition,
    const IO::ConnectivitiesContainerType& rCondConnectivities,
    const IO::ConnectivitiesContainerType& rElemConnectivities,
    std::vector<idxtype>& rCondPartition)
{
    const SizeType num_elements  = rElemConnectivities.size();
    const SizeType num_conditions = rCondConnectivities.size();

    rCondPartition.resize(num_conditions, -1);

    // Pre-sort element connectivities for O(n log n) subset checks
    IO::ConnectivitiesContainerType elements_sorted(rElemConnectivities);
    for (SizeType i = 0; i < num_elements; ++i) {
        std::sort(elements_sorted[i].begin(), elements_sorted[i].end());
    }

    auto it_cond = rCondConnectivities.begin();
    for (auto it_part = rCondPartition.begin(); it_part != rCondPartition.end(); ++it_part, ++it_cond) {
        const SizeType nodes_in_cond = it_cond->size();
        std::vector<int> neighbour_partitions(nodes_in_cond, -1);
        std::vector<int> neighbour_weights(nodes_in_cond, 0);
        SizeType found = 0;

        for (std::size_t node_id : *it_cond) {
            const int np = rNodePartition[node_id - 1];
            SizeType i = 0;
            for (; i < found; ++i) {
                if (np == neighbour_partitions[i]) { neighbour_weights[i]++; break; }
            }
            if (i == found) {
                neighbour_weights[found]   = 1;
                neighbour_partitions[found] = np;
                ++found;
            }
        }
        *it_part = neighbour_partitions[FindMax(found, neighbour_weights)];

        // Synchronise: if this condition is a face of an element, assign it to that element's partition
        IO::ConnectivitiesContainerType::value_type cond_sorted(*it_cond);
        std::sort(cond_sorted.begin(), cond_sorted.end());

        for (std::size_t node_id : *it_cond) {
            for (std::size_t elem_idx : mNodeConnectivities[node_id - 1]) {
                if (std::includes(
                        elements_sorted[elem_idx].begin(), elements_sorted[elem_idx].end(),
                        cond_sorted.begin(), cond_sorted.end()))
                {
                    *it_part = rElemPartition[elem_idx];
                    break;
                }
            }
        }
    }

    PrintDebugData("Condition Partition", rCondPartition);
}

// ─── PartitionConstraintsSynchronous ─────────────────────────────────────────

void KaHIPDivideHeterogeneousInputProcess::PartitionConstraintsSynchronous(
    const std::vector<idxtype>& rNodePartition,
    const IO::ConnectivitiesContainerType& rConstraintConnectivities,
    std::vector<idxtype>& rConstraintPartition)
{
    const SizeType num_constraints = rConstraintConnectivities.size();

    mNodeConnectivities.assign(mNumNodes, std::unordered_set<std::size_t>());
    for (std::size_t i = 0; i < num_constraints; ++i) {
        for (std::size_t node_id : rConstraintConnectivities[i]) {
            mNodeConnectivities[node_id - 1].insert(i);
        }
    }

    rConstraintPartition.resize(num_constraints, -1);

    // First pass: constraints fully interior to one partition
    auto it_const = rConstraintConnectivities.begin();
    for (auto it_part = rConstraintPartition.begin(); it_part != rConstraintPartition.end(); ++it_part, ++it_const) {
        const int my_partition = rNodePartition[(*it_const)[0] - 1];
        SizeType same = 1;
        for (auto it_node = it_const->begin() + 1; it_node != it_const->end(); ++it_node) {
            if (rNodePartition[*it_node - 1] == my_partition) ++same; else break;
        }
        if (same == it_const->size()) *it_part = my_partition;
    }

    // Second pass: boundary constraints → majority partition
    it_const = rConstraintConnectivities.begin();
    for (auto it_part = rConstraintPartition.begin(); it_part != rConstraintPartition.end(); ++it_part, ++it_const) {
        if (*it_part != -1) continue;

        const SizeType nodes_in_const = it_const->size();
        std::vector<int> neighbour_partitions(nodes_in_const, -1);
        std::vector<int> neighbour_weights(nodes_in_const, 0);
        SizeType found = 0;

        for (std::size_t node_id : *it_const) {
            const int np = rNodePartition[node_id - 1];
            SizeType i = 0;
            for (; i < found; ++i) {
                if (np == neighbour_partitions[i]) { neighbour_weights[i]++; break; }
            }
            if (i == found) {
                neighbour_weights[found]   = 1;
                neighbour_partitions[found] = np;
                ++found;
            }
        }
        *it_part = neighbour_partitions[FindMax(found, neighbour_weights)];
    }

    PrintDebugData("Constraint Partition", rConstraintPartition);
}

// ─── RedistributeHangingNodes ─────────────────────────────────────────────────

void KaHIPDivideHeterogeneousInputProcess::RedistributeHangingNodes(
    std::vector<idxtype>& rNodePartition,
    const std::vector<idxtype>& rGeometryPartition,
    const IO::ConnectivitiesContainerType& rGeometryConnectivities,
    const std::vector<idxtype>& rElementPartition,
    const IO::ConnectivitiesContainerType& rElementConnectivities,
    const std::vector<idxtype>& rConditionPartition,
    const IO::ConnectivitiesContainerType& rConditionConnectivities,
    const std::vector<idxtype>& rConstraintPartition,
    const IO::ConnectivitiesContainerType& rConstraintConnectivities)
{
    // Count how many times each node is used locally within its partition.
    std::vector<int> node_use_counts(rNodePartition.size(), 0);

    auto count_local_uses = [&](
        const IO::ConnectivitiesContainerType& connectivities,
        const std::vector<idxtype>& entity_partition)
    {
        for (std::size_t i = 0; i < entity_partition.size(); ++i) {
            for (std::size_t node_id : connectivities[i]) {
                if (rNodePartition[node_id - 1] == entity_partition[i]) {
                    node_use_counts[node_id - 1]++;
                }
            }
        }
    };

    count_local_uses(rGeometryConnectivities,  rGeometryPartition);
    count_local_uses(rElementConnectivities,   rElementPartition);
    count_local_uses(rConditionConnectivities, rConditionPartition);
    count_local_uses(rConstraintConnectivities, rConstraintPartition);

    // Identify hanging nodes (not locally used in their own partition)
    std::vector<std::size_t> hanging_nodes;
    for (std::size_t i = 0; i < node_use_counts.size(); ++i) {
        if (node_use_counts[i] == 0)
            hanging_nodes.push_back(i + 1); // 1-indexed node ID
    }

    KRATOS_INFO_IF("KaHIPDivideHeterogeneousInputProcess", mVerbosity > 0)
        << (hanging_nodes.empty()
            ? "No isolated (hanging) nodes found."
            : "Relocating " + std::to_string(hanging_nodes.size()) + " isolated nodes.")
        << std::endl;

    // For each hanging node, find the partition that references it most and move it there.
    for (std::size_t node_id : hanging_nodes) {
        std::vector<int> usage_by_partition(mNumberOfPartitions, 0);

        auto count_references = [&](
            const IO::ConnectivitiesContainerType& connectivities,
            const std::vector<idxtype>& entity_partition)
        {
            for (std::size_t i = 0; i < entity_partition.size(); ++i) {
                for (std::size_t nid : connectivities[i]) {
                    if (nid == node_id) usage_by_partition[entity_partition[i]]++;
                }
            }
        };

        count_references(rGeometryConnectivities,  rGeometryPartition);
        count_references(rElementConnectivities,   rElementPartition);
        count_references(rConditionConnectivities, rConditionPartition);
        count_references(rConstraintConnectivities, rConstraintPartition);

        const std::size_t destination = FindMax(mNumberOfPartitions, usage_by_partition);

        KRATOS_INFO_IF("KaHIPDivideHeterogeneousInputProcess", mVerbosity > 0)
            << "Sending node " << node_id << " to partition " << destination << std::endl;

        rNodePartition[node_id - 1] = static_cast<idxtype>(destination);
    }

    KRATOS_INFO_IF("KaHIPDivideHeterogeneousInputProcess", mVerbosity > 0 && !hanging_nodes.empty())
        << "Relocated " << hanging_nodes.size() << " isolated nodes." << std::endl;
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

KaHIPDivideHeterogeneousInputProcess::SizeType
KaHIPDivideHeterogeneousInputProcess::FindMax(SizeType NumTerms, const std::vector<int>& rVect)
{
    int max_val = rVect[0];
    SizeType imax = 0;
    for (SizeType i = 1; i < NumTerms; ++i) {
        if (rVect[i] > max_val) {
            max_val = rVect[i];
            imax = i;
        }
    }
    return imax;
}

void KaHIPDivideHeterogeneousInputProcess::PrintDebugData(
    const std::string& rLabel,
    const std::vector<idxtype>& rPartitionData)
{
    if (mVerbosity <= 1) return;

    KRATOS_INFO("KaHIPDivideHeterogeneousInputProcess") << rLabel << std::endl;
    for (SizeType p = 0; p < mNumberOfPartitions; ++p) {
        int count = 0;
        std::cout << "Partition " << p << ": ";
        for (SizeType i = 0; i < rPartitionData.size(); ++i) {
            if (rPartitionData[i] == static_cast<idxtype>(p)) {
                ++count;
                if (mVerbosity > 2) std::cout << (i + 1) << ",";
            }
        }
        std::cout << count << " objects." << std::endl;
    }
}

} // namespace Kratos
