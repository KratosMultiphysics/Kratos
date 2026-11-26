//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "mpi/includes/mpi_communicator.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_io/hdf5_connectivities_data.h"
#include "custom_io/hdf5_points_data.h"
#include "custom_utilities/container_io_utils.h"
#include "custom_utilities/factor_elements_and_conditions_utility.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"

// Include base h
#include "hdf5_partitioned_model_part_io.h"

namespace Kratos
{
namespace HDF5
{
PartitionedModelPartIO::PartitionedModelPartIO(
    Parameters Settings,
    File::Pointer pFile)
    : ModelPartIO(Settings, pFile)
{
    KRATOS_TRY;

    Check();

    KRATOS_CATCH("");
}

bool PartitionedModelPartIO::ReadNodes(NodesContainerType& rNodes)
{
    KRATOS_TRY;

    rNodes.clear();

    // read local nodes
    BaseType::ReadNodes(rNodes);

    // Read ghost nodes.
    Internals::PointsData<Internals::NodesIO> ghost_points(mPrefix + "/Nodes/Ghost", mpFile);
    ghost_points.Read(rNodes, Internals::NodesIO{});

    return true;

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::WriteNodes(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    NodesContainerType local_nodes;
    NodesContainerType ghost_nodes;

    const unsigned num_nodes = rNodes.size();
    local_nodes.reserve(num_nodes);
    ghost_nodes.reserve(0.1 * num_nodes);

    // Divide nodes into local and global containers.
    int my_pid = mpFile->GetPID();
    for (auto it_node = rNodes.begin(); it_node != rNodes.end(); ++it_node) {
        if (it_node->FastGetSolutionStepValue(PARTITION_INDEX) == my_pid) {
            local_nodes.push_back(*it_node.base());
        } else {
            ghost_nodes.push_back(*it_node.base());
        }
    }

    // Write local nodes.
    BaseType::WriteNodes(local_nodes);

    // Write ghost nodes.
    Internals::PointsData<Internals::NodesIO> ghost_points(mPrefix + "/Nodes/Ghost", mpFile);
    ghost_points.Write(ghost_nodes, Internals::NodesIO{}, Parameters(R"({})"));
    WritePartitionIndex(mPrefix + "/Nodes/Ghost", ghost_nodes);

    // Write partition index
    WritePartitionIndex(mPrefix + "/Nodes/Ghost", ghost_nodes);

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::Check()
{
    if (mpFile->GetTotalProcesses() == 1) {
        KRATOS_ERROR << "Using PartitionedModelPartIO with single process file access." << std::endl;
    }
}

void PartitionedModelPartIO::ReadParitionIndices(ModelPart& rModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(PARTITION_INDEX))
        << "PARTITION_INDEX variable not found in nodal solution step variables list of "
        << rModelPart.FullName() << ".";

    // get the communicator
    auto& r_communicator = rModelPart.GetCommunicator();
    const auto my_pid = r_communicator.MyPID();

    auto& r_nodes = rModelPart.Nodes();

    // set the partition index of all nodes first to the current rank
    // since we don't know which nodes belongs to ghost or local yet.
    block_for_each(r_nodes, [my_pid](auto& rNode) {
        rNode.FastGetSolutionStepValue(PARTITION_INDEX) = my_pid;
    });

    // now we read the ghost node partition index from the H5 file.
    const auto& r_ghost_node_path = mPrefix + "/Nodes/Ghost";
    const auto [start_index, block_size] = StartIndexAndBlockSize(*mpFile, r_ghost_node_path);
    Vector<int> partition_ids, node_ids;
    mpFile->ReadDataSet(r_ghost_node_path + "/PARTITION_INDEX", partition_ids, start_index, block_size);
    mpFile->ReadDataSet(r_ghost_node_path + "/Ids", node_ids, start_index, block_size);

    // now we correct PARTITION_INDEX of ghost nodes.
    IndexPartition<IndexType>(node_ids.size()).for_each([&r_nodes, &node_ids, &partition_ids](const auto Index) {
        r_nodes[node_ids[Index]].FastGetSolutionStepValue(PARTITION_INDEX) = partition_ids[Index];
    });

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::SetCommunicator(ModelPart& rModelPart) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mpFile->GetDataCommunicator().IsDistributed())
        << "Trying to use PartitionedModelPartIO with an HDF5File which is based on serial communicator.";

    // set the mpi communicator
    if (!rModelPart.GetCommunicator().IsDistributed()) {
        rModelPart.SetCommunicator(Kratos::make_shared<MPICommunicator>(rModelPart.pGetNodalSolutionStepVariablesList().get(), mpFile->GetDataCommunicator()));
    }

    // get the communicator
    auto& r_communicator = rModelPart.GetCommunicator();
    const auto my_pid = r_communicator.MyPID();

    // now we have set all the PARTITION_INDEX values correctly for all nodes,
    // now we set the local mesh and ghost mesh nodes correctly.
    auto& r_local_nodes = r_communicator.LocalMesh().Nodes();
    auto& r_ghost_nodes = r_communicator.GhostMesh().Nodes();

    for (auto p_itr = rModelPart.Nodes().ptr_begin(); p_itr != rModelPart.Nodes().ptr_end(); ++p_itr) {
        if ((*p_itr)->FastGetSolutionStepValue(PARTITION_INDEX) == my_pid) {
            r_local_nodes.insert(r_local_nodes.end(), *p_itr);
        } else {
            r_ghost_nodes.insert(r_ghost_nodes.end(), *p_itr);
        }
    }

    // now populate sub model parts
    for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
        this->SetCommunicator(r_sub_model_part);
    }

    KRATOS_CATCH("");
}

void PartitionedModelPartIO::WritePartitionIndex(
    const std::string& rPath,
    const NodesContainerType& rGhostNodes)
{
    KRATOS_TRY;

    Vector<int> partition_ids(rGhostNodes.size());

    IndexPartition<IndexType>(rGhostNodes.size()).for_each([&rGhostNodes, &partition_ids](const auto Index) {
        partition_ids[Index] = (rGhostNodes.begin() + Index)->FastGetSolutionStepValue(PARTITION_INDEX);
    });

    WriteInfo info;
    mpFile->WriteDataSet(rPath + "/PARTITION_INDEX", partition_ids, info);

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.