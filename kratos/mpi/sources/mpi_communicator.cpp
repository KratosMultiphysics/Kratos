//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//                   Jordi Cotela Dalmau
//

// System includes

// External includes

// Project includes
#include "mpi/includes/mpi_communicator.h"

namespace Kratos
{

// Helper function to retrieve dofs from nodes
Node::DofType::Pointer GetDof(
    VariableData::KeyType Key,
    Node& rNode)
{
    auto& r_node_dofs = rNode.GetDofs();
    for (auto i_dof = r_node_dofs.begin(); i_dof != r_node_dofs.end(); ++i_dof) {
        if ((*i_dof)->GetVariable().Key() == Key) {
            return i_dof->get();
        }
    }
    KRATOS_ERROR << "Dof for variable " << Key << " missing for node " << rNode.Id();
}

}

namespace Kratos {

template<>
void MPICommunicator::FillBuffer(
    std::vector<std::size_t>& rBuffer,
    MPICommunicator::MeshType& rSourceMesh,
    MPIInternals::DofSetAccess& rAccess)
{
    const auto& r_dof_data = rAccess.GetDofSet().GetContainer();
    auto& r_node_map = rAccess.GetNodeMap();
    const auto& r_source_nodes = rSourceMesh.Nodes();

    const std::size_t send_size_guess = (r_node_map.size() > 0) ? r_node_map[0].second + 1: 0;
    rBuffer.clear();
    rBuffer.reserve(r_source_nodes.size() * send_size_guess);

    for (auto iter = r_source_nodes.begin(); iter != r_source_nodes.end(); ++iter) {
        const auto& dof_data = r_node_map[iter->Id()];
        const std::size_t offset = dof_data.first;
        const std::size_t count = dof_data.second;

        // how many keys to read for this dof
        rBuffer.push_back(count);

        for (unsigned int i = 0; i < count; i++) {
            rBuffer.push_back(r_dof_data[offset + i]->GetVariable().Key());
        }
    }
}

template<>
void MPICommunicator::UpdateValues(
    const std::vector<size_t>& rBuffer,
    MPICommunicator::MeshType& rSourceMesh,
    MPIInternals::DofSetAccess& rAccess,
    MPICommunicator::Operation<OperationType::Combine>)
{
    auto& r_container = rSourceMesh.Nodes();
    auto& r_dof_set = rAccess.GetDofSet();
    std::size_t dof_count, key, position = 0;

    for (auto iter = r_container.begin(); iter != r_container.end(); ++iter) {
        dof_count = rBuffer[position++];
        for (std::size_t i = 0; i < dof_count; i++) {
            key = rBuffer[position++];
            r_dof_set.push_back(GetDof(key, *iter));
        }
    }

    KRATOS_WARNING_IF_ALL_RANKS("MPICommunicator", position > rBuffer.size())
    << GetDataCommunicator()
    << "Error in estimating receive buffer size." << std::endl;
}

bool MPICommunicator::SynchronizeDofSet(MPICommunicator::DofSetType& rDofSet)
{
    constexpr MeshAccess<DistributedType::Interface> interface_meshes;
    MPIInternals::DofSetAccess dof_set_access(rDofSet);
    MPICommunicator::Operation<OperationType::Combine> combine_dofs;

    TransferDistributedValuesUnknownSize(interface_meshes, interface_meshes, dof_set_access, combine_dofs);

    rDofSet.Unique();
    return true;
}

}