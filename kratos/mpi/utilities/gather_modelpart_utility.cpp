//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Michael Andre, https://github.com/msandre
//                   Jordi Cotela Dalmau
//

// System includes

// External includes

// Project includes
#include "mpi/utilities/gather_modelpart_utility.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "includes/data_communicator.h"
#include "mpi/includes/mpi_communicator.h"

namespace Kratos
{

GatherModelPartUtility::GatherModelPartUtility(
    const int GatherRank,
    ModelPart& rOriginModelPart,
    const int MeshId,
    ModelPart& rDestinationModelPart
    ) : mrModelPart(rDestinationModelPart), mGatherRank(GatherRank)
{
    KRATOS_TRY;

    const DataCommunicator& r_comm = rOriginModelPart.GetCommunicator().GetDataCommunicator();
    const int mpi_rank = r_comm.Rank();
    const int mpi_size = r_comm.Size();

    rDestinationModelPart.GetNodalSolutionStepVariablesList() = rOriginModelPart.GetNodalSolutionStepVariablesList();
    if (r_comm.IsDistributed()) {
      VariablesList* pVariablesList = &rDestinationModelPart.GetNodalSolutionStepVariablesList();
      rDestinationModelPart.SetCommunicator(Communicator::Pointer(new MPICommunicator(pVariablesList, r_comm)));
    }
    rDestinationModelPart.SetBufferSize(rOriginModelPart.GetBufferSize());

    // Copy the mesh of interest to rDestinationModelPart
    // be careful to push back the pointer and not copy
    // construct the object
    auto& r_mesh = rOriginModelPart.GetMesh(MeshId);
    for (auto it = r_mesh.NodesBegin(); it != r_mesh.NodesEnd(); ++it) {
        rDestinationModelPart.Nodes().push_back(*it.base());
    }

    for (auto it = r_mesh.ElementsBegin();  it != r_mesh.ElementsEnd(); ++it) {
        rDestinationModelPart.Elements().push_back(*it.base());
    }

    for (auto it = r_mesh.ConditionsBegin(); it != r_mesh.ConditionsEnd(); ++it) {
        rDestinationModelPart.Conditions().push_back(*it.base());
    }

    // send everything to node with id "GatherRank"
    // transfer nodes
    std::vector<NodesContainerType> SendNodes(mpi_size);
    std::vector<NodesContainerType> RecvNodes(mpi_size);
    SendNodes[GatherRank].reserve(rDestinationModelPart.Nodes().size());
    if (r_comm.IsDistributed()) {
      for (auto it = rDestinationModelPart.NodesBegin(); it != rDestinationModelPart.NodesEnd(); ++it) {
          // only send the nodes owned by this partition
          if (it->FastGetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
              SendNodes[GatherRank].push_back(*it.base());
      }
    } else {
      for (auto it = rDestinationModelPart.NodesBegin(); it != rDestinationModelPart.NodesEnd(); ++it) {
          SendNodes[GatherRank].push_back(*it.base());
      }
    }

    rDestinationModelPart.GetCommunicator().TransferObjects(SendNodes, RecvNodes);
    for (unsigned int i = 0; i < RecvNodes.size(); i++) {
        for (auto it = RecvNodes[i].begin(); it != RecvNodes[i].end(); ++it) {
            if (rDestinationModelPart.Nodes().find(it->Id()) == rDestinationModelPart.Nodes().end())
                rDestinationModelPart.Nodes().push_back(*it.base());
        }
    }
    int temp = rDestinationModelPart.Nodes().size();
    rDestinationModelPart.Nodes().Unique();
    KRATOS_ERROR_IF(temp != int(rDestinationModelPart.Nodes().size())) << "The rDestinationModelPart has repeated nodes" << std::endl;
    SendNodes.clear();
    RecvNodes.clear();

    // Transfer elements
    std::vector<ElementsContainerType> SendElements(mpi_size);
    std::vector<ElementsContainerType> RecvElements(mpi_size);
    SendElements[GatherRank].reserve(rDestinationModelPart.Elements().size());
    for (auto it = rDestinationModelPart.ElementsBegin(); it != rDestinationModelPart.ElementsEnd(); ++it) {
        SendElements[GatherRank].push_back(*it.base());
    }
    rDestinationModelPart.GetCommunicator().TransferObjects(SendElements, RecvElements);
    for (unsigned int i = 0; i < RecvElements.size(); i++) {
        for (auto it = RecvElements[i].begin(); it != RecvElements[i].end(); ++it) {
            // Replace the nodes copied with the element by nodes in the model part
            Element::GeometryType& rGeom = it->GetGeometry();
            unsigned int NumNodes = rGeom.PointsNumber();
            for (unsigned int iNode = 0; iNode < NumNodes; iNode++) {
                auto itNode = rDestinationModelPart.Nodes().find(rGeom(iNode)->Id());
                if (itNode != rDestinationModelPart.Nodes().end())
                    rGeom(iNode) = *itNode.base();
            }
            rDestinationModelPart.Elements().push_back(*it.base());
        }
    }
    SendElements.clear();
    RecvElements.clear();

    // Transfer conditions
    std::vector<ConditionsContainerType> SendConditions(mpi_size);
    std::vector<ConditionsContainerType> RecvConditions(mpi_size);
    SendConditions[GatherRank].reserve(rDestinationModelPart.Conditions().size());
    for (auto it = rDestinationModelPart.ConditionsBegin(); it != rDestinationModelPart.ConditionsEnd(); ++it) {
        SendConditions[GatherRank].push_back(*it.base());
    }
    rDestinationModelPart.GetCommunicator().TransferObjects(SendConditions, RecvConditions);
    for (unsigned int i = 0; i < RecvConditions.size(); i++) {
        for (auto it = RecvConditions[i].begin(); it != RecvConditions[i].end(); ++it) {
            // Replace the nodes copied with the condition by nodes in the model part
            Condition::GeometryType& rGeom = it->GetGeometry();
            unsigned int NumNodes = rGeom.PointsNumber();
            for (unsigned int iNode = 0; iNode < NumNodes; iNode++) {
                auto itNode = rDestinationModelPart.Nodes().find(rGeom(iNode)->Id());
                if (itNode != rDestinationModelPart.Nodes().end())
                    rGeom(iNode) = *itNode.base();
            }
            rDestinationModelPart.Conditions().push_back(*it.base());
        }
    }
    SendConditions.clear();
    RecvConditions.clear();

    if (r_comm.IsDistributed()) {
        ParallelFillCommunicator(rDestinationModelPart, r_comm).Execute();
    }

    KRATOS_CATCH("");
}

void GatherModelPartUtility::GatherOnMaster()
{
    KRATOS_TRY;
    mrModelPart.GetCommunicator().SynchronizeNodalSolutionStepsData();
    KRATOS_CATCH("");
}

template <class TDataType>
void GatherModelPartUtility::GatherOnMaster(const Variable<TDataType>& ThisVariable)
{
    KRATOS_TRY;
    mrModelPart.GetCommunicator().SynchronizeVariable(ThisVariable);
    KRATOS_CATCH("");
}

template <class TDataType>
void GatherModelPartUtility::ScatterFromMaster(const Variable<TDataType>& ThisVariable)
{
    KRATOS_TRY;

    Communicator& r_comm = mrModelPart.GetCommunicator();

    if (mGatherRank != r_comm.GetDataCommunicator().Rank()) {
        for (auto it = mrModelPart.NodesBegin(); it != mrModelPart.NodesEnd(); ++it)
            it->FastGetSolutionStepValue(ThisVariable) = ThisVariable.Zero();
    }
    r_comm.AssembleCurrentData(ThisVariable);

    KRATOS_CATCH("");
}

template void GatherModelPartUtility::GatherOnMaster(const Variable<double>&);
template void GatherModelPartUtility::GatherOnMaster(const Variable<array_1d<double, 3>>&);
template void GatherModelPartUtility::ScatterFromMaster(const Variable<double>&);
template void GatherModelPartUtility::ScatterFromMaster(const Variable<array_1d<double, 3>>&);

std::string GatherModelPartUtility::Info() const
{
    std::stringstream buffer;
    buffer << "GatherModelPartUtility";
    return buffer.str();
}

void GatherModelPartUtility::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "GatherModelPartUtility" << std::endl;
}

void GatherModelPartUtility::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
