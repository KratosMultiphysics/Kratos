#include "mpi/utilities/gather_modelpart_utility.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "includes/data_communicator.h"
#include "mpi/includes/mpi_communicator.h"

namespace Kratos
{

GatherModelPartUtility::GatherModelPartUtility(int gather_rank,
                                               ModelPart& origin_model_part,
                                               int mesh_id,
                                               ModelPart& destination_model_part)
    : mr_model_part(destination_model_part), mgather_rank(gather_rank)
{
    KRATOS_TRY;

    const DataCommunicator& r_comm =
        origin_model_part.GetCommunicator().GetDataCommunicator();
    const int mpi_rank = r_comm.Rank();
    const int mpi_size = r_comm.Size();

    destination_model_part.GetNodalSolutionStepVariablesList() =
        origin_model_part.GetNodalSolutionStepVariablesList();
    if (r_comm.IsDistributed())
    {
      VariablesList* pVariablesList =
          &destination_model_part.GetNodalSolutionStepVariablesList();
      destination_model_part.SetCommunicator(
          Communicator::Pointer(new MPICommunicator(pVariablesList, r_comm)));
    }
    destination_model_part.SetBufferSize(origin_model_part.GetBufferSize());

    // copy the mesh of interest to destination_model_part
    // be careful to push back the pointer and not copy
    // construct the object
    for (NodesContainerType::iterator it = origin_model_part.GetMesh(mesh_id).NodesBegin();
         it != origin_model_part.GetMesh(mesh_id).NodesEnd(); ++it)
    {
        destination_model_part.Nodes().push_back(*it.base());
    }

    for (ElementsContainerType::iterator it =
             origin_model_part.GetMesh(mesh_id).ElementsBegin();
         it != origin_model_part.GetMesh(mesh_id).ElementsEnd(); ++it)
    {
        destination_model_part.Elements().push_back(*it.base());
    }

    for (ConditionsContainerType::iterator it =
             origin_model_part.GetMesh(mesh_id).ConditionsBegin();
         it != origin_model_part.GetMesh(mesh_id).ConditionsEnd(); ++it)
    {
        destination_model_part.Conditions().push_back(*it.base());
    }

    // send everything to node with id "gather_rank"
    // transfer nodes
    std::vector<NodesContainerType> SendNodes(mpi_size);
    std::vector<NodesContainerType> RecvNodes(mpi_size);
    SendNodes[gather_rank].reserve(destination_model_part.Nodes().size());
    if (r_comm.IsDistributed())
    {
      for (NodesContainerType::iterator it = destination_model_part.NodesBegin();
           it != destination_model_part.NodesEnd(); ++it)
      {
          // only send the nodes owned by this partition
          if (it->FastGetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
              SendNodes[gather_rank].push_back(*it.base());
      }
    }
    else
    {
      for (NodesContainerType::iterator it = destination_model_part.NodesBegin();
           it != destination_model_part.NodesEnd(); ++it)
      {
          SendNodes[gather_rank].push_back(*it.base());
      }
    }

    destination_model_part.GetCommunicator().TransferObjects(SendNodes, RecvNodes);
    for (unsigned int i = 0; i < RecvNodes.size(); i++)
    {
        for (NodesContainerType::iterator it = RecvNodes[i].begin();
             it != RecvNodes[i].end(); ++it)
            if (destination_model_part.Nodes().find(it->Id()) ==
                destination_model_part.Nodes().end())
                destination_model_part.Nodes().push_back(*it.base());
    }
    int temp = destination_model_part.Nodes().size();
    destination_model_part.Nodes().Unique();
    KRATOS_ERROR_IF(temp != int(destination_model_part.Nodes().size()))
        << "the destination_model_part has repeated nodes";
    SendNodes.clear();
    RecvNodes.clear();

    // transfer elements
    std::vector<ElementsContainerType> SendElements(mpi_size);
    std::vector<ElementsContainerType> RecvElements(mpi_size);
    SendElements[gather_rank].reserve(destination_model_part.Elements().size());
    for (ElementsContainerType::iterator it = destination_model_part.ElementsBegin();
         it != destination_model_part.ElementsEnd(); ++it)
    {
        SendElements[gather_rank].push_back(*it.base());
    }
    destination_model_part.GetCommunicator().TransferObjects(SendElements, RecvElements);
    for (unsigned int i = 0; i < RecvElements.size(); i++)
    {
        for (ElementsContainerType::iterator it = RecvElements[i].begin();
             it != RecvElements[i].end(); ++it)
        {
            // replace the nodes copied with the element by nodes
            // in the model part
            Element::GeometryType& rGeom = it->GetGeometry();
            unsigned int NumNodes = rGeom.PointsNumber();
            for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
            {
                NodesContainerType::iterator itNode =
                    destination_model_part.Nodes().find(rGeom(iNode)->Id());
                if (itNode != destination_model_part.Nodes().end())
                    rGeom(iNode) = *itNode.base();
            }
            destination_model_part.Elements().push_back(*it.base());
        }
    }
    SendElements.clear();
    RecvElements.clear();

    // transfer conditions
    std::vector<ConditionsContainerType> SendConditions(mpi_size);
    std::vector<ConditionsContainerType> RecvConditions(mpi_size);
    SendConditions[gather_rank].reserve(destination_model_part.Conditions().size());
    for (ConditionsContainerType::iterator it = destination_model_part.ConditionsBegin();
         it != destination_model_part.ConditionsEnd(); ++it)
    {
        SendConditions[gather_rank].push_back(*it.base());
    }
    destination_model_part.GetCommunicator().TransferObjects(SendConditions, RecvConditions);
    for (unsigned int i = 0; i < RecvConditions.size(); i++)
    {
        for (ConditionsContainerType::iterator it = RecvConditions[i].begin();
             it != RecvConditions[i].end(); ++it)
        {
            // replace the nodes copied with the condition by nodes
            // in the model part
            Condition::GeometryType& rGeom = it->GetGeometry();
            unsigned int NumNodes = rGeom.PointsNumber();
            for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
            {
                NodesContainerType::iterator itNode =
                    destination_model_part.Nodes().find(rGeom(iNode)->Id());
                if (itNode != destination_model_part.Nodes().end())
                    rGeom(iNode) = *itNode.base();
            }
            destination_model_part.Conditions().push_back(*it.base());
        }
    }
    SendConditions.clear();
    RecvConditions.clear();

    if (r_comm.IsDistributed())
    {
      ParallelFillCommunicator(destination_model_part, r_comm).Execute();
    }

    KRATOS_CATCH("");
}

void GatherModelPartUtility::GatherOnMaster()
{
    KRATOS_TRY;
    mr_model_part.GetCommunicator().SynchronizeNodalSolutionStepsData();
    KRATOS_CATCH("");
}

template <class TDataType>
void GatherModelPartUtility::GatherOnMaster(const Variable<TDataType>& ThisVariable)
{
    KRATOS_TRY;
    mr_model_part.GetCommunicator().SynchronizeVariable(ThisVariable);
    KRATOS_CATCH("");
}

template <class TDataType>
void GatherModelPartUtility::ScatterFromMaster(const Variable<TDataType>& ThisVariable)
{
    KRATOS_TRY;

    Communicator& r_comm = mr_model_part.GetCommunicator();

    if (mgather_rank != r_comm.GetDataCommunicator().Rank())
    {
        for (NodesContainerType::iterator it = mr_model_part.NodesBegin();
             it != mr_model_part.NodesEnd(); ++it)
            it->FastGetSolutionStepValue(ThisVariable) = ThisVariable.Zero();
    }
    r_comm.AssembleCurrentData(ThisVariable);

    KRATOS_CATCH("");
}

template void GatherModelPartUtility::GatherOnMaster(const Variable<double>&);
template void GatherModelPartUtility::GatherOnMaster(const Variable<array_1d<double, 3>>&);
template void GatherModelPartUtility::ScatterFromMaster(const Variable<double>&);
template void GatherModelPartUtility::ScatterFromMaster(const Variable<array_1d<double, 3>>&);

namespace Internals
{

}

} // namespace Kratos.
