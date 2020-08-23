//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors: Aditya Ghantasala
//

// System includes

/* Project includes */
#include "custom_utilities/gather_modelpart_on_all_ranks.h"
#ifdef KRATOS_USING_MPI
#include "mpi/includes/mpi_communicator.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#endif

namespace Kratos {

void GatherModelPartOnAllRanksUtility::GatherModelPartOnAllRanks(
    ModelPart &rModelPartToGather, ModelPart &rGatheredModelPart) {
  typedef ModelPart::NodesContainerType NodesContainerType;
  typedef ModelPart::ElementsContainerType ElementsContainerType;
  typedef ModelPart::ConditionsContainerType ConditionsContainerType;

  // GatherModelPartOnAllRanksUtility::TransferNodeGeometries(rModelPartToGather, rGatheredModelPart);

  const DataCommunicator &r_comm =
      rModelPartToGather.GetCommunicator().GetDataCommunicator();
  const int mpi_size = r_comm.Size();
  const int mpi_rank = r_comm.Rank();

  rGatheredModelPart.GetNodalSolutionStepVariablesList() =
      rModelPartToGather.GetNodalSolutionStepVariablesList();

  rGatheredModelPart.SetBufferSize(rModelPartToGather.GetBufferSize());

  // send everything to node with id "gather_rank"
  // transfer nodes
  std::vector<NodesContainerType> SendNodes(mpi_size);
  std::vector<NodesContainerType> RecvNodes(mpi_size);

  for (int dest_rank = 0; dest_rank < mpi_size; ++dest_rank) {
    SendNodes[dest_rank].reserve(rModelPartToGather.Nodes().size());
    for (NodesContainerType::iterator it = rModelPartToGather.NodesBegin();
         it != rModelPartToGather.NodesEnd(); ++it) {
      // only send the nodes owned by this partition
      if (it->FastGetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
        SendNodes[dest_rank].push_back(*it.base());
    }
  }
  rModelPartToGather.GetCommunicator().TransferObjects(SendNodes, RecvNodes);
  for (unsigned int i = 0; i < RecvNodes.size(); i++) {
    for (NodesContainerType::iterator it = RecvNodes[i].begin();
         it != RecvNodes[i].end(); ++it){
      // if (rGatheredModelPart.Nodes().find(it->Id()) ==
      //     rGatheredModelPart.Nodes().end())
      // KRATOS_INFO("Received node is : ")<<*it<<std::endl;
      rGatheredModelPart.Nodes().push_back(*it.base());
    }
  }
  int temp = rGatheredModelPart.Nodes().size();
  KRATOS_ERROR_IF(temp != int(rGatheredModelPart.Nodes().size()))
      << "the rGatheredModelPart has repeated nodes";
  SendNodes.clear();
  RecvNodes.clear();

  for (NodesContainerType::iterator it = rModelPartToGather.NodesBegin();
       it != rModelPartToGather.NodesEnd(); ++it) {
    rGatheredModelPart.Nodes().push_back(*it.base());
  }
  rGatheredModelPart.Nodes().Unique();

  // transfer elements
  std::vector<ElementsContainerType> SendElements(mpi_size);
  std::vector<ElementsContainerType> RecvElements(mpi_size);
  for (int dest_rank = 0; dest_rank < mpi_size; ++dest_rank) {
    SendElements[dest_rank].reserve(rModelPartToGather.Elements().size());
    for (ElementsContainerType::iterator it =
             rModelPartToGather.ElementsBegin();
         it != rModelPartToGather.ElementsEnd(); ++it) {
      SendElements[dest_rank].push_back(*it.base());
    }
  }
  rModelPartToGather.GetCommunicator().TransferObjects(SendElements,
                                                       RecvElements);

  for (unsigned int i = 0; i < RecvElements.size(); i++) {
    for (ElementsContainerType::iterator it = RecvElements[i].begin();
         it != RecvElements[i].end(); ++it) {
      // replace the nodes copied with the element by nodes
      // in the model part
      Element::GeometryType &rGeom = it->GetGeometry();
      unsigned int NumNodes = rGeom.PointsNumber();
      for (unsigned int iNode = 0; iNode < NumNodes; iNode++) {
        NodesContainerType::iterator itNode =
            rGatheredModelPart.Nodes().find(rGeom(iNode)->Id());
        if (itNode != rGatheredModelPart.Nodes().end())
          rGeom(iNode) = *itNode.base();
      }
      rGatheredModelPart.Elements().push_back(*it.base());
    }
  }
  SendElements.clear();
  RecvElements.clear();
  for (ElementsContainerType::iterator it = rModelPartToGather.ElementsBegin();
       it != rModelPartToGather.ElementsEnd(); ++it) {
    rGatheredModelPart.Elements().push_back(*it.base());
  }
  // rGatheredModelPart.Elements().Unique();

  // transfer conditions
  std::vector<ConditionsContainerType> SendConditions(mpi_size);
  std::vector<ConditionsContainerType> RecvConditions(mpi_size);
  for (int dest_rank = 0; dest_rank < mpi_size; ++dest_rank) {
    SendConditions[dest_rank].reserve(rModelPartToGather.Conditions().size());
    for (ConditionsContainerType::iterator it =
             rModelPartToGather.ConditionsBegin();
         it != rModelPartToGather.ConditionsEnd(); ++it) {
      SendConditions[dest_rank].push_back(*it.base());
    }
  }
  rModelPartToGather.GetCommunicator().TransferObjects(SendConditions,
                                                       RecvConditions);
  for (unsigned int i = 0; i < RecvConditions.size(); i++) {
    for (ConditionsContainerType::iterator it = RecvConditions[i].begin();
         it != RecvConditions[i].end(); ++it) {
      // replace the nodes copied with the condition by nodes
      // in the model part
      Condition::GeometryType &rGeom = it->GetGeometry();
      unsigned int NumNodes = rGeom.PointsNumber();
      for (unsigned int iNode = 0; iNode < NumNodes; iNode++) {
        NodesContainerType::iterator itNode =
            rGatheredModelPart.Nodes().find(rGeom(iNode)->Id());
        if (itNode != rGatheredModelPart.Nodes().end())
          rGeom(iNode) = *itNode.base();
      }
      rGatheredModelPart.Conditions().push_back(*it.base());
    }
  }
  SendConditions.clear();
  RecvConditions.clear();

  for (ConditionsContainerType::iterator it =
           rModelPartToGather.ConditionsBegin();
       it != rModelPartToGather.ConditionsEnd(); ++it) {
    rGatheredModelPart.Conditions().push_back(*it.base());
  }

  IndexType cond_id = 1;
  for (auto &condition : rGatheredModelPart.Conditions()) {
    condition.SetId(cond_id);
    cond_id++;
  }

  // rGatheredModelPart.Conditions().Unique();

  // Now check and delete all the double conditions.
  // The double conditions are basically the processor interfaces.
  hashmap n_faces_map;
  const int num_conditions =
      static_cast<int>(rGatheredModelPart.NumberOfConditions());
  const auto conditions_begin = rGatheredModelPart.ConditionsBegin();
// Fill map that counts number of faces for given set of nodes
#pragma omp parallel for
  for (int i_c = 0; i_c < num_conditions; ++i_c) {
    auto i_condition = conditions_begin + i_c;
    auto &cond_geo = i_condition->GetGeometry();

    vector<IndexType> ids(i_condition->GetGeometry().size());

    // Store node ids
    int i = 0;
    for (auto &node : cond_geo)
      ids[i++] = node.Id();

    //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
    std::sort(ids.begin(), ids.end());

// Fill the map
#pragma omp critical
    n_faces_map[ids] += 1;
  }

  std::vector<IndexType> condition_ids_to_remove;
  std::vector<IndexType> node_ids_to_remove;
  std::vector<IndexType> node_ids_to_keep;
  std::vector<IndexType> all_node_ids;
  all_node_ids.reserve(rGatheredModelPart.NumberOfNodes());
  node_ids_to_keep.reserve(rGatheredModelPart.NumberOfNodes());

  for (int i_c = 0; i_c < num_conditions; ++i_c) {
    auto i_condition = conditions_begin + i_c;
    auto &cond_geo = i_condition->GetGeometry();
    vector<IndexType> ids(i_condition->GetGeometry().size());

    // Store node ids
    int i = 0;
    for (auto &node : cond_geo)
      ids[i++] = node.Id();

    //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
    std::sort(ids.begin(), ids.end());

    if (n_faces_map[ids] > 1)
      condition_ids_to_remove.push_back(i_condition->Id());
  }

  // Now find out how to remove the nodes connected "only" to those conditions.
  for (const auto &cond : rGatheredModelPart.Conditions())
    for (const auto &node : cond.GetGeometry())
      all_node_ids.push_back(node.Id());
  std::set<IndexType> s(all_node_ids.begin(), all_node_ids.end());
  all_node_ids.assign(s.begin(), s.end());

  for (const int &cond_id : condition_ids_to_remove)
    rGatheredModelPart.RemoveCondition(cond_id);

  for (const auto &cond : rGatheredModelPart.Conditions())
    for (const auto &node : cond.GetGeometry())
      node_ids_to_keep.push_back(node.Id());

  std::set<IndexType> s1(node_ids_to_keep.begin(), node_ids_to_keep.end());
  node_ids_to_keep.assign(s1.begin(), s1.end());

  // now all_node_ids only has ids to remove
  // https://stackoverflow.com/questions/21195217/elegant-way-to-remove-all-elements-of-a-vector-that-are-contained-in-another-vec
  all_node_ids.erase(remove_if(begin(all_node_ids), end(all_node_ids),
                               [&](IndexType x) {
                                 return binary_search(begin(node_ids_to_keep),
                                                      end(node_ids_to_keep), x);
                               }),
                     end(all_node_ids));

  for (const int &node_id : all_node_ids)
    rGatheredModelPart.RemoveNode(node_id);

#ifdef KRATOS_USING_MPI
  if(r_comm.IsDistributed()){
    ModelPart &r_root_mp = rModelPartToGather.GetRootModelPart();
    Communicator::Pointer p_new_comm = Kratos::make_shared<MPICommunicator>(
        &rModelPartToGather.GetNodalSolutionStepVariablesList(),
        r_root_mp.GetCommunicator().GetDataCommunicator());
    rGatheredModelPart.SetCommunicator(p_new_comm);
    // ParallelFillCommunicator(rGatheredModelPart).Execute();
  }
#endif
}

void GatherModelPartOnAllRanksUtility::TransferNodeGeometries(
    ModelPart &rModelPartToGather, ModelPart &rGatheredModelPart)
    {
        KRATOS_INFO_ALL_RANKS("Start TransferNodeGeometries : ")<<std::endl;
        typedef ModelPart::NodesContainerType NodesContainerType;
        typedef std::vector<int> IdsVectorType;
        typedef std::vector<int> PidsVectorType;
        typedef std::vector<double> CoordsVectorType;

        const DataCommunicator &r_comm =
            rModelPartToGather.GetCommunicator().GetDataCommunicator();
        const int mpi_size = r_comm.Size();
        const int mpi_rank = r_comm.Rank();

        std::vector<IdsVectorType> NodesIds(mpi_size);
        std::vector<PidsVectorType> NodesPids(mpi_size);
        std::vector<CoordsVectorType> NodesCoords(mpi_size);


        NodesIds[mpi_rank].reserve(rModelPartToGather.NumberOfNodes());
        NodesPids[mpi_rank].reserve(rModelPartToGather.NumberOfNodes());
        NodesCoords[mpi_rank].reserve(3*rModelPartToGather.NumberOfNodes());
        for (NodesContainerType::iterator it = rModelPartToGather.NodesBegin();
            it != rModelPartToGather.NodesEnd(); ++it) {
                NodesIds[mpi_rank].push_back((*it).Id());
                NodesCoords[mpi_rank].push_back((*it).X());
                NodesCoords[mpi_rank].push_back((*it).Y());
                NodesCoords[mpi_rank].push_back((*it).Z());
                NodesPids[mpi_rank].push_back((*it).FastGetSolutionStepValue(PARTITION_INDEX));
        }

        for (int dest_rank = 0; dest_rank < mpi_size; ++dest_rank){
            r_comm.Barrier();
            r_comm.Broadcast(NodesIds[dest_rank], dest_rank);
        }
        r_comm.Barrier();
        KRATOS_INFO_ALL_RANKS("End TransferNodeIds : ")<<std::endl;

        for (int dest_rank = 0; dest_rank < mpi_size; ++dest_rank){
            r_comm.Barrier();
            r_comm.Broadcast(NodesPids[dest_rank], dest_rank);
        }
        r_comm.Barrier();
        KRATOS_INFO_ALL_RANKS("End TransferNodePids : ")<<std::endl;

        // for (int dest_rank = 0; dest_rank < mpi_size; ++dest_rank){
        //     r_comm.Barrier();
        //     r_comm.Broadcast(NodesCoords[dest_rank], dest_rank);
        // }
        // KRATOS_INFO_ALL_RANKS("End TransferNodeCoords : ")<<std::endl;
        // r_comm.Barrier();
    }
}