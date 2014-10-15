//   Project Name:        Kratos
//   Last Modified by:    $Author: Pablo $
//   Date:                $Date: 2010-05-12 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_GATHER_MODELPART_UTILITY)
#define  KRATOS_GATHER_MODELPART_UTILITY


#ifdef _OPENMP
#include <omp.h>
#endif




// System includes

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "parallel_fill_communicator.h"



namespace Kratos
{
class GatherModelPartUtility
{
public:

    ///This function is designed to obtain data from "origin_model_part.GetMesh(mesh_id)", copy it to a new model part
    ///and make rank "gather_rank" to have a copy of it. Transferred nodes will be treated as ghost on the gather_rank
    ///@param gather_rank --> mpi rank to which the model part is gathered
    ///@param origin_model_part --> model part on which the origin mesh is contained
    ///@param mesh_id --> id of the mesh which contains the data
    ///@param destination_model_part --> model part to which we gather the data
    GatherModelPartUtility(int gather_rank, ModelPart& origin_model_part, int mesh_id, ModelPart& destination_model_part)
        :mr_model_part(destination_model_part), mgather_rank(gather_rank)
    {
        KRATOS_TRY

        int mpi_rank;
        int mpi_size;

        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        //fill the communicator
        VariablesList * mVariables_List = &origin_model_part.GetNodalSolutionStepVariablesList();
        destination_model_part.SetCommunicator(Communicator::Pointer(new MPICommunicator(mVariables_List)));

	destination_model_part.GetNodalSolutionStepVariablesList() = *mVariables_List;
	destination_model_part.SetBufferSize(origin_model_part.GetBufferSize());

        //*********************************************************************************
        //copy the mesh of interest to destination_model_part (on each mpi processor)

        for(ModelPart::NodesContainerType::iterator it = origin_model_part.GetMesh(mesh_id).NodesBegin();
                it != origin_model_part.GetMesh(mesh_id).NodesEnd(); it++)
        {
            destination_model_part.Nodes().push_back(*it.base());
        }

        for(ModelPart::ElementsContainerType::iterator it = origin_model_part.GetMesh(mesh_id).ElementsBegin();
                it != origin_model_part.GetMesh(mesh_id).ElementsEnd(); it++)
        {
            destination_model_part.Elements().push_back(*it.base());
        }

        for(ModelPart::ConditionsContainerType::iterator it = origin_model_part.GetMesh(mesh_id).ConditionsBegin();
                it != origin_model_part.GetMesh(mesh_id).ConditionsEnd(); it++)
        {
            destination_model_part.Conditions().push_back(*it.base());
        }


        //*********************************************************************************
        //send everything to node with id "gather_rank"

        //transfer nodes --> NOTE THAT WE ALSO TRANSFER RANK 0 to RANK 0!!
        std::vector< ModelPart::NodesContainerType > SendNodes(mpi_size);
        std::vector< ModelPart::NodesContainerType > RecvNodes(mpi_size);
        SendNodes[gather_rank].reserve(destination_model_part.Nodes().size());
        for(ModelPart::NodesContainerType::iterator it = destination_model_part.NodesBegin();
                it != destination_model_part.NodesEnd(); it++)
        {
            if(it->FastGetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
                SendNodes[gather_rank].push_back(*it.base());
        }
        destination_model_part.GetCommunicator().TransferObjects(SendNodes,RecvNodes);
        for(unsigned int i=0; i<RecvNodes.size(); i++)
        {
            for(ModelPart::NodesContainerType::iterator it = RecvNodes[i].begin(); it != RecvNodes[i].end(); it++)
	      if (destination_model_part.Nodes().find((*it).Id()) == destination_model_part.Nodes().end())
                destination_model_part.Nodes().push_back(*it);
        }
	int temp = destination_model_part.Nodes().size();
	destination_model_part.Nodes().Unique();
	if(temp != static_cast<int>(destination_model_part.Nodes().size()))
	  KRATOS_ERROR(std::logic_error,"the destination_model_part has repeated nodes","");
        SendNodes.clear();
        RecvNodes.clear();

        //transfer elements
        std::vector<ModelPart::ElementsContainerType > SendElements(mpi_size);
        std::vector<ModelPart::ElementsContainerType > RecvElements(mpi_size);
        SendElements[gather_rank].reserve(destination_model_part.Elements().size());
        for(ModelPart::ElementsContainerType::iterator it = destination_model_part.ElementsBegin();
                it != destination_model_part.ElementsEnd(); it++)
        {
            SendElements[gather_rank].push_back(*it.base());
        }
        destination_model_part.GetCommunicator().TransferObjects(SendElements,RecvElements);
        for(unsigned int i=0; i<RecvElements.size(); i++)
        {
            for(ModelPart::ElementsContainerType::iterator it = RecvElements[i].begin(); it != RecvElements[i].end(); it++)
                destination_model_part.Elements().push_back(*it);
        }
        SendElements.clear();
        RecvElements.clear();

        //transfer conditions
        std::vector< ModelPart::ConditionsContainerType > SendConditions(mpi_size);
        std::vector< ModelPart::ConditionsContainerType > RecvConditions(mpi_size);
        SendConditions[gather_rank].reserve(destination_model_part.Conditions().size());
        for(ModelPart::ConditionsContainerType::iterator it = destination_model_part.ConditionsBegin();
                it != destination_model_part.ConditionsEnd(); it++)
        {
            SendConditions[gather_rank].push_back(*it.base());
        }
        destination_model_part.GetCommunicator().TransferObjects(SendConditions,RecvConditions);
        for(unsigned int i=0; i<RecvConditions.size(); i++)
        {
            for(ModelPart::ConditionsContainerType::iterator it = RecvConditions[i].begin(); it != RecvConditions[i].end(); it++)
                destination_model_part.Conditions().push_back(*it);
        }
        SendConditions.clear();
        RecvConditions.clear();


        ParallelFillCommunicator(destination_model_part).Execute();


        KRATOS_CATCH("")
    }

    void GatherOnMaster()
    {
        KRATOS_TRY
        mr_model_part.GetCommunicator().SynchronizeNodalSolutionStepsData();
        KRATOS_CATCH("")
    }

    template< class TDataType >
    void GatherOnMaster(Variable<TDataType>& ThisVariable)
    {
      KRATOS_TRY
	mr_model_part.GetCommunicator().SynchronizeVariable(ThisVariable);
      KRATOS_CATCH("")
    }

    template< class TDataType >
    void ScatterFromMaster(Variable<TDataType>& ThisVariable)
    {
        KRATOS_TRY

        if(mgather_rank != mr_model_part.GetCommunicator().MyPID() )
        {
            for(ModelPart::NodesContainerType::iterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); it++)
                it->FastGetSolutionStepValue(ThisVariable) = ThisVariable.Zero();
        }
        mr_model_part.GetCommunicator().AssembleCurrentData(ThisVariable);

        KRATOS_CATCH("")
    }


    virtual ~GatherModelPartUtility()
    {
    }

private:
    ModelPart& mr_model_part;
    int mgather_rank;

};

} // namespace Kratos.

#endif // KRATOS_GATHER_MODELPART_UTILITY  defined 


