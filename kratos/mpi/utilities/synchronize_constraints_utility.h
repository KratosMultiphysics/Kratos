//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala

#if !defined(MPI_CONSTRAINTS_UTILITY)
#define MPI_CONSTRAINTS_UTILITY

// System includes
#include <vector>


// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/global_pointer_utilities.h"
#include "includes/linear_master_slave_constraint.h"
#include "utilities/pointer_communicator.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Utility for calculating the Distance on a given modelpart
class MpiConstraintsUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MpiConstraintsUtility
    KRATOS_CLASS_POINTER_DEFINITION(MpiConstraintsUtility);
    typedef typename ModelPart::DofType DofType;
    typedef typename MasterSlaveConstraint::DofPointerVectorType DofPointerVectorType;
    typedef typename MasterSlaveConstraint::DofGlobalPointerType DofGlobalPointerType;
    typedef typename MasterSlaveConstraint::DofGlobalPointerVectorType DofGlobalPointerVectorType;
    typedef typename ModelPart::NodeType NodeType;
    typedef typename ModelPart::NodesContainerType NodesContainerType;
    typedef GlobalPointer<NodeType> GlobalPointerNodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MpiConstraintsUtility(ModelPart& rModelPart):mrModelPart(rModelPart)
    {

    }

    /// Destructor.
    ~MpiConstraintsUtility() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    template <class TVariableType>
    void AddConstraint(TVariableType &rVariable,
                        const IndexType ConstraintId,
                        const IndexType SlaveNodeId,
                        const IndexType MasterNodeId,
                        const double Weight,
                        const double Constant)
    {
        if(mrModelPart.GetCommunicator().LocalMesh().HasNode(SlaveNodeId))
            mMasterSlaveDetailsVector.push_back(
                Kratos::make_shared<MasterSlaveDetails<TVariableType>>(rVariable, ConstraintId, SlaveNodeId, MasterNodeId, Weight, Constant)
                );
    }


    // This is where all the nodes are transferred and actual constraints are created locally where slave exists.
    void SynchronizeAndCreateConstraints()
    {
        GetNodeIdsToSynchronize();
        CreateConstraints();
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name private Member Variables
    ///@{

    class BaseMasterSlaveDetails
    {
    public:
        BaseMasterSlaveDetails(IndexType ConstraintId, IndexType SlaveNodeId, IndexType  MasterNodeId, IndexType Weight, IndexType Constant)
        {
            mSlaveNodeId = SlaveNodeId;
            mMasterNodeId = MasterNodeId;
            mWeight = Weight;
            mConstant = Constant;
            mId = ConstraintId;
        }
        virtual IndexType GetID(){return mId;}
        virtual IndexType SlaveNodeId(){return mSlaveNodeId;}
        virtual IndexType MasterNodeId(){return mMasterNodeId;}
        virtual double Constant(){return mConstant;}
        virtual double Weight(){return mWeight;}
        virtual IndexType GetVariableKey()=0;
    private:
        IndexType mId;
        IndexType mSlaveNodeId;
        IndexType mMasterNodeId;
        double mWeight;
        double mConstant;
    };


    template <class TVariableType>
    class MasterSlaveDetails : public BaseMasterSlaveDetails
    {
    public:
        typedef TVariableType VaribleType;
        MasterSlaveDetails(TVariableType& rVariable, IndexType ConstraintId, IndexType SlaveNodeId, IndexType  MasterNodeId, IndexType Weight, IndexType Constant) 
        :BaseMasterSlaveDetails(ConstraintId,SlaveNodeId,MasterNodeId,Weight,Constant), mrVariable(rVariable)
        {
        }

        IndexType GetVariableKey() override {return mrVariable.Key();}

        TVariableType& mrVariable;
    };

    std::vector<Kratos::shared_ptr<BaseMasterSlaveDetails>> mMasterSlaveDetailsVector;
    std::vector<int> mNodeIdsToSync;
    ModelPart& mrModelPart;
    DofPointerVectorType mDofPointersVector;

    ///@}
    ///@name private operations
    ///@{

    void GetNodeIdsToSynchronize()
    {
        for(const auto& constraint_info : mMasterSlaveDetailsVector)
        {
            const auto& r_slave_node_id = constraint_info->SlaveNodeId();
            const auto& r_master_node_id = constraint_info->MasterNodeId();
            mNodeIdsToSync.push_back(r_slave_node_id);
            mNodeIdsToSync.push_back(r_master_node_id);
        }
    }



    void CreateConstraints()
    {
        typedef std::pair<IndexType, IndexType> DofIdVarKeyPairType;
        typedef GlobalPointer<DofIdVarKeyPairType> GlobalPointerDofIdVarKeyPairType;
        typedef GlobalPointersVector<DofIdVarKeyPairType> GlobalPointerVectorDofIdVarKeyPairType;
        GlobalPointerVectorDofIdVarKeyPairType global_ptrs_vectors_dof_id_key;
        auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();
        IndexType current_rank = r_data_communicator.Rank();
        IndexType comm_size = r_data_communicator.Size();
        auto remote_nodes_gps_map = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(mrModelPart.Nodes(), mNodeIdsToSync, r_data_communicator);
        LinearMasterSlaveConstraint::DofPointerVectorType slave_dof_vector;
        LinearMasterSlaveConstraint::DofPointerVectorType master_dof_vector;
        LinearMasterSlaveConstraint::MatrixType relation_matrix(1,1);
        LinearMasterSlaveConstraint::VectorType constant_vector(1);
        GlobalPointersVector<NodeType> vector_of_node_global_pointers;
        vector_of_node_global_pointers.reserve(remote_nodes_gps_map.size());
        for(auto rank_gp_pair : remote_nodes_gps_map )
        {
            vector_of_node_global_pointers.push_back(rank_gp_pair.second);
        }

        GlobalPointerCommunicator<Node<3>> pointer_comm(r_data_communicator, vector_of_node_global_pointers);
        auto result_proxy = pointer_comm.Apply(
                [](GlobalPointer<Node<3>>& gp){
                    std::vector<DofIdVarKeyPairType> vec_to_return;
                    for(const auto& dof : gp->GetDofs())
                        vec_to_return.push_back( std::make_pair(dof.Id(), dof.GetVariable().Key()) );
                    return vec_to_return;
                    }
        );

        for(const auto& constraint_info : mMasterSlaveDetailsVector)
        {
            const auto& r_slave_node_id = constraint_info->SlaveNodeId();
            const auto& r_master_node_id = constraint_info->MasterNodeId();
            auto slave_dof_id_var_key_pairs = result_proxy.Get(remote_nodes_gps_map[r_slave_node_id]);
            auto master_dof_id_var_key_pairs = result_proxy.Get(remote_nodes_gps_map[r_master_node_id]);

            auto slave_pair = std::find(std::begin(slave_dof_id_var_key_pairs), std::end(slave_dof_id_var_key_pairs),
                                        std::make_pair(r_slave_node_id, constraint_info->GetVariableKey()));
            global_ptrs_vectors_dof_id_key.push_back(Kratos::make_shared<std::pair<IndexType, IndexType>>(*slave_pair));


            relation_matrix(0,0) = constraint_info->Weight();
            constant_vector(0) = constraint_info->Constant();
        }


        auto remote_dofs_gps_map = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(global_ptrs_vectors_dof_id_key, mNodeIdsToSync, r_data_communicator);
        // std::cout<<"my rank :: "<<current_rank<<"  size :: "<<mMasterSlaveDetailsVector.size()<<std::endl;
    }



    ///@}

}; // Class MpiConstraintsUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // MPI_CONSTRAINTS_UTILITY  defined
