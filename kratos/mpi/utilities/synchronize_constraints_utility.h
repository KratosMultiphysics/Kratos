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
#include <typeinfo>


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
    typedef GlobalPointersVector<DofType> GlobalPointerVectorDofsType;

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
        virtual IndexType Id(){return mId;}
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
    GlobalPointerVectorDofsType mVectorOfGlobalPointers;

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
        auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();
        //IndexType current_rank = r_data_communicator.Rank();
        auto remote_nodes_gps_map = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(mrModelPart.Nodes(), mNodeIdsToSync, r_data_communicator);
        GlobalPointersVector<NodeType> vector_of_node_global_pointers;
        vector_of_node_global_pointers.reserve(remote_nodes_gps_map.size());
        for(auto rank_gp_pair : remote_nodes_gps_map )
        {
            vector_of_node_global_pointers.push_back(rank_gp_pair.second);
        }

        GlobalPointerCommunicator<Node<3>> nodes_pointer_comm(r_data_communicator, vector_of_node_global_pointers);
        auto nodes_result_proxy = nodes_pointer_comm.Apply(
                [](GlobalPointer<Node<3>>& gp){
                    typedef std::pair<std::size_t, GlobalPointer<DofType>> VarKeyDofGpPairType;
                    std::vector<VarKeyDofGpPairType> vec_of_dofs_global_ptr;
                    for (auto& dof : gp->GetDofs())
                    {
                            auto dof_shr_ptr = gp->pGetDof(dof.GetVariable());
                            vec_of_dofs_global_ptr.push_back(
                                std::make_pair(dof.GetVariable().Key(), GlobalPointer<DofType>(dof_shr_ptr, gp->FastGetSolutionStepValue(PARTITION_INDEX)) )
                                );
                    }
                    return vec_of_dofs_global_ptr;
                }
        );
        for(const auto& constraint_info : mMasterSlaveDetailsVector)
        {
            LinearMasterSlaveConstraint::DofPointerVectorType slave_dof_vector;
            LinearMasterSlaveConstraint::DofPointerVectorType master_dof_vector;
            LinearMasterSlaveConstraint::MatrixType relation_matrix(1,1);
            LinearMasterSlaveConstraint::VectorType constant_vector(1);
            const auto& r_slave_node_id = constraint_info->SlaveNodeId();
            const auto& r_master_node_id = constraint_info->MasterNodeId();
            auto slave_node_var_dofs_gps_map = nodes_result_proxy.Get(remote_nodes_gps_map[r_slave_node_id]);
            auto master_node_var_dofs_gps_map = nodes_result_proxy.Get(remote_nodes_gps_map[r_master_node_id]);
            // Get the global pointers for the dofs of the slave and masters -> construct them using the PARTITION_INDEX of the node
            for(auto& slave_gp_pair : slave_node_var_dofs_gps_map)
            {
                if(slave_gp_pair.first == constraint_info->GetVariableKey()){
                    mVectorOfGlobalPointers.push_back(slave_gp_pair.second);
                    slave_dof_vector.push_back( DofType::Pointer( (*(mVectorOfGlobalPointers.ptr_end()-1)).get() ) );
                    break;
                }
            }
            for(auto& master_gp_pair : master_node_var_dofs_gps_map)
            {
                if(master_gp_pair.first == constraint_info->GetVariableKey()){
                    mVectorOfGlobalPointers.push_back(master_gp_pair.second);
                    //master_dof_vector.push_back( DofType::Pointer( (*(mVectorOfGlobalPointers.ptr_end()-1)).get() ) );
                    break;
                }
            }
            // Then use the communicator to make the constraints.
            relation_matrix(0,0) = constraint_info->Weight();
            constant_vector(0) = constraint_info->Constant();

            mrModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_info->Id(), master_dof_vector, slave_dof_vector, relation_matrix, constant_vector);
        }
        std::cout<<"Num constraints :: "<<mVectorOfGlobalPointers.size()<<std::endl;
        // // Use the vector of the global pointers to make the pointer_communicator
        // GlobalPointerCommunicator<DofType> dofs_pointer_comm(r_data_communicator, mVectorOfGlobalPointers);
        // auto dofs_result_proxy = dofs_pointer_comm.Apply(
        //         [](GlobalPointer<DofType>& gp){
        //             return gp;
        //         }
        // );

    }

    void GetDofsVector()
    {

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
