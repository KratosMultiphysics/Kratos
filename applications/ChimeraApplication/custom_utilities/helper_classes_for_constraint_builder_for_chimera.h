//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//
//
#if !defined (KRATOS_HELPER_CLASSES_FOR_CONSTRAINT_BUILDER_FOR_CHIMERA_H_INCLUDED)
#define KRATOS_HELPER_CLASSES_FOR_CONSTRAINT_BUILDER_FOR_CHIMERA_H_INCLUDED

#include "utilities/helper_classes_for_constraint_builder.h"

// System includes
#include <vector>
#include <unordered_map>

// project includes
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/node.h"
#include "includes/lock_object.h"

namespace Kratos
{
namespace Internals
{
///@name Internals Globals
///@{

///@}
///@name Type Definitions
///@{


/// Matrix and vector definition
typedef Kratos::Matrix MatrixType;
typedef Kratos::Vector VectorType;

/// Indexes definition
typedef IndexedObject::IndexType IndexType;
typedef std::vector<IndexType> VectorIndexType;


///@}
///@name Internal Classes
///@{

/**
 * @class ConstraintImposerForChimera
 * @ingroup KratosCore
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          > // Made template to include the possibility to work with both local and global matrices for imposing the constraints.
class ConstraintImposerForChimera 
{
public:
    ///@name Type Definitions
    ///@{
    typedef Internals::AuxiliaryGlobalMasterSlaveConstraint AuxiliaryGlobalMasterSlaveRelationType;
    typedef std::unordered_map< IndexType, unique_ptr< AuxiliaryGlobalMasterSlaveRelationType > > GlobalMasterSlaveRelationContainerType;
    typedef std::vector<Dof<double>::Pointer> DofsVectorType;
    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;
    typedef Internals::LocalIndices LocalIndicesType;
    typedef Kratos::Matrix MatrixType;
    typedef Kratos::Vector VectorType;
    typedef std::vector<IndexType> VectorIndexType;
    typedef std::vector<IndexType> EquationIdVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    explicit ConstraintImposerForChimera(GlobalMasterSlaveRelationContainerType& rGlobalMasterSlaveRelations)
                    : mrGlobalMasterSlaveConstraints(rGlobalMasterSlaveRelations)
    {
    }

    ~ConstraintImposerForChimera()
    {
    }

    ConstraintImposerForChimera( const ConstraintImposerForChimera &OtherObject) 
                        : mrGlobalMasterSlaveConstraints (OtherObject.mrGlobalMasterSlaveConstraints) // copy constructor
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief   This adds the equation IDs of masters of all the slaves corresponding to pCurrentElement to EquationIds
     * @details Here cannot use the pure Geometry because, we would need the dof list from the element/geometry.
     * @param   rCurrentContainer the element or condition where the rEquationIds to be modified for master-slave constraints
     * @param   rEquationIds the equation id vector for the above element or condition
     * @param   rCurrentProcessInfo the current process info
     */
    template <typename TContainerType>
    void ApplyConstraints(TContainerType& rCurrentContainer,
                          typename TContainerType::EquationIdVectorType& rEquationIds,
                          ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        this->Reset();
        // If no slave is found for this container , no need of going on
        if (! Internals::HasSlaveNode(rCurrentContainer.GetGeometry()))
        {
            return;
        }
        DofsVectorType ContainerDofs;
        rCurrentContainer.GetDofList(ContainerDofs, rCurrentProcessInfo);
        IndexType slave_equation_id;
        // For each node check if it is ac slave or not If it is .. we change the Transformation matrix
        for (IndexType j = 0; j < ContainerDofs.size(); j++)
        {
            slave_equation_id = ContainerDofs[j]->EquationId(); // consider everything as a slave.
            // Get the global constraint equation for this slave.
            auto global_master_slave_constraint = mrGlobalMasterSlaveConstraints.find(slave_equation_id);
            if (global_master_slave_constraint != mrGlobalMasterSlaveConstraints.end())
            { // if a equation exists for this slave
                global_master_slave_constraint->second->EquationIdsVector(slave_equation_id, mMasterEquationIds); // get the slave and master equation ids for this slave.
                rEquationIds.reserve(mMasterEquationIds.size());
                for (auto &master_eq_id : mMasterEquationIds)
                {
                    // Add the current slaves master eq ids to the equation ids
                    rEquationIds.push_back(master_eq_id);
                }
            }
        }
        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints::ApplyConstraints failed ..");
    }

    /**
     * @brief   This function modifies the LHS and RHS of the rCurrentContainer to account for any master-slave constraints its nodes/dofs
     *          are carrying.
     * @details Here cannot use the pure Geometry because, we would need the dof list from the element/geometry.
     * @param   rCurrentContainer the element or condition where the rEquationIds to be modified for master-slave constraints
     * @param   rLHSContribution the LHS contribution of the rCurrentContainer
     * @param   rRHSContribution the RHS contribution of the rCurrentContainer
     * @param   rEquationIds the equation id vector for the above element or condition
     * @param   rCurrentProcessInfo the current process info
     */
    template <typename TContainerType>
    void ApplyConstraints(TContainerType& rCurrentContainer,
                          LocalSystemMatrixType& rLHSContribution,
                          LocalSystemVectorType& rRHSContribution,
                          typename TContainerType::EquationIdVectorType& rEquationIds,
                          ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        // If no slave is found for this container , no need of going on
        if (! Internals::HasSlaveNode(rCurrentContainer.GetGeometry()))
            return;
        this->Reset();
        // Saving th original system size
        const IndexType initial_sys_size = rLHSContribution.size1();

        // first fill in the rEquationIds using the above function (overloaded one)
        ApplyConstraints<TContainerType>(rCurrentContainer, rEquationIds, rCurrentProcessInfo); // now rEquationIds has all the slave equation ids appended to it.
        IndexType total_number_of_masters = rEquationIds.size() - initial_sys_size;
        // Calculating the local indices corresponding to internal, master, slave dofs of this container
        CalculateLocalIndices(rEquationIds, mLocalIndices, total_number_of_masters);
        // resizing the matrices to the new required length
        ResizeAndInitializeLocalMatrices(rLHSContribution, rRHSContribution, rEquationIds.size());

        // Calculating the F = T'*(F-K*g) which is local to this container
        ModifyRHSForConstraints(rLHSContribution, rRHSContribution, rEquationIds);
        // Calculating the K = T' * K *T which is local to this container
        ModifyLHSForConstraints(rLHSContribution, rRHSContribution, rEquationIds);

        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints:: Applying Multipoint constraints failed ..");
    }
    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    GlobalMasterSlaveRelationContainerType& mrGlobalMasterSlaveConstraints;
    // For Formulating which are the internal, slave indices locally.
    LocalIndicesType mLocalIndices;
    // container's transformation matrix and constant vector
    MatrixType mTransformationMatrixLocal;
    VectorType mConstantVectorLocal;
    // containers for holding equation ids and container dofs
    EquationIdVectorType mMasterEquationIds;
    DofsVectorType mContainerDofs;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief   This function does two operations : K = T' * K * T . This operations are done in place.
     *          Meaning that there is no memory duplication and no explicit matrix and matrix or matrix vector multiplication.
     *          Individual entries of K and F are modified to achieve the result.
     * @param   rLHSContribution The lhs matrix of the container
     * @param   rRHSContribution The rhs vector of the container
     * @param   rEquationIds the list of equation ids (extended with the masters).
     */
    void ModifyLHSForConstraints(MatrixType &rLHSContribution, VectorType& rRHSContribution, EquationIdVectorType &rEquationIds)
    {
        mLocalIndices.container_master_weights.reserve(mLocalIndices.master_index_vector.size());
        mLocalIndices.container_master_slaves.reserve(mLocalIndices.master_index_vector.size());
        mLocalIndices.processed_master_indices.reserve(mLocalIndices.master_index_vector.size());
        IndexType slave_equation_id;
        EquationIdVectorType master_equation_ids;
        VectorType master_weights_vector;
        double slave_constant;

        for (auto& slave_index : mLocalIndices.slave_index_vector) { // Loop over all the slaves for this container
            // Get the global equation for this constraint
            auto global_master_slave_constraint = mrGlobalMasterSlaveConstraints.find(rEquationIds[slave_index]);
            // Get the tranformation matrix and constant_vector from the current slave
            global_master_slave_constraint->second->EquationIdsVector(slave_equation_id, master_equation_ids);
            global_master_slave_constraint->second->CalculateLocalSystem(master_weights_vector, slave_constant);

            IndexType master_index = 0;
            double master_weight = 0.0;
            IndexType i_master = 0;
            for (auto&  master_eq_id : master_equation_ids)
            { // Loop over all the masters the slave has
                master_index = std::distance(rEquationIds.begin(), std::find(rEquationIds.begin(), rEquationIds.end(), master_eq_id));
                //master_weight = mTransformationMatrixLocal(slave_index,master_index);
                master_weight = master_weights_vector(i_master);
                for (auto& internal_index : mLocalIndices.internal_index_vector) {
                    // For K(m,u) and K(u,m)
                    rLHSContribution(internal_index, master_index) += rLHSContribution(internal_index, slave_index) * master_weight;
                    //rLHSContribution(master_index, internal_index) += rLHSContribution(slave_index, internal_index) * master_weight;
                }

                mLocalIndices.container_master_weights.push_back( master_weight );
                mLocalIndices.container_master_slaves.push_back( slave_index );
                mLocalIndices.processed_master_indices.push_back( master_index );
                i_master++;
            } // Loop over all the masters the slave has
        }

        //Adding contribution from slave to Kmm
        IndexType master_i = 0;
        for (auto& master_index : mLocalIndices.processed_master_indices) {
            IndexType master_i_other = 0;
            for (auto& master_index_other : mLocalIndices.processed_master_indices) {
                //rLHSContribution(master_index, master_index_other) += mLocalIndices.container_master_weights[master_i] *
                                                                        rLHSContribution(mLocalIndices.container_master_slaves[master_i], mLocalIndices.container_master_slaves[master_i_other])
                                                                        * mLocalIndices.container_master_weights[master_i_other];
                master_i_other++;
            }
            master_i++;
        }

        // For K(u,s) and K(s,u). This is to be done at the end only
        for (auto& slave_index : mLocalIndices.slave_index_vector) {
            for (auto& internal_index : mLocalIndices.internal_index_vector) {
                rLHSContribution(slave_index, internal_index) = 0.0;
                rLHSContribution(internal_index, slave_index) = 0.0;
            }
        }
    }

    /**
     * @brief   This function does two operation : F = T'*(F-K*b). This operation is done in place.
     *          Meaning that there is no memory duplication and no explicit matrix and matrix or matrix vector multiplication.
     *          Individual entries of K and F are modified to achieve the result. 
     * @param   rLHSContribution The lhs matrix of the container
     * @param   rRHSContribution The rhs vector of the container
     * @param   rEquationIds the list of equation ids (extended with the masters).
     */
    void ModifyRHSForConstraints(MatrixType &rLHSContribution, VectorType& rRHSContribution, EquationIdVectorType &rEquationIds)
    {
        IndexType slave_equation_id;
        EquationIdVectorType master_equation_ids;
        VectorType master_weights_vector;
        double slave_constant;
        VectorType master_weights_vector_other;
        double constant_other;

        for (auto& slave_index : mLocalIndices.slave_index_vector) { // Loop over all the slaves for this container
            // Get the global equation for this constraint
            auto global_master_slave_constraint = mrGlobalMasterSlaveConstraints.find(rEquationIds[slave_index]);
            // Get the tranformation matrix and constant_vector from the current slave
            global_master_slave_constraint->second->EquationIdsVector(slave_equation_id, master_equation_ids);
            global_master_slave_constraint->second->CalculateLocalSystem(master_weights_vector, slave_constant);

            IndexType master_index = 0;
            double master_weight = 0.0;
            IndexType i_master = 0;
            for (auto&  master_eq_id : master_equation_ids)
            { // Loop over all the masters the slave has
                master_index = std::distance(rEquationIds.begin(), std::find(rEquationIds.begin(), rEquationIds.end(), master_eq_id));
                //master_weight = mTransformationMatrixLocal(slave_index,master_index);
                master_weight = master_weights_vector(i_master);
                for (auto& internal_index : mLocalIndices.internal_index_vector) {
                    rRHSContribution(internal_index) -= rLHSContribution(internal_index, slave_index) * slave_constant;
                }
                // For RHS(m) += A'*LHS(s,s)*B
                for (auto& slave_index_other : mLocalIndices.slave_index_vector) {
                    auto global_master_slave_constraint_other = mrGlobalMasterSlaveConstraints.find(rEquationIds[slave_index_other]);
                    global_master_slave_constraint_other->second->CalculateLocalSystem(master_weights_vector_other, constant_other);
                    //rRHSContribution(master_index) -= rLHSContribution(slave_index, slave_index_other) * master_weight * constant_other;
                }
                // Changing the RHS side of the equation
                //rRHSContribution(master_index) += master_weight * rRHSContribution(slave_index);

                i_master++;
            } // Loop over all the masters the slave has

            rRHSContribution(slave_index) = 0.0;
        }
    }




    /**
     * @brief   Resets the member vectors and matrices to zero and zero size
     */
    void Reset()
    {
        mLocalIndices.Reset();
        mTransformationMatrixLocal.resize(0,0, false);
        mConstantVectorLocal.resize(0, false);
        mMasterEquationIds.clear();
        mContainerDofs.clear();
    }

    /**
     * @brief   This function calculates the local indices of a given element or condition
     * @param   rEquationIds vector of the equation ids
     * @param   rLocalIndexStructure reference to the structure of LocalIndicesType
     */
    void CalculateLocalIndices(EquationIdVectorType& rEquationIds, LocalIndicesType& rLocalIndexStructure, IndexType rTotalNumberOfMasters)
    {
        CalculateLocalSlaveIndices(rEquationIds, rLocalIndexStructure);
        CalculateLocalInternalIndices(rEquationIds, rLocalIndexStructure);
        CalculateLocalMasterIndices(rEquationIds, rLocalIndexStructure, rTotalNumberOfMasters);
    }



    /**
     * @brief   This function calculates the local slave indices of a given element or condition
     * @param   rEquationIds vector of the equation ids
     * @param   rLocalSlaveIndexVector reference to the vector of slave indices
     */
    void CalculateLocalSlaveIndices(EquationIdVectorType& rEquationIds, LocalIndicesType& rLocalIndexStructure)
    {
        KRATOS_TRY
        int index = 0;
        for (auto &eq_id : rEquationIds)
        {
            auto global_master_slave_constraint = mrGlobalMasterSlaveConstraints.find(eq_id);
            if (global_master_slave_constraint != mrGlobalMasterSlaveConstraints.end())
                rLocalIndexStructure.slave_index_vector.push_back(index);

            index++;
        }
        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints::CalculateLocalSlaveIndices failed ..");
    }

    /**
     * @brief   This function calculates the local internal indices of a given element or condition
     * @param   rEquationIds vector of the equation ids
     * @param   rLocalIndexStructure reference to the vector of slave indices
     */
    void CalculateLocalInternalIndices(EquationIdVectorType& rEquationIds, LocalIndicesType& rLocalIndexStructure)
    {
        KRATOS_TRY
        VectorIndexType local_index_vector(rEquationIds.size());
        for (IndexType i = 0; i<rEquationIds.size(); ++i)
            local_index_vector[i] = i;

        std::sort(local_index_vector.begin(), local_index_vector.end());
        std::sort(rLocalIndexStructure.slave_index_vector.begin(), rLocalIndexStructure.slave_index_vector.end());

        std::set_difference(local_index_vector.begin(), local_index_vector.end(),
                            rLocalIndexStructure.slave_index_vector.begin(), rLocalIndexStructure.slave_index_vector.end(),
                            std::back_inserter(rLocalIndexStructure.internal_index_vector));

        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints::CalculateLocalInternalIndices failed ..");
    }

    /**
     * @brief   This function calculates the local internal indices of a given element or condition
     * @param   rEquationIds vector of the equation ids
     * @param   rLocalIndexStructure reference to the vector of slave indices
     * @param   rTotalNumberOfMasters total number of masters for the given element or condition.
     */
    void CalculateLocalMasterIndices(EquationIdVectorType& rEquationIds, LocalIndicesType& rLocalIndexStructure, IndexType rTotalNumberOfMasters)
    {
        // Get number of master indices for this current container
        rLocalIndexStructure.master_index_vector.reserve(rTotalNumberOfMasters + rEquationIds.size() );
        for (IndexType i = rEquationIds.size()-1; i < rEquationIds.size() -rTotalNumberOfMasters; --i)
            rLocalIndexStructure.master_index_vector.push_back(i);
    }

    ///@}
};

} // namespace Internals

} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
