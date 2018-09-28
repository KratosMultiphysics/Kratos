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

#if !defined(AUXILIARY_GLOBAL_MASTER_SLAVE_RELATION)
#define AUXILIARY_GLOBAL_MASTER_SLAVE_RELATION
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

/// Geometric definitions
typedef Node<3> NodeType;
typedef Geometry<NodeType> GeometryType;

/// Matrix and vector definition
typedef Kratos::Matrix MatrixType;
typedef Kratos::Vector VectorType;

/// Indexes definition
typedef IndexedObject::IndexType IndexType;
typedef std::vector<IndexType> VectorIndexType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

/**
 * @brief this method checks if any of the nodes of the given rGeometry is marked SLAVE.
 * @param rGeometry The geometry to check for.
 */
bool HasSlaveNode(GeometryType& rGeometry)
{
    for(auto& node : rGeometry)
        if (node.IsDefined(SLAVE))
            return node.Is(SLAVE);

    return false;
}

/**
 * @brief   This function resizes the given matrix and vector pair to the new size provided.
 *          And Initializes the extra part added to zero.
 * @param   rMatrix matrix to be resized
 * @param   rVector vector to be resized
 * @param   FinalSize the final size of the resized quantities.
 */
void ResizeAndInitializeLocalMatrices(MatrixType& rMatrix, VectorType& rVector,
                                        IndexType FinalSize)
{
    KRATOS_TRY
    // storing the initial matrix and vector and their properties
    KRATOS_ERROR_IF(rMatrix.size1() != rVector.size())<<"ResizeAndInitializeLocalMatrices :: Dimension of the matrix and vector passed are not the same !"<<std::endl;
    const IndexType initial_sys_size = rMatrix.size1();
    MatrixType matrix(initial_sys_size, initial_sys_size);
    noalias(matrix) = rMatrix;
    VectorType vector(initial_sys_size);
    noalias(vector) = rVector;

    rMatrix.resize(FinalSize, FinalSize, false);
    rVector.resize(FinalSize, false);
    // reassigning the original part of the matrix
    for (IndexType m = 0; m < initial_sys_size; ++m)
    {
        for (IndexType n = 0; n < initial_sys_size; ++n)
        {
            rMatrix(m,n) = matrix(m,n);
        }
        rVector(m) = vector(m);
    }
    // Making the extra part of matrix zero
    for (IndexType m = initial_sys_size; m < FinalSize; ++m)
    {
        for (IndexType n = 0; n < FinalSize; ++n)
        {
            rMatrix(m, n) = 0.0;
            rMatrix(n, m) = 0.0;
        }
        rVector(m) = 0.0;
    }
    KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints::ResizeAndInitializeLocalMatrices failed ..");
}

///@}
///@name Internals Classes
///@{

/**
 * @class AuxiliaryGlobalMasterSlaveConstraint
 * @ingroup KratosCore
 * @brief This class stores the information regarding the AuxiliaryGlobalMasterSlaveConstraint equation.
 *         Naming convention is defined like this. (each object of this class will store one equation in the given form
 *
 *   SlaveEquationId = w_1*MasterEquationId_1 + w_2*MasterEquationId_2 + ..... + w_n*MasterEquationId_n
 *
 *   This stores the condensed form of the MasterSlaveConstraint objects into one object. if only one relation for a slave is added as
 *   MasterSlaveConstraint then there will only be one entry for master for its corresponding AuxiliaryGlobalMasterSlaveConstraint.
 *   Currently this class is designed to hold only one equation. There is only one unique object of this class for each slave.
 *
 *   Future plan is to also make it possible to work with matrices (T) and vectors (for slave and master equation ids and constants)
 *
 *
 *  IMPORTANT : This is not seen by the user. This is a helper data structure which is exists only in the builder and solver.
 *
 * @author Aditya Ghantasala
 */
class AuxiliaryGlobalMasterSlaveConstraint : public IndexedObject
{
public:
    ///@name Type Definitions
    ///@{

    typedef IndexedObject BaseType;
    typedef Internals::IndexType IndexType;
    typedef Internals::MatrixType MatrixType;
    typedef Internals::VectorType VectorType;
    typedef std::vector<IndexType> EquationIdVectorType;

    /// Pointer definition of AuxiliaryGlobalMasterSlaveConstraint
    KRATOS_CLASS_POINTER_DEFINITION(AuxiliaryGlobalMasterSlaveConstraint);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor of the class
     * @param SlaveEquationId the slave equation id for which this class is being constructed.
     */
    explicit AuxiliaryGlobalMasterSlaveConstraint(IndexType SlaveEquationId = 0) : IndexedObject(SlaveEquationId),
                                                                                    mLhsValue(0.0),
                                                                                    mRhsValue(0.0)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to get the slave equation Id corresponding to this constraint.
     * @param Constant the value of the constant to be assigned.
     */
    IndexType SlaveEquationId() const { return this->Id(); }

    /**
     * @brief Function to set the lefthand side of the constraint (the slave dof value)
     * @param LhsValue the value of the lhs (the slave dof value)
     */
    void SetLeftHandSide(const double LhsValue)
    {
        mLockObject.SetLock();
        mLhsValue = LhsValue;
        mLockObject.UnSetLock();
    }

    /**
     * @brief Function to update the righthand side of the constraint (the combination of all the master dof values and constants)
     * @param RHSValue the value of the lhs (the slave dof value)
     */
    void SetRightHandSide(const double RhsValue)
    {
        mRhsValue = RhsValue;
    }
    void UpdateRightHandSide(const double RhsValueUpdate)
    {
        mLockObject.SetLock();
        mRhsValue = mRhsValue + RhsValueUpdate;
        mLockObject.UnSetLock();
    }

    // Get number of masters for this slave
    IndexType NumberOfMasters() const
    {
        return mMasterEquationIdVector.size();
    }

    /**
     * @brief this determines the master equation IDs connected to this constraint
     * @param rResult the elemental equation ID vector
     */
    virtual void EquationIdsVector(IndexType& rSlaveEquationId,
                                  EquationIdVectorType& rMasterEquationIds)
    {
        if (rMasterEquationIds.size() == 0)
            rMasterEquationIds.resize(this->NumberOfMasters(), false);

        rSlaveEquationId = this->SlaveEquationId();
        rMasterEquationIds = mMasterEquationIdVector;
    }

    /**
     * @brief   this is called during the assembling process in order
     *          to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rMasterWeightsVector the elemental left hand side matrix
     * @param rConstant the elemental right hand side
     */
    virtual void CalculateLocalSystem(VectorType &rMasterWeightsVector,
                                      double &rConstant)
    {
        if (rMasterWeightsVector.size() == 0)
            rMasterWeightsVector.resize(this->NumberOfMasters(), false);

        for (IndexType i = 0; i < this->NumberOfMasters(); ++i)
            rMasterWeightsVector(i) = mMasterWeightsVector[i];

        /// Here this is required because, when in the builder and solver , we are actually imposing the constraint on the update
        /// of the DOF value (residual formulation), this does not necessarily guarantee the DOFs themselves follow the constraint equation.
        /// So, we calculate the LHS value and RHS value of the constraint equation (with DOF values) and if they are not
        /// satisfying the constraint, we use the residual as the constant.
        rConstant = mRhsValue - mLhsValue;

    }

    /**
     * @brief This method clears the equations ids
     */
    void Clear()
    {
        //clearing the contents
        mMasterEquationIdVector.clear();
        mMasterWeightsVector.clear();
        //shrinking the memory
        mMasterEquationIdVector.shrink_to_fit();
        mMasterWeightsVector.shrink_to_fit();
    }

    /**
     * @brief This method adds a new master
     */
    void AddMaster(const IndexType MasterEquationId, const double Weight)
    {
        const int index = GetMasterEquationIdPosition(MasterEquationId);
        if (index >= 0) {
            #pragma omp atomic
            mMasterWeightsVector[index] += Weight;
        } else {
            mLockObject.SetLock(); // locking for exclusive access to the vectors mMasterEquationIdVector and mMasterWeightsVectors
            mMasterEquationIdVector.push_back(MasterEquationId);
            mMasterWeightsVector.push_back(Weight);
            mLockObject.UnSetLock(); // unlocking
        }
    }

    /**
     * @brief This method resers the LHS/RHS relationship
     */
    void Reset()
    {
        this->mLhsValue = 0.0;
        this->mRhsValue = 0.0;
    }

    /**
     * @brief This method returns the correspondin EquationId for the master
     */
    int GetMasterEquationIdPosition(const IndexType MasterEquationId) const
    {
        auto it = find(mMasterEquationIdVector.begin(), mMasterEquationIdVector.end(), MasterEquationId);
        if (it != mMasterEquationIdVector.end())
            return it - mMasterEquationIdVector.begin();
        else
            return -1;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "AuxiliaryGlobalMasterSlaveConstraint # " << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    double mLhsValue;
    double mRhsValue;

    EquationIdVectorType mMasterEquationIdVector;
    std::vector<double> mMasterWeightsVector;

    LockObject mLockObject;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer &rSerializer) const override
    {
        // No need to save anything from this class as they will be reconstructed
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject);
    }

    void load(Serializer &rSerializer) override
    {
        // No need to load anything from this class as they will be reconstructed
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject);
    }

    ///@}
}; // End of ConstraintEquation class

/**
 * @struct LocalIndices
 * @ingroup KratosCore
 * @brief This class stores the stores three different vectors of local internal, slave, master indices
 *          which are used in constraint builder and solver.
 *
 * @author Aditya Ghantasala
 */
struct LocalIndices
{
    typedef Internals::IndexType IndexType;
    typedef Internals::VectorIndexType VectorIndexType;

    void Reset()
    {
        internal_index_vector.resize(0);
        master_index_vector.resize(0);
        slave_index_vector.resize(0);
    }

    VectorIndexType internal_index_vector; // indicies corresponding to internal DOFs
    VectorIndexType master_index_vector; // indicies corresponding to master DOFs
    VectorIndexType slave_index_vector; // indicies corresponding to slave DOFs
};

///@}
///@name Type Definitions
///@{

/// AuxiliaryGlobalMasterSlaveConstraint definitions
typedef Internals::AuxiliaryGlobalMasterSlaveConstraint AuxiliaryGlobalMasterSlaveConstraintType;
//typedef PointerVectorSet<AuxiliaryGlobalMasterSlaveConstraint, IndexedObject> GlobalMasterSlaveRelationContainerType;
typedef std::unordered_map< IndexType, unique_ptr< AuxiliaryGlobalMasterSlaveConstraintType > > GlobalMasterSlaveRelationContainerType;

///@}
///@name Internal Classes
///@{

/**
 * @class ConstraintImposer
 * @ingroup KratosCore
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          > // Made template to include the possibility to work with both local and global matrices for imposing the constraints.
class ConstraintImposer {
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

    explicit ConstraintImposer(GlobalMasterSlaveRelationContainerType& rGlobalMasterSlaveRelations)
        : mrGlobalMasterSlaveConstraints(rGlobalMasterSlaveRelations)
    {
    }

    ~ConstraintImposer()
    {
    }

    ConstraintImposer( const ConstraintImposer &OtherObject) :
                mrGlobalMasterSlaveConstraints (OtherObject.mrGlobalMasterSlaveConstraints) // copy constructor
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
        typename TContainerType::EquationIdVectorType equation_ids = rEquationIds;
        // Saving th original system size
        const IndexType initial_sys_size = rLHSContribution.size1();

        // first fill in the rEquationIds using the above function (overloaded one)
        ApplyConstraints<TContainerType>(rCurrentContainer, rEquationIds, rCurrentProcessInfo); // now rEquationIds has all the slave equation ids appended to it.
        IndexType total_number_of_masters = rEquationIds.size() - initial_sys_size;
        // Calculating the local indices corresponding to internal, master, slave dofs of this container
        CalculateLocalIndices(rEquationIds, mLocalIndices, total_number_of_masters);

        // resizing the matrices to the new required length
        ResizeAndInitializeLocalMatrices(rLHSContribution, rRHSContribution, rEquationIds.size());

        // Calculating the T and C which are local to this container
        ModifyMatrices(rLHSContribution, rRHSContribution, rEquationIds);

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
     * @brief   This function does two operations : K = T' * K * T and F = T'*(F-K*b). Both these operations are done in place.
     *          Meaning that there is no memory duplication and no explicit matrix and matrix or matrix vector multiplication.
     *          Individual entries of K and F are modified to achieve the result. 
     * @param   rLHSContribution The lhs matrix of the container
     * @param   rRHSContribution The rhs vector of the container
     * @param   rEquationIds the list of equation ids (extended with the masters).
     */
    void ModifyMatrices(MatrixType &rLHSContribution, VectorType& rRHSContribution, EquationIdVectorType &rEquationIds)
    {
        std::vector<double> container_master_weights;
        container_master_weights.reserve(mLocalIndices.master_index_vector.size());
        std::vector<IndexType> container_master_slaves;
        container_master_slaves.reserve(mLocalIndices.master_index_vector.size());
        std::vector<IndexType> processed_master_indices;
        processed_master_indices.reserve(mLocalIndices.master_index_vector.size());
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
                    rRHSContribution(internal_index) += -rLHSContribution(internal_index, slave_index) * slave_constant;
                    // For K(m,u) and K(u,m)
                    rLHSContribution(internal_index, master_index) += rLHSContribution(internal_index, slave_index) * master_weight;
                    rLHSContribution(master_index, internal_index) += rLHSContribution(slave_index, internal_index) * master_weight;
                }
                // For RHS(m) += A'*LHS(s,s)*B
                for (auto& slave_index_other : mLocalIndices.slave_index_vector) {
                    auto global_master_slave_constraint_other = mrGlobalMasterSlaveConstraints.find(rEquationIds[slave_index_other]);
                    global_master_slave_constraint_other->second->CalculateLocalSystem(master_weights_vector_other, constant_other);
                    rRHSContribution(master_index) -= rLHSContribution(slave_index, slave_index_other) * master_weight * constant_other;
                }
                // Changing the RHS side of the equation
                rRHSContribution(master_index) += master_weight * rRHSContribution(slave_index);

                container_master_weights.push_back( master_weight );
                container_master_slaves.push_back( slave_index );
                processed_master_indices.push_back( master_index );
                i_master++;
            } // Loop over all the masters the slave has

            rRHSContribution(slave_index) = 0.0;
        }

        //Adding contribution from slave to Kmm
        IndexType master_i = 0;
        for (auto& master_index : processed_master_indices) {
            IndexType master_i_other = 0;
            for (auto& master_index_other : processed_master_indices) {
                rLHSContribution(master_index, master_index_other) += container_master_weights[master_i] *
                                                                        rLHSContribution(container_master_slaves[master_i], container_master_slaves[master_i_other])
                                                                        * container_master_weights[master_i_other];
                master_i_other++;
            }
            master_i++;
        }

        // For K(u,s) and K(s,u)
        for (auto& slave_index : mLocalIndices.slave_index_vector) {
            for (auto internal_index : mLocalIndices.internal_index_vector) {
                rLHSContribution(slave_index, internal_index) = 0.0;
                rLHSContribution(internal_index, slave_index) = 0.0;
            }
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
     * @brief   This function calculates the local transformation matrix and the constant vector for each
     *          each element or condition. The T matrix and C vector for each element or condition for the slaves they contain .
     * @param   rLocalIndices object of Struct LocalIndicesType containing the local internal, master and slave indices
     * @param   rTransformationMatrixLocal reference to the tranformation matrix which is to be calculated.
     * @param   rEquationIds the list of equation ids.
     */
    void CalculateLocalTransformationMatrix(LocalIndicesType& rLocalIndices,
                                            MatrixType& rTransformationMatrixLocal,
                                            EquationIdVectorType& rEquationIds)
    {
        KRATOS_TRY
        IndexType slave_equation_id;
        EquationIdVectorType master_equation_ids;
        VectorType mMasterWeightsVector;
        double slave_constant;
        int i_masters_total = rEquationIds.size();
        for (const auto &slave_index : rLocalIndices.slave_index_vector)
        {
            auto global_master_slave_constraint = mrGlobalMasterSlaveConstraints.find(rEquationIds[slave_index]);
            KRATOS_DEBUG_ERROR_IF (global_master_slave_constraint == mrGlobalMasterSlaveConstraints.end()) <<
                             "No master slave constraint equation found for atleast one of the dofs .. !" << std::endl;
            global_master_slave_constraint->second->EquationIdsVector(slave_equation_id, master_equation_ids);
            global_master_slave_constraint->second->CalculateLocalSystem(mMasterWeightsVector, slave_constant);
            for (IndexType i_master = 0; i_master < master_equation_ids.size(); ++i_master)
            {
                rTransformationMatrixLocal(slave_index, i_masters_total) += mMasterWeightsVector(i_master);
                i_masters_total++;
            }
        }
        for (const auto &master_index : rLocalIndices.master_index_vector)
            rTransformationMatrixLocal(master_index, master_index) = 1.0;

        for (const auto &internal_index : rLocalIndices.internal_index_vector)
            rTransformationMatrixLocal(internal_index, internal_index) = 1.0;

        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints::CalculateLocalTransformationMatrix failed ..");
    }

    /**
     * @brief   This function calculates the local constant vector for each
     *          each element or condition. C vector for each element or condition for the slaves they contain .
     * @param   rLocalSlaveIndexVector vector of slave indices
     * @param   rConstantVectorLocal reference to the constant vector to be calculated
     * @param   rEquationIds the list of equation ids.
     */
    void CalculateLocalConstantVector(LocalIndicesType& rLocalIndexStructure,
                                      VectorType& rConstantVectorLocal,
                                      EquationIdVectorType& rEquationIds)
    {
        KRATOS_TRY
        VectorType mMasterWeightsVector;
        double slave_constant;

        for (const auto &slave_index : rLocalIndexStructure.slave_index_vector)
        {
            auto global_master_slave_constraint = mrGlobalMasterSlaveConstraints.find(rEquationIds[slave_index]);
            if (global_master_slave_constraint != mrGlobalMasterSlaveConstraints.end())
            {
                global_master_slave_constraint->second->CalculateLocalSystem(mMasterWeightsVector, slave_constant);
                rConstantVectorLocal(slave_index) = slave_constant;
            }
            else
                KRATOS_ERROR << "No master slave constraint equation found for atleast one of the dofs .. !" << std::endl;
        }
        KRATOS_CATCH("ResidualBasedBlockBuilderAndSolverWithConstraints::CalculateLocalConstantVector failed ..");
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
