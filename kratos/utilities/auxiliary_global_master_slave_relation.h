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
// project includes
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/node.h"
#include "includes/lock_object.h"


namespace Kratos
{

/**
 * @class AuxiliaryGlobalMasterSlaveRelation
 * @ingroup KratosCore
 * @brief This class stores the information regarding the AuxiliaryGlobalMasterSlaveRelation equation.
 *         Naming convention is defined like this. (each object of this class will store one equation in the given form
 *
 *   SlaveEquationId = w_1*MasterEquationId_1 + w_2*MasterEquationId_2 + ..... + w_n*MasterEquationId_n
 *
 *   This stores the condensed form of the MasterSlaveConstraint objects into one object. if only one relation for a slave is added as
 *   MasterSlaveConstraint then there will only be one entry for master for its corresponding AuxiliaryGlobalMasterSlaveRelation.
 *   Currently this class is designed to hold only one equation. There is only one unique object of this class for each slave. 
 *
 *   Future plan is to also make it possible to work with matrices (T) and vectors (for slave and master equation ids and constants)
 *
 *
 *  IMPORTANT : This is not seen by the user. This is a helper data structure which is exists only in the builder and solver. 
 *
 * @author Aditya Ghantasala
 */
class AuxiliaryGlobalMasterSlaveRelation : public IndexedObject
{
  public:
    typedef IndexedObject BaseType;
    typedef BaseType::IndexType IndexType;
    typedef Matrix MatrixType;
    typedef Vector VectorType;
    typedef std::vector<IndexType> EquationIdVectorType;
    KRATOS_CLASS_POINTER_DEFINITION(AuxiliaryGlobalMasterSlaveRelation);

    /**
     * @brief Constructor of the class
     * @param SlaveEquationId the slave equation id for which this class is being constructed.
     */
    explicit AuxiliaryGlobalMasterSlaveRelation(IndexType SlaveEquationId = 0) : IndexedObject(SlaveEquationId), mConstant(0.0)
    {
        mLhsValue = 0.;
        mRhsValue = 0.;
        mConstant = 0.;
    }

    /**
     * @brief Function to get the slave equation Id corresponding to this constraint.
     * @param Constant the value of the constant to be assigned.
     */
    IndexType SlaveEquationId() const { return this->Id(); }

    /**
     * @brief Function to set the lefthand side of the constraint (the slave dof value)
     * @param LhsValue the value of the lhs (the slave dof value)
     */
    void SetLHSValue(const double LhsValue)
    {
        mLockObject.SetLock();
            mLhsValue = LhsValue;
        mLockObject.UnSetLock();
    }

    /**
     * @brief Function to update the righthand side of the constraint (the combination of all the master dof values and constants)
     * @param RHSValue the value of the lhs (the slave dof value)
     */
    void SetRHSValue(const double RhsValue)
    {
        mRhsValue = RhsValue;
    }
    void UpdateRHSValue(const double RhsValueUpdate)
    {
        mLockObject.SetLock();
            mRhsValue = mRhsValue + RhsValueUpdate;
        mLockObject.UnSetLock();
    }

    // Get number of masters for this slave
    IndexType GetNumberOfMasters() const
    {
        return mMasterEquationIdVector.size();
    }

    /**
     * @brief this determines the master equation IDs connected to this constraint
     * @param rResult the elemental equation ID vector
     */
    virtual void EquationIdVector(IndexType& rSlaveEquationId,
                                  EquationIdVectorType& rMasterEquationIds)
    {
        if (rMasterEquationIds.size() != 0)
            rMasterEquationIds.resize(this->GetNumberOfMasters(), false);

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
        if (rMasterWeightsVector.size() != 0)
            rMasterWeightsVector.resize(this->GetNumberOfMasters(), false);

        for (IndexType i = 0; i < this->GetNumberOfMasters(); ++i)
            rMasterWeightsVector(i) = mMasterWeightsVector[i];


        mLockObject.SetLock();
            mConstant = mRhsValue - mLhsValue;
        mLockObject.UnSetLock();
        rConstant = mConstant;
    }

    void PrintInfo() const
    {
        KRATOS_INFO("GlobalMasterSlaveRelation")<<std::endl;
        KRATOS_INFO("GlobalMasterSlaveRelation") << "------------------------------" << std::endl;
        KRATOS_INFO("GlobalMasterSlaveRelation") << "SlaveEquationId :: " << SlaveEquationId() << std::endl;
        KRATOS_INFO("GlobalMasterSlaveRelation") << "Constant :: " << mConstant << std::endl;
        KRATOS_INFO("GlobalMasterSlaveRelation") << "------------------------------ :: Masters = " << mMasterEquationIdVector.size()<< std::endl;
        for (IndexType i = 0; i< mMasterEquationIdVector.size(); i++)
        {
            KRATOS_INFO("GlobalMasterSlaveRelation") << i << " Master  equation id :: " <<mMasterEquationIdVector[i] << ", weight :: " << mMasterWeightsVector[i] << std::endl;
        }
        KRATOS_INFO("GlobalMasterSlaveRelation") << "------------------------------" << std::endl;
    }

    void Clear()
    {
        mLockObject.SetLock(); // locking for exclusive access to the vectors mMasterEquationIdVector and mMasterWeightsVectors
            //clearing the contents
            mMasterEquationIdVector.clear();
            mMasterWeightsVector.clear();
            //shrinking the memory
            mMasterEquationIdVector.shrink_to_fit();
            mMasterWeightsVector.shrink_to_fit();
        mLockObject.UnSetLock(); // unlocking
    }

    void AddMaster(IndexType MasterEquationId, double Weight)
    {
        mLockObject.SetLock(); // locking for exclusive access to the vectors mMasterEquationIdVector and mMasterWeightsVectors
            int index = MasterEquationIdExists(MasterEquationId);
            if (index > 0)
            {
                mMasterWeightsVector[index] += Weight;
            } else
            {
                mMasterEquationIdVector.push_back(MasterEquationId);
                mMasterWeightsVector.push_back(Weight);
            }
        mLockObject.UnSetLock(); // unlocking
    }

  private:
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

    int MasterEquationIdExists(IndexType MasterEquationId)
    {
        auto it = find(mMasterEquationIdVector.begin(), mMasterEquationIdVector.end(), MasterEquationId);
        if (it != mMasterEquationIdVector.end())
            return it - mMasterEquationIdVector.begin();
        else
            return -1;
    }

    ///@}
    double mLhsValue;
    double mRhsValue;

    EquationIdVectorType mMasterEquationIdVector;
    std::vector<double> mMasterWeightsVector;

    double mConstant;
    LockObject mLockObject;

}; // End of ConstraintEquation class

/**
 * @class LocalIndices
 * @ingroup KratosCore
 * @brief This class stores the stores three different vectors of local internal, slave, master indices
 *          which are used in constraint builder and solver.
 *
 * @author Aditya Ghantasala
 */
class LocalIndices
{
    public:
    KRATOS_CLASS_POINTER_DEFINITION(LocalIndices);
    typedef std::size_t IndexType;
    typedef std::vector<IndexType> VectorIndexType;

    /**
     * Structure "LocalIndices" to be used as an encapsulation for set of local elemental/conditional indices
    */
    LocalIndices()
    {
        internal_index_vector.resize(0);
        master_index_vector.resize(0);
        slave_index_vector.resize(0);
    }

    ~LocalIndices()
    {
        internal_index_vector.resize(0);
        master_index_vector.resize(0);
        slave_index_vector.resize(0);

        internal_index_vector.shrink_to_fit();
        master_index_vector.shrink_to_fit();
        slave_index_vector.shrink_to_fit();
    }

    /*
    * Access functions for the index vectors
    */
    VectorIndexType& InternalIndices(){return internal_index_vector;}
    VectorIndexType& MasterIndices(){return master_index_vector;}
    VectorIndexType& SlaveIndices(){return slave_index_vector;}

    private:
    VectorIndexType internal_index_vector; // indicies corresponding to internal DOFs
    VectorIndexType master_index_vector; // indicies corresponding to master DOFs
    VectorIndexType slave_index_vector; // indicies corresponding to slave DOFs
};

} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
