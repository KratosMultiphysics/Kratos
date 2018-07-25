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

#if !defined(AUXILARY_GLOBAL_MASTER_SLAVE_RELATION)
#define AUXILARY_GLOBAL_MASTER_SLAVE_RELATION
// System includes
#include <vector>

// project includes
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/node.h"
#include "includes/process_info.h"

namespace Kratos
{

/**
 * @class AuxilaryGlobalMasterSlaveRelation
 * @ingroup KratosCore
 * @brief This class stores the information regarding the AuxilaryGlobalMasterSlaveRelation equation.
 *         Naming convenction is defined like this. (each object of this class will store one equation in the given form
 *
 *   SlaveEquationId = w_1*MasterEquationId_1 + w_2*MasterEquationId_2 + ..... + w_n*MasterEquationId_n
 *
 *   This stores the condensed form of the MasterSlaveConstraint objects into one object. if only one relation for a slave is added as
 *   MasterSlaveConstraint then there will only be one entry for master for its corresponding AuxilaryGlobalMasterSlaveRelation.
 *   Currently this class is designed to hold only one equation. There is only one unique object of this class for each slave. 
 *
 *   Future plan is to also make it possible to work with matrices (T) and vectors (for slave and master equation ids and constants)
 *
 *
 *  IMPORTANT : This is not seen by the user. This is a helper data structure which is exists only in the builder and solver. 
 *
 * @author Aditya Ghantasala
 */
class AuxilaryGlobalMasterSlaveRelation : public IndexedObject
{
  public:
    typedef IndexedObject BaseType;
    typedef BaseType::IndexType IndexType;
    typedef Matrix MatrixType;
    typedef Vector VectorType;
    typedef std::vector<IndexType> EquationIdVectorType;
    KRATOS_CLASS_POINTER_DEFINITION(AuxilaryGlobalMasterSlaveRelation);

    /**
     * @brief Constructor of the class
     * @param SlaveEquationId the slave equation id for which this class is being constructed.
     */
    AuxilaryGlobalMasterSlaveRelation(IndexType SlaveEquationId = 0) : IndexedObject(SlaveEquationId), mConstant(0.0)
    {
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
    void SetLHSValue(double const &LhsValue)
    {
        mLhsValue = LhsValue;
    }

    /**
     * @brief Function to update the righthand side of the constraint (the combination of all the master dof values and constants)
     * @param RHSValue the value of the lhs (the slave dof value)
     */
    void SetRHSValue(double const &RhsValue) { mRhsValue = RhsValue; }
    void UpdateRHSValue(double const &RhsValueUpdate)
    {
        mRhsValue += RhsValueUpdate;
    }

    // Get number of masters for this slave
    IndexType GetNumberOfMasters() const
    {
        return mMasterEquationIdVector.size();
    }

    /**
     * @brief this determines the master equation IDs connected to this constraint
     * @param rResult the elemental equation ID vector
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void EquationIdVector(IndexType &rSlaveEquationId,
                                  EquationIdVectorType &rMasterEquationIds,
                                  ProcessInfo &rCurrentProcessInfo)
    {
        if (rMasterEquationIds.size() != 0 || rMasterEquationIds.size() == 0)
            rMasterEquationIds.resize(this->GetNumberOfMasters(), false);

        rSlaveEquationId = this->SlaveEquationId();

        for (IndexType i = 0; i < this->GetNumberOfMasters(); ++i)
        {
            rMasterEquationIds[i] = mMasterEquationIdVector[i];
        }
    }

    /**
     * @brief   this is called during the assembling process in order
     *          to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rMasterWeightsVector the elemental left hand side matrix
     * @param rConstant the elemental right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void CalculateLocalSystem(VectorType &rMasterWeightsVector,
                                      double &rConstant,
                                      ProcessInfo &rCurrentProcessInfo)
    {
        if (rMasterWeightsVector.size() != 0 || rMasterWeightsVector.size() == 0)
        {
            rMasterWeightsVector.resize(this->GetNumberOfMasters(), false);
        }

        for (IndexType i = 0; i < this->GetNumberOfMasters(); ++i)
            rMasterWeightsVector(i) = mMasterWeightsVector[i];


        mConstant = mRhsValue - mLhsValue;
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
        mMasterEquationIdVector.clear();
        mMasterWeightsVector.clear();
    }

    void AddMaster(IndexType MasterEquationId, double Weight)
    {
        int index = MasterEquationIdExists(MasterEquationId);
        if (index > 0)
        {
            mMasterWeightsVector[index] += Weight;
        } else
        {
            mMasterEquationIdVector.push_back(MasterEquationId);
            mMasterWeightsVector.push_back(Weight);
        }
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

    std::vector<IndexType> mMasterEquationIdVector;
    std::vector<double> mMasterWeightsVector;

    double mConstant;

}; // End of ConstraintEquation class

} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
