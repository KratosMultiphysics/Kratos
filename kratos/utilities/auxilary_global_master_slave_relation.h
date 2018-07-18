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
#include <unordered_set>
#include <iostream>
#include <assert.h>

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
    typedef std::size_t IndexType;
    typedef Matrix MatrixType;
    typedef Vector VectorType;
    typedef std::vector<IndexType> EquationIdVectorType;
    KRATOS_CLASS_POINTER_DEFINITION(AuxilaryGlobalMasterSlaveRelation);

    // empty constructor and methods to add master and slave independently.
    AuxilaryGlobalMasterSlaveRelation() : IndexedObject(0)
    {
        mConstant = 0.0;
    }

    AuxilaryGlobalMasterSlaveRelation(IndexType const &rSlaveEquationId) : IndexedObject(rSlaveEquationId)
    {
        mConstant = 0.0;
    }

    /**
     * Function to get the slave equation Id corresponding to this constraint.
     * @param Constant the value of the constant to be assigned.
     */
    IndexType SlaveEquationId() const { return this->Id(); }

    /**
     * Function to set the lefthand side of the constraint (the slave dof value)
     * @param LhsValue the value of the lhs (the slave dof value)
     */
    void SetLHSValue(double const &LhsValue)
    {
        mLhsValue = LhsValue;
    }

    /**
     * Function to update the righthand side of the constraint (the combination of all the master dof values and constants)
     * @param RhsValue the value of the lhs (the slave dof value)
     */
    void SetRhsValue(double const &RhsValue) { mRhsValue = RhsValue; }
    void UpdateRhsValue(double const &RhsValueUpdate)
    {
        mRhsValue += RhsValueUpdate;
    }

    // Get number of masters for this slave
    IndexType GetNumberOfMasters() const
    {
        return mMasterEquationIdVector.size();
    }

    /**
     * this determines the master equation IDs connected to this constraint
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

        for (IndexType i = 0; i < this->GetNumberOfMasters(); i++)
        {
            rMasterEquationIds[i] = mMasterEquationIdVector[i];
        }
    }

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
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

        for (IndexType i = 0; i < this->GetNumberOfMasters(); i++)
            rMasterWeightsVector(i) = mMasterWeightsVector[i];


        mConstant = mRhsValue - mLhsValue;
        rConstant = mConstant;
    }

    void PrintInfo() const
    {
        std::cout<<std::endl;
        std::cout << "------------------------------" << std::endl;
        std::cout << "SlaveEquationId :: " << SlaveEquationId() << std::endl;
        std::cout << "Constant :: " << mConstant << std::endl;
        std::cout << "------------------------------ :: Masters = " << mMasterEquationIdVector.size()<< std::endl;
        for (IndexType i = 0; i< mMasterEquationIdVector.size(); i++)
        {
            std::cout << i << " Master  equation id :: " <<mMasterEquationIdVector[i] << ", weight :: " << mMasterWeightsVector[i] << std::endl;
        }
        std::cout << "------------------------------" << std::endl;
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

    virtual void save(Serializer &rSerializer) const override
    {
    }

    virtual void load(Serializer &rSerializer) override
    {
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
