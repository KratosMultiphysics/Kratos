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

#if !defined(CONSTRAINT_EQUATION_H)
#define CONSTRAINT_EQUATION_H
// System includes
#include <vector>
#include <unordered_set>
#include <iostream>
#include <assert.h>

// project includes
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/node.h"
#include <boost/functional/hash.hpp>

namespace Kratos
{


/**
 * @class MasterSlaveRelation
 * @ingroup KratosCore
 * @brief This class stores the information regarding the MasterSlaveRelation equation. Naming convenction is defined like this. (each object of this class will store one equation in the given form
 * @details *   SlaveEquationId = w_1*MasterEquationId_1 + w_2*MasterEquationId_2 + ..... + w_n*MasterEquationId_n
 *
 *   each slaveDOF and each of its masterDOF have the following attributes
 *   a. dof ID
 *   b. dof KEY
 *   c. equation ID
 *
 *   one should be able to access the constraint equation both with (dofID , dofKEY) pair and/or equationID of the slave.
 *   This equation is imposed on the linear system of equations either element wise or on the global system.
 * @author Aditya Ghantasala
 */


/*
*

*
*/
class MasterSlaveRelation : public IndexedObject
{
    typedef std::size_t IndexType;
  public:
    KRATOS_CLASS_POINTER_DEFINITION(MasterSlaveRelation);
    typedef MasterData::Pointer MasterDataPointerType;

  public:

    // empty constructor and methods to add master and slave independently.
    MasterSlaveRelation() : IndexedObject(0), mSlaveDofId(0)
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }


    MasterSlaveRelation(IndexType const &rSlaveEquationId) : IndexedObject(rSlaveEquationId), mSlaveEquationId(rSlaveEquationId)
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }

    void SetConstant(double Constant) { mConstant = Constant; }
    void SetConstantUpdate(double ConstantUpdate) { mConstantUpdate = ConstantUpdate; }
    double Constant() const { return mConstant; }
    double ConstantUpdate() const { return mConstantUpdate; }
    IndexType SlaveEquationId() const { return mSlaveEquationId; }

    // This is only for serializer. Not to be used outside
    void AddMaster(std::size_t MasterEquationId, double Weight)
    {
        auto res = mMasterDataSet.find(MasterEquationId);
        if (res != mMasterDataSet.end())
        {
            (*res)->MasterWeight() += Weight;
        }
        else
        {
            mMasterDataSet[MasterEquationId] = Weight;
        }
    }

    // Get number of masters for this slave
    std::size_t GetNumberOfMasters() const
    {
        return mMasterDataSet.size();
    }

    void PrintInfo(std::ostream& Output) const override
    {
        Output << "##############################" << std::endl;
        Output << "SlaveDofID :: " << SlaveDofId() << std::endl;
        Output << "SlaveDofKey :: " << SlaveDofKey() << std::endl;
        Output << "SlaveEquationId :: " << SlaveEquationId() << std::endl;
        Output << "Constant :: " << Constant() << std::endl;
        int index = 0;
        Output << "############################## :: Masters" << std::endl;
        for (auto &master : mMasterDataSet)
        {
            Output << index << " Master  ID :: " << (*master).MasterDofId() << ", weight :: " << (*master).MasterWeight() << std::endl;
            index++;
        }
        Output << "##############################" << std::endl;
    }

    void Clear()
    {
        mMasterDataSet.clear();
    }


    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer &rSerializer) const override
    {

    }

    virtual void load(Serializer &rSerializer) override
    {

    }
    ///@}

  private:
    std::unordered_set<IndexType, double> mMasterDataSet;

    IndexType mSlaveEquationId;
    double mConstant;
    double mConstantUpdate;

}; // End of ConstraintEquation class


} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED