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

// MasterSlaveRelation class start
/*
*   This class stores the information regarding the MasterSlaveRelation equation. Naming convenction is defined like this.
*   (each object of this class will store one equation in the given form)
*
*   slaveDOF = w_1*masterDOF_1 + w_2*masterDOF_2 + ..... + w_n*masterDOF_n
*   
*   each slaveDOF and each of its masterDOF have the following attributes 
*   a. dof ID
*   b. dof KEY
*   c. equation ID
*
*   one should be able to access the constraint equation both with (dofID , dofKEY) pair and/or equationID of the slave.
*   This equation is imposed on the linear system of equations either element wise or on the global system.
*
*/
class MasterSlaveRelation
{
  private:
    typedef Dof<double> DofType;
    struct MasterData
    {
        KRATOS_CLASS_POINTER_DEFINITION(MasterData);
        MasterData(DofType const &rMasterDof, double Weight = 0.0) : mMasterWeight(Weight), mId(rMasterDof.Id()), mKey(rMasterDof.GetVariable().Key())
        {
        }
        // This is only for serializer. This is not meant to be used anywhere else
        MasterData(std::size_t Id, std::size_t Key, double Weight = 0.0) : mMasterWeight(Weight), mId(Id), mKey(Key)
        {
        }
        std::size_t MasterDofKey() { return mKey; }
        double MasterWeight() const { return mMasterWeight; }
        double &MasterWeight() { return mMasterWeight; }
        std::size_t MasterDofId() { return mId; }
        std::size_t MasterEqId() { return mEquationId; }
        void SetMasterEqId(std::size_t Id) { mEquationId = Id; }

      private:
        double mMasterWeight;
        const std::size_t mId;
        const std::size_t mKey;
        std::size_t mEquationId;
    };

  public:
    KRATOS_CLASS_POINTER_DEFINITION(MasterSlaveRelation);
    typedef MasterData::Pointer MasterDataPointerType;

  private:
    // Custom hash function to store MasterData objects in a unordered_set
    struct MasterHasher
    {
        std::size_t
        operator()(const MasterDataPointerType &rObj) const
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, rObj->MasterDofId());
            boost::hash_combine(seed, rObj->MasterDofKey());
            return seed;
        }
    };

    // Custom comparator that compares the MasterData objects by their key and ID
    struct MasterComparator
    {
        bool
        operator()(const MasterDataPointerType &rObj1, const MasterDataPointerType &rObj2) const
        {
            return (rObj1->MasterDofId() == rObj2->MasterDofId()) && (rObj1->MasterDofKey() == rObj2->MasterDofKey());
        }
    };

  public:
    typedef std::unordered_set<MasterDataPointerType, MasterHasher, MasterComparator>::iterator iterator;  
    ConstraintEquation(DofType const &rSlaveDof) : mId(rSlaveDof.Id()), mKey(rSlaveDof.GetVariable().Key())
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }

    // This is only for serializer. This is not meant to be used anywhere else
    ConstraintEquation(std::size_t Id, std::size_t Key) : mId(Id), mKey(Key)
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }

    void SetConstant(double Constant) { mConstant = Constant; }
    void SetConstantUpdate(double ConstantUpdate) { mConstantUpdate = ConstantUpdate; }
    double Constant() { return mConstant; }
    double ConstantUpdate() { return mConstantUpdate; }
    std::size_t SlaveDofId() { return mId; }
    std::size_t SlaveDofKey() { return mKey; }
    std::size_t SlaveEquationId() { return mEquationId; }
    void SetSlaveEquationId(std::size_t Id) { mEquationId = Id; }

    // Add a master or update a master(if already present) to this slave given are the masterDofId, masterDofKey, weight
    void AddMaster(DofType const &rMasterDof, double Weight)
    {
        MasterDataPointerType master_data = Kratos::make_shared<MasterData>(rMasterDof, Weight);
        auto res = mMasterDataSet.find(master_data);
        if (res != mMasterDataSet.end())
        {
            (*res)->MasterWeight() += Weight;
        }
        else
        {
            mMasterDataSet.insert(master_data);
        }
    }

    // This is only for serializer. Not to be used outside
    void AddMaster(std::size_t MasterDofId, std::size_t MasterDofKey, double Weight)
    {
        MasterDataPointerType master_data = Kratos::make_shared<MasterData>(MasterDofId, MasterDofKey, Weight);
        auto res = mMasterDataSet.find(master_data);
        if (res != mMasterDataSet.end())
        {
            (*res)->MasterWeight() += Weight;
        }
        else
        {
            mMasterDataSet.insert(master_data);
        }
    }

    // Get number of masters for this slave
    std::size_t NumberOfMasters()
    {
        return mMasterDataSet.size();
    }

    void PrintInfo()
    {
        std::cout << "##############################" << std::endl;
        std::cout << "SlaveDofID :: " << SlaveDofId() << std::endl;
        std::cout << "SlaveDofKey :: " << SlaveDofKey() << std::endl;
        std::cout << "SlaveEquationId :: " << SlaveEquationId() << std::endl;
        std::cout << "Constant :: " << Constant() << std::endl;
        int index = 0;
        std::cout << "############################## :: Masters" << std::endl;
        for (auto &master : mMasterDataSet)
        {
            std::cout << index << " Master  ID :: " << (*master).MasterDofId() << ", weight :: " << (*master).MasterWeight() << std::endl;
            index++;
        }
        std::cout << "##############################" << std::endl;
    }

  private:
    std::unordered_set<MasterDataPointerType, MasterHasher, MasterComparator> mMasterDataSet;

    const std::size_t mId;
    const std::size_t mKey;
    std::size_t mEquationId;
    double mConstant;
    double mConstantUpdate;

}; // End of ConstraintEquation class


} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED