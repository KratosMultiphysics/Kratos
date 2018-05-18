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

    // empty constructor and methods to add master and slave independently.
    MasterSlaveRelation() : mId(0), mKey(0)
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }


    MasterSlaveRelation(DofType const &rSlaveDof) : mId(rSlaveDof.Id()), mKey(rSlaveDof.GetVariable().Key())
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }

    // This is only for serializer. This is not meant to be used anywhere else
    MasterSlaveRelation(std::size_t Id, std::size_t Key) : mId(Id), mKey(Key)
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }

    void SetConstant(double Constant) { mConstant = Constant; }
    void SetConstantUpdate(double ConstantUpdate) { mConstantUpdate = ConstantUpdate; }
    double Constant() const { return mConstant; }
    double ConstantUpdate() const { return mConstantUpdate; }
    std::size_t SlaveDofId() const { return mId; }
    std::size_t SlaveDofKey() const { return mKey; }
    std::size_t SlaveEquationId() const { return mEquationId; }
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
    std::size_t GetNumberOfMasters() const
    {
        return mMasterDataSet.size();
    }

    void PrintInfo(std::ostream& Output) const
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
 
    virtual void save(Serializer &rSerializer) const
    {
        rSerializer.save("slave_id", mId);            // saving the vector of the slave id
        rSerializer.save("slave_key", mKey);          // saving the vector of the slave key
        rSerializer.save("constant", mConstant);              // saving the id of the master
        rSerializer.save("constant_update", mConstantUpdate); // saving the id of the master
        rSerializer.save("num_masters", GetNumberOfMasters());    // Writint number of masters for this slave
        for (const auto &master_data : mMasterDataSet)
        {
            rSerializer.save("master_id", master_data->MasterDofId());   // saving the id of the master
            rSerializer.save("master_key", master_data->MasterDofKey()); // saving the id of the master
            rSerializer.save("weight", master_data->MasterWeight());     // saving the id of the master
        }
    }

    virtual void load(Serializer &rSerializer)
    {
        std::size_t slave_id(0), slave_key(0), num_masters(0);
        double constant(0.0), constant_update(0.0);
        rSerializer.load("slave_id", slave_id);
        rSerializer.load("slave_key", slave_key);
        rSerializer.load("constant", constant);
        rSerializer.load("constant_update", constant_update);
        rSerializer.load("num_masters", num_masters);

        for (std::size_t j = 0; j < num_masters; j++)
        {
            std::size_t master_id(0), master_key(0);
            double weight(0);
            rSerializer.load("master_id", master_id);
            rSerializer.load("master_key", master_key);
            rSerializer.load("weight", weight);
            this->AddMaster(master_id, master_key, weight);
        }
    }
    ///@}

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