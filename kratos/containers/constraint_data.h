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
//
//

#if !defined(CONSTRAINT_DATA_H)
#define CONSTRAINT_DATA_H
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

// ConstraintEquation class start
/*
*   This class stores the information regarding the constraint equation. Naming convenction is defined like this.
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
class ConstraintEquation
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
    KRATOS_CLASS_POINTER_DEFINITION(ConstraintEquation);
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

    // To make this class object iterate over all the master data vector.
    iterator begin() { return mMasterDataSet.begin(); }
    iterator end() { return mMasterDataSet.end(); }

  private:
    std::unordered_set<MasterDataPointerType, MasterHasher, MasterComparator> mMasterDataSet;

    const std::size_t mId;
    const std::size_t mKey;
    std::size_t mEquationId;
    double mConstant;
    double mConstantUpdate;

}; // End of ConstraintEquation class

/** \brief ConstraintEquationContainer
 * 
 *  Each object of this class stores a set of constraint equations
 * 
 *  These constraint equations are of the form :
 *  
 *  slaveDOF = w_1*masterDOF_1 + w_2*masterDOF_2 + ..... + w_n*masterDOF_n + c
 * 
 *  It contains a data structure which will allow to access the constraint equations stored either with slaveDOF itself or equationID of the slaveDOF.
 *  This class also provides interface to access the data of the constraint equations.
 *  
 */
class ConstraintEquationContainer
{

  public:
    /// Pointer definition of ConstraintEquationContainer
    KRATOS_CLASS_POINTER_DEFINITION(ConstraintEquationContainer);

    typedef Dof<double> DofType;
    typedef std::size_t IndexType;
    typedef Node<3> NodeType;
    typedef ConstraintEquation::Pointer ConstraintEquationPointerType;

    // Custom hash function to store Constraint equation objects in a unordered_set
    struct EquationHasher
    {
        std::size_t
        operator()(const ConstraintEquationPointerType &rObj) const
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, rObj->SlaveDofId());
            boost::hash_combine(seed, rObj->SlaveDofKey());
            return seed;
        }
    };

    // Custom comparator that compares the Constraint equation objects by their SlaveDOF key and SlaveDOF ID
    struct EquationComparator
    {
        bool
        operator()(const ConstraintEquationPointerType &rObj1, const ConstraintEquationPointerType &rObj2) const
        {
            return (rObj1->SlaveDofId() == rObj2->SlaveDofId()) && (rObj1->SlaveDofKey() == rObj2->SlaveDofKey());
        }
    };

    typedef std::unordered_set<ConstraintEquationPointerType, EquationHasher, EquationComparator> ConstraintEquationSetType;

    typedef ConstraintEquationSetType::iterator iterator;

    ///@name Life Cycle
    ///@{

    ///@}

    ///@name Access
    ///@{

    /**
		Clears the maps contents
		*/
    void Clear()
    {
        mConstraintEquationsSet.clear();
    }
    iterator begin() { return mConstraintEquationsSet.begin(); }
    iterator end() { return mConstraintEquationsSet.end(); }

    /**
		Get the Constraint equation for this slave
		@return Data vector for this slave
		*/
    ConstraintEquation &GetConstraintEquation(DofType &rSlaveDof)
    {
        //ConstraintEquationPointerType equation = Kratos::make_shared<ConstraintEquation>(rSlaveDof);
        auto res = mConstraintEquationsSet.find(Kratos::make_shared<ConstraintEquation>(rSlaveDof));
        if (res != mConstraintEquationsSet.end())
        {
            return *(*res);
        }
        else
        {
            KRATOS_ERROR<<"No constraint equation found for given slave dof.\nThere is something worng"<<std::endl;
        }
    }

    /**
		Get the Total number of MasterDOFs for a given slave dof
		@return Total number of MasterDOFs for a given slave dof
		 */
    int GetNumberOfMasterDofsForSlave(DofType const &rSlaveDof)
    {
        //ConstraintEquationPointerType equation = Kratos::make_shared<ConstraintEquation>(rSlaveDof);
        auto res = mConstraintEquationsSet.find(Kratos::make_shared<ConstraintEquation>(rSlaveDof));
        int numMasters = -1;
        if (res != mConstraintEquationsSet.end())
        {
            numMasters = (*res)->NumberOfMasters();
        }
        return numMasters;
    }

    /**
		Get the Data for this slave
		@return Data vector for this slave
		*/
    ConstraintEquationSetType &GetData()
    {
        return mConstraintEquationsSet;
    }

    /**
		Adds a constraints between the given slave and master with a weight. 		
		*/

    // Takes in a slave dof and a master dof
    void AddConstraint(DofType const &rSlaveDof, DofType const &rMasterDof, double Weight, double Constant = 0.0)
    {
        auto res = mConstraintEquationsSet.find(Kratos::make_shared<ConstraintEquation>(rSlaveDof));
        if (res != mConstraintEquationsSet.end())
        { // Equation already exists
            ((*res))->AddMaster(rMasterDof, Weight);
        }
        else
        {
            ConstraintEquationPointerType equation = Kratos::make_shared<ConstraintEquation>(rSlaveDof);
            equation->AddMaster(rMasterDof, Weight);
            equation->SetConstant(Constant);
            mConstraintEquationsSet.insert(equation);
        }
    }
    ///@

    ///@name Static Operations
    ///
    //@{
    /**
		 * Returns the string containing a detailed description of this object.
		 * @return the string with informations
		 */
    virtual std::string GetInfo() const
    {
        return "Object of ConstraintEquationContainer";
    }

    ///@}
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << " ConstraintEquationContainer object " << std::endl;
        rOStream << " Number of constraint equations : " << mConstraintEquationsSet.size() << std::endl;
    }

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer &rSerializer) const
    {
        rSerializer.save("NumConstraintEquations", mConstraintEquationsSet.size());
        for (const auto &equation : mConstraintEquationsSet)
        {
            rSerializer.save("slave_id", equation->SlaveDofId());            // saving the vector of the slave id
            rSerializer.save("slave_key", equation->SlaveDofKey());          // saving the vector of the slave key
            rSerializer.save("constant", equation->Constant());              // saving the id of the master
            rSerializer.save("constant_update", equation->ConstantUpdate()); // saving the id of the master
            rSerializer.save("num_masters", equation->NumberOfMasters());    // Writint number of masters for this slave
            for (const auto &master_data : *equation)
            {
                rSerializer.save("master_id", master_data->MasterDofId());   // saving the id of the master
                rSerializer.save("master_key", master_data->MasterDofKey()); // saving the id of the master
                rSerializer.save("weight", master_data->MasterWeight());     // saving the id of the master
            }
        }
    }

    virtual void load(Serializer &rSerializer)
    {
        std::size_t num_constraint_equations = 0;
        rSerializer.load("NumConstraintEquations", num_constraint_equations);
        for (std::size_t i = 0; i < num_constraint_equations; i++)
        {
            std::size_t slave_id(0), slave_key(0), num_masters(0);
            double constant(0.0), constant_update(0.0);
            rSerializer.load("slave_id", slave_id);
            rSerializer.load("slave_key", slave_key);
            rSerializer.load("constant", constant);
            rSerializer.load("constant_update", constant_update);
            rSerializer.load("num_masters", num_masters);

            ConstraintEquationPointerType equation = Kratos::make_shared<ConstraintEquation>(slave_id, slave_key);

            for (std::size_t j = 0; j < num_masters; j++)
            {
                std::size_t master_id(0), master_key(0);
                double weight(0);
                rSerializer.load("master_id", master_id);
                rSerializer.load("master_key", master_key);
                rSerializer.load("weight", weight);
                equation->AddMaster(master_id, master_key, weight);
            }
            mConstraintEquationsSet.insert(equation);
        }
    }

  private:
    ///@name Member Variables
    ///@{
    ConstraintEquationSetType mConstraintEquationsSet;
    ///@}
};

///@name Input/Output funcitons
///@{

///@}

} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
