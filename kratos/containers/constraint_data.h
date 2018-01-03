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
#include <unordered_map>
#include <iostream>
#include <tuple>
#include <utility>
#include <assert.h>

// project includes
#include <boost/functional/hash.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/range/iterator_range.hpp>
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/node.h"

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
    struct MasterData
    {

        MasterData()
        {
            masterDofId = -1;
            masterKey = -1;
            masterWeight = 0.0;
            masterEqId = -1;
        }
        MasterData(unsigned int iMasterDofId, unsigned int iMasterKey, double iWeight)
        {
            masterDofId = iMasterDofId;
            masterKey = iMasterKey;
            masterWeight = iWeight;
            masterEqId = -1;
        }
        MasterData(unsigned int iequationId, double iWeight)
        {
            masterDofId = -1;
            masterKey = -1;
            masterEqId = iequationId;
            masterWeight = iWeight;
        }
        unsigned int masterDofId;
        unsigned int masterKey;
        double masterWeight;
        unsigned int masterEqId;
    };

    std::vector<MasterData> mMasterData;

    unsigned int slaveDofId;
    std::size_t slaveDofKey;
    mutable unsigned int slaveEquationId;
    mutable double constant;
    mutable double constantUpdate;       

  public:
    KRATOS_CLASS_POINTER_DEFINITION(ConstraintEquation);
    typedef std::vector<MasterData>::const_iterator const_iterator;
    friend class ConstraintEquationContainer; 

    ConstraintEquation(unsigned int iDofId, std::size_t iDofKey) : slaveDofId(iDofId), slaveDofKey(iDofKey)
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }
    ConstraintEquation(unsigned int iDofId, std::size_t iDofKey, unsigned int iEquationId) : slaveDofId(iDofId), slaveDofKey(iDofKey), slaveEquationId(iEquationId)
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }

/*     ConstraintEquation(const ConstraintEquation &iOther)
    {
        this->slaveDofId = iOther.slaveDofId;
        this->slaveDofKey = iOther.slaveDofKey;
        this->slaveEquationId = iOther.slaveEquationId;
        this->constant = iOther.constant;
        this->constantUpdate = iOther.constantUpdate;

        this->mMasterData = std::move(iOther.mMasterData);
    } */

    void SetConstant(double iConstant) { constant = iConstant; }
    void SetConstantUpdate(double iConstantUpdate) { constantUpdate = iConstantUpdate; }
    void SetSlaveEquationId(unsigned int iEquationId) { slaveEquationId = iEquationId; }
    double Constant() const {return constant;}
    double ConstantUpdate() const {return constantUpdate;}
    unsigned int SlaveDofId() const { return slaveDofId; }
    unsigned int SlaveDofKey() const { return slaveDofKey; }
    unsigned int SlaveEquationId() const { return slaveEquationId; }

    // Add a master to this slave given are the masterDofId, masterDofKey, weight
    void AddMasterData(unsigned int masterDofId, std::size_t masterDofKey, double weight)
    {
        auto it = std::find_if(mMasterData.begin(), mMasterData.end(), [&masterDofId](const MasterData &obj) { return obj.masterDofId == masterDofId; });
        if (it == mMasterData.end()) // No entry for this master dof ID. So we add the data
        {
            mMasterData.push_back(MasterData(masterDofId, masterDofKey, weight));
        }
        else
        { // The dof is perviously added(probably multiple times), but not sure if with the same key.
            auto it = mMasterData.begin();
            auto condLambda = [&masterDofId](const MasterData &obj) { return obj.masterDofId == masterDofId; };
            while ((it != mMasterData.end())) // We check the key for all the occurances of this dofId
            {
                if (it->masterKey == masterDofKey) // Same pair of dofId and dofKey already exists. So add the weights.
                {
                    it->masterWeight += weight;
                    break;
                }
                it = std::find_if(std::next(it), mMasterData.end(), condLambda);
            }

            if (it == mMasterData.end()) // Here the dofId is added but not with the given key. So we add the data to the vectors.
            {
                mMasterData.push_back(MasterData(masterDofId, masterDofKey, weight));
            }
        }
    }

    // For this slave, add a masterEquationId to a master given by masterDofID and masterDofKey
    void SetMasterEquationId(unsigned int masterDofId, std::size_t masterDofKey, unsigned int masterEquationId)
    {
        auto it = mMasterData.begin();
        auto condLambda = [&masterDofId](const MasterData &obj) { return obj.masterDofId == masterDofId; };
        while ((it != mMasterData.end())) // We check the key for all the occurances of this dofId
        {
            if (it->masterKey == masterDofKey) // Same pair of dofId and dofKey already exists. So add the weights.
            {
                it->masterEqId = masterEquationId;
                break;
            }
            it = std::find_if(std::next(it), mMasterData.end(), condLambda);
        }
    }

    // Get number of masters for this slave
    int NumberOfMasters()
    {
        return mMasterData.size();
    }

    void PrintInfo()
    {
        std::cout << "SlaveDofID :: " << slaveDofId << std::endl;
        std::cout << "SlaveDofKey :: " << slaveDofKey << std::endl;
        std::cout << "SlaveEquationId :: " << slaveEquationId << std::endl;
        std::cout << "Constant :: " << constant << std::endl;
        int index = 0;
        std::cout << "##############################" << std::endl;
        for (auto &master : mMasterData)
        {
            std::cout << index << " Master  ID :: " << master.masterDofId << ", equationID :: " << master.masterEqId << ", weight :: " << master.masterWeight << std::endl;
            index++;
        }
        std::cout << "##############################" << std::endl;
    }

    // To make this class object iterate over all the master data vector.
    const_iterator begin() const { return mMasterData.begin(); }
    const_iterator end() const { return mMasterData.end(); }

}; // End of ConstraintEquation class

/** \brief ConstraintEquationContainer
 * 
 *  Each object of this class stores a set of constraint equations
 * 
 *  These constraint equations are of the form :
 *  
 *  slaveDOF = w_1*masterDOF_1 + w_2*masterDOF_2 + ..... + w_n*masterDOF_n
 * 
 *  It contains a boost multi index which will allow to access the constraint equations stored either with (dofID , dofKEY) pair or equationID of the slaveDOF.
 *  This class also provides interface to access the data of the constraint equations.
 *  
 */
class ConstraintEquationContainer
{

  public:
    /// Pointer definition of ConstraintEquationContainer
    KRATOS_CLASS_POINTER_DEFINITION(ConstraintEquationContainer);

    typedef Dof<double> DofType;
    typedef unsigned int IndexType;
    typedef Node<3> NodeType;
    typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;

    typedef ConstraintEquation::Pointer ConstraintEquationPointerType;
    struct SlaveDofId_Key
    {
    };
    struct SlaveEquationId
    {
    };
  
    typedef boost::multi_index_container<
        ConstraintEquationPointerType, boost::multi_index::indexed_by<
                                           boost::multi_index::hashed_unique<
                                               boost::multi_index::tag<SlaveDofId_Key>, boost::multi_index::composite_key<
                                                                                            ConstraintEquation, boost::multi_index::member<ConstraintEquation, unsigned int, &ConstraintEquation::slaveDofId>, boost::multi_index::member<ConstraintEquation, std::size_t, &ConstraintEquation::slaveDofKey>>>,
                                           boost::multi_index::hashed_non_unique<
                                               boost::multi_index::tag<SlaveEquationId>, boost::multi_index::member<ConstraintEquation, unsigned int, &ConstraintEquation::slaveEquationId>>>>
        ConstraintEquationMultiMapType;

    typedef ConstraintEquationMultiMapType::const_iterator const_iterator;

    ///@name Life Cycle
    ///@{

    /**
		Creates a MPC data object
		*/
    ConstraintEquationContainer() : mDataContainer()
    {
    }
    /// Destructor.
    virtual ~ConstraintEquationContainer(){};

    ///@}

    ///@name Access
    ///@{

    /**
		Clears the maps contents
		*/
    void Clear()
    {
        mDataContainer.clear();
    }

    const_iterator begin() const { return mDataContainer.begin(); }
    const_iterator end() const { return mDataContainer.end(); }

    /**
		Get the Data for this slave
		@return Data vector for this slave
		*/
    const ConstraintEquation &GetConstraintEquation(const DofType &SlaveDof)
    {
        auto &index = mDataContainer.get<SlaveDofId_Key>();
        auto pos = index.find(boost::make_tuple(SlaveDof.Id(), SlaveDof.GetVariable().Key()));
        return *(pos->get());
    }

    const ConstraintEquation &GetConstraintEquation(const unsigned int slaveEqutionId)
    {
        auto &index = mDataContainer.get<SlaveEquationId>();
        auto pos = index.find(slaveEqutionId);
        return *(pos->get());
    }

    void AddEquationIdToSlave(const DofType &SlaveDof, unsigned int iEquationId)
    {
        auto &index = mDataContainer.get<SlaveDofId_Key>();
        auto pos = index.find(boost::make_tuple(SlaveDof.Id(), SlaveDof.GetVariable().Key()));
        ConstraintEquationPointerType dummy11 = ConstraintEquationPointerType(new ConstraintEquation(*(*pos)));
        dummy11->slaveEquationId = iEquationId;
        index.replace(pos, dummy11);
    }

    /**
		Get the Total number of MasterDOFs for a given slave dof
		@return Total number of MasterDOFs for a given slave dof
		 */
    unsigned int GetNumbeOfMasterDofsForSlave(const DofType &SlaveDof)
    {
        auto &index = mDataContainer.get<SlaveDofId_Key>();
        auto pos = index.find(boost::make_tuple(SlaveDof.Id(), SlaveDof.GetVariable().Key()));
        int numMasters = -1;
        if (pos != index.end())
            numMasters = (pos->get())->NumberOfMasters();
        return numMasters;
    }
    unsigned int GetNumbeOfMasterDofsForSlave(unsigned int slaveEqutionId)
    {
        auto &index = mDataContainer.get<SlaveEquationId>();
        auto pos = index.find(slaveEqutionId);
        int numMasters = -1;
        if (pos != index.end())
            numMasters = (*pos)->NumberOfMasters();
        return numMasters;
    }

    /**
		Get the Data for this slave
		@return Data vector for this slave
		*/
    const ConstraintEquationMultiMapType &GetData()
    {
        return mDataContainer;
    }

    /**
		Adds a constraints between the given slave and master with a weight. 		
		*/

    // Takes in a slave dof and a master dof
    void AddConstraint(DofType &SlaveDof, DofType &MasterDof, double weight, double constant = 0.0)
    {
        //here we can get the dof since we are sure that such dof exist
        //auto &slave_dof = mp_model_part.Nodes(SlaveNodeId).GetDof(SlaveVariable);
        IndexType MasterNodeId = MasterDof.Id();
        unsigned int MasterVariableKey = (MasterDof).GetVariable().Key();
        unsigned int slaveVariableKey = SlaveDof.GetVariable().Key();

        ConstraintEquationPointerType dummy = ConstraintEquationPointerType(new ConstraintEquation(SlaveDof.Id(), slaveVariableKey));
        dummy->SetConstant(constant);
        dummy->AddMasterData(MasterNodeId, MasterVariableKey, weight);
        dummy->SetSlaveEquationId(0);

        std::pair<ConstraintEquationMultiMapType::iterator, bool> ret = mDataContainer.insert(dummy);
        ConstraintEquationMultiMapType::iterator pos = ret.first;
        if (!ret.second)
        {
            (pos->get())->AddMasterData(MasterNodeId, MasterVariableKey, weight);
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
    virtual void GetInfo() const
    {
    }

    ///@}
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << " ConstraintEquationContainer object " << std::endl;
        rOStream << " Number of constraint equations : " << mDataContainer.size() << std::endl;
    }

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer &rSerializer) const
    {
    }

    virtual void load(Serializer &rSerializer)
    {
    }

  private:
    ///@name Member Variables
    ///@{
    ConstraintEquationMultiMapType mDataContainer;
    std::string mName;
    ///@}
};

///@name Input/Output funcitons
///@{

///@}

} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
