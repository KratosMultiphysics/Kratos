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
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"

namespace Kratos
{

// SlaveData class start
class SlaveData
{
  public:
    KRATOS_CLASS_POINTER_DEFINITION(SlaveData);

    SlaveData(unsigned int iDofId, std::size_t iDofKey) : dofId(iDofId), dofKey(iDofKey)
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }
    SlaveData(unsigned int iDofId, std::size_t iDofKey, unsigned int iEquationId) : dofId(iDofId), dofKey(iDofKey), equationId(iEquationId)
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }

    SlaveData(const SlaveData &iOther)
    {
        this->dofId = iOther.dofId;
        this->dofKey = iOther.dofKey;
        this->equationId = iOther.equationId;
        this->constant = iOther.constant;
        this->constantUpdate = iOther.constant;

        this->masterDofIds = std::move(iOther.masterDofIds);
        this->masterDofKeys = std::move(iOther.masterDofKeys);
        this->masterWeights = std::move(iOther.masterWeights);
        this->masterEquationIds = std::move(iOther.masterEquationIds);
    }

    void SetConstant(double iConstant) { constant = iConstant; }
    void SetConstantUpdate(double iConstantUpdate) { constantUpdate = iConstantUpdate; }
    void SetEquationId(unsigned int iEquationId) { equationId = iEquationId; }

    // Add a master to this slave given are the masterDofId, masterDofKey, weight
    void AddMasterData(unsigned int masterDofId, std::size_t masterDofKey, double weight)
    {
        auto it = std::find(masterDofIds.begin(), masterDofIds.end(), masterDofId);
        if (it == masterDofIds.end()) // No entry for this master dof ID. So we add the data
        {
            masterDofIds.push_back(masterDofId);
            masterDofKeys.push_back(masterDofKey);
            masterWeights.push_back(weight);
            masterEquationIds.push_back(0);
        }
        else
        { // The dof is perviously added(probably multiple times), but not sure if with the same key.
            it = masterDofIds.begin();
            while ((it = std::find(it, masterDofIds.end(), masterDofId)) != masterDofIds.end()) // We check the key for all the occurances of this dofId
            {
                int masterIndex = std::distance(masterDofIds.begin(), it);
                auto masterKey = masterDofKeys[masterIndex];
                if (masterKey == masterDofKey) // Same pair of dofId and dofKey already exists. So add the weights.
                {
                    masterWeights[masterIndex] += weight;
                    break;
                }
                it++;
            }
            if (it == masterDofIds.end()) // Here the dofId is added but not with the given key. So we add the data to the vectors.
            {
                masterDofIds.push_back(masterDofId);
                masterDofKeys.push_back(masterDofKey);
                masterWeights.push_back(weight);
                masterEquationIds.push_back(0);
            }
        }
    }

    // For this slave, add a masterEquationId to a master given by masterDofID and masterDofKey
    void SetMasterEquationId(unsigned int masterDofId, std::size_t masterDofKey, unsigned int masterEquationId)
    {
        auto it = masterDofIds.begin();
        while ((it = std::find(it, masterDofIds.end(), masterDofId)) != masterDofIds.end()) // We check the key for all the occurances of this dofId
        {
            int masterIndex = std::distance(masterDofIds.begin(), it);
            auto masterKey = masterDofKeys[masterIndex];
            if (masterKey == masterDofKey) // Same pair of dofId and dofKey already exists. So add the weights.
            {
                masterEquationIds[masterIndex] = masterEquationId;
                break;
            }
            it++;
        }
    }

    // Get number of masters for this slave
    int NumberOfMasters()
    {
        return masterDofIds.size();
    }

    void PrintInfo()
    {
        std::cout << "SlaveDofID :: " << dofId << std::endl;
        std::cout << "SlaveDofKey :: " << dofKey << std::endl;
        std::cout << "SlaveEquationId :: " << equationId << std::endl;
        std::cout << "Constant :: " << constant << std::endl;
        int index = 0;
        std::cout << "##############################" << std::endl;
        for (auto &master : masterDofIds)
        {
            std::cout << index << " Master  ID :: " << master << ", equationID :: " << masterEquationIds[index] << ", weight :: " << masterWeights[index] << constant << std::endl;
            index++;
        }
        std::cout << "##############################" << std::endl;
    }

    unsigned int dofId;
    std::size_t dofKey;
    unsigned int equationId;
    double constant;
    double constantUpdate;

    std::vector<unsigned int> masterDofIds;
    std::vector<std::size_t> masterDofKeys;
    std::vector<double> masterWeights;
    std::vector<unsigned int> masterEquationIds;
}; // End of SlaveData class

/** \brief MpcData
	* A class that implements the data structure needed for applying constraints.
    */
class ConstraintData
{

  public:
    /// Pointer definition of ConstraintData
    KRATOS_CLASS_POINTER_DEFINITION(ConstraintData);

    typedef Dof<double> DofType;
    typedef VariableData VariableDataType;
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>> VariableComponentType;
    typedef unsigned int IndexType;
    typedef std::vector<Dof<double>::Pointer> DofsVectorType;
    typedef std::unordered_map<unsigned int, double> MasterIdWeightMapType;
    typedef std::pair<unsigned int, unsigned int> SlavePairType;
    typedef std::tuple<unsigned int, unsigned int, double> key_tupple;
    typedef Kratos::Variable<double> VariableType;
    typedef Node<3> NodeType;
    typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;
    typedef std::pair<unsigned int, unsigned int> MasterIdKeyPairType;

    typedef SlaveData::Pointer SlaveDataPointerType;
    struct SlaveDofId_Key
    {
    };
    struct SlaveEquationId
    {
    };
    typedef boost::multi_index_container<
        SlaveDataPointerType, boost::multi_index::indexed_by<
                                  boost::multi_index::hashed_unique<
                                      boost::multi_index::tag<SlaveDofId_Key>, boost::multi_index::composite_key<
                                                                                   SlaveData, boost::multi_index::member<SlaveData, unsigned int, &SlaveData::dofId>, boost::multi_index::member<SlaveData, std::size_t, &SlaveData::dofKey>>>,
                                  boost::multi_index::hashed_non_unique<
                                      boost::multi_index::tag<SlaveEquationId>, boost::multi_index::member<SlaveData, unsigned int, &SlaveData::equationId>>>>
        SlaveDataMultiMapType;

    typedef SlaveDataMultiMapType::iterator iterator;
    typedef SlaveDataMultiMapType::const_iterator const_iterator;

    ///@name Life Cycle
    ///@{

    /**
		Creates a MPC data object
		*/
    ConstraintData() : mDataContainer()
    {
    }
    /// Destructor.
    virtual ~ConstraintData(){};

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

    iterator begin() { return mDataContainer.begin(); }
    const_iterator begin() const { return mDataContainer.begin(); }
    iterator end() { return mDataContainer.end(); }
    const_iterator end() const { return mDataContainer.end(); }

    /**
		Get the Data for this slave
		@return Data vector for this slave
		*/
    const SlaveData &GetSlaveData(DofType &SlaveDof)
    {
        auto &index = mDataContainer.get<SlaveDofId_Key>();
        auto pos = index.find(boost::make_tuple(SlaveDof.Id(), SlaveDof.GetVariable().Key()));
        return *(pos->get());
    }

    const SlaveData &GetSlaveData(const unsigned int slaveEqutionId)
    {
        auto &index = mDataContainer.get<SlaveEquationId>();
        auto pos = index.find(slaveEqutionId);
        return *(pos->get());
    }

    void AddEquationIdToSlave(const DofType &SlaveDof, unsigned int iEquationId)
    {
        auto &index = mDataContainer.get<SlaveDofId_Key>();
        auto pos = index.find(boost::make_tuple(SlaveDof.Id(), SlaveDof.GetVariable().Key()));
		SlaveDataPointerType dummy11 = SlaveDataPointerType(new SlaveData(*(*pos)));
        dummy11->equationId = iEquationId;
        index.replace(pos,dummy11); 
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
        numMasters = (pos->get())->NumberOfMasters();
        return numMasters;
    }
    unsigned int GetNumbeOfMasterDofsForSlave(unsigned int slaveEqutionId)
    {
        auto &index = mDataContainer.get<SlaveEquationId>();
        auto pos = index.find(slaveEqutionId);
        int numMasters = -1;
        if(pos != index.end())
            numMasters = (*pos)->NumberOfMasters();
        else 
            numMasters = 0;
        return numMasters;
    }

    /**
		Get the Data for this slave
		@return Data vector for this slave
		*/
    const SlaveDataMultiMapType &GetData()
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

        SlaveDataPointerType dummy = SlaveDataPointerType(new SlaveData(SlaveDof.Id(), slaveVariableKey));
        dummy->SetConstant(constant);
        dummy->AddMasterData(MasterNodeId, MasterVariableKey, weight);
        dummy->SetEquationId(0);

        std::pair<SlaveDataMultiMapType::iterator, bool> ret = mDataContainer.insert(dummy);
        SlaveDataMultiMapType::iterator pos = ret.first;
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
        rOStream << " MpcData object " << std::endl;
        rOStream << " Number of slaves : " << mDataContainer.size() << std::endl;
    }

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer &rSerializer) const
    {
        /*         rSerializer.save("MpcDataName", mName);
        rSerializer.save("NumConstraints", mDofConstraints.size());
        for (const auto &slaveMasterrelation : mDofConstraints)
        {

            rSerializer.save("slaveID", (slaveMasterrelation.first).first);   // saving the vector of the slave id
            rSerializer.save("slaveKey", (slaveMasterrelation.first).second); // saving the vector of the slave key

            rSerializer.save("numMasters", (slaveMasterrelation.second).size()); // Writint number of masters for this slave
            for (const auto &masterIdKeyConstant : (slaveMasterrelation.second))
            {
                rSerializer.save("masterID", std::get<0>(masterIdKeyConstant.first));  // saving the id of the master
                rSerializer.save("masterKey", std::get<1>(masterIdKeyConstant.first)); // saving the id of the master
                rSerializer.save("constant", std::get<2>(masterIdKeyConstant.first));  // saving the id of the master

                rSerializer.save("weight", masterIdKeyConstant.second); // saving the id of the master
            }
        } */
    }

    virtual void load(Serializer &rSerializer)
    {
        /*         rSerializer.load("MpcDataName", mName);
        int numConstraints = 0;
        rSerializer.load("NumConstraints", numConstraints);
        for (int i = 0; i < numConstraints; i++)
        {
            int slaveID(0), slaveKey(0), numMasters(0);
            rSerializer.load("slaveID", slaveID);
            rSerializer.load("slaveKey", slaveKey);
            rSerializer.load("numMasters", numMasters);
            for (int j = 0; j < numMasters; j++)
            {
                int masterID(0), masterKey(0);
                double constant(0), weight(0);

                rSerializer.load("masterID", masterID);
                rSerializer.load("masterKey", masterKey);
                rSerializer.load("constant", constant);
                rSerializer.load("weight", weight);

                mDofConstraints[std::make_pair(slaveID, slaveKey)][std::tie(masterID, masterKey, constant)] += weight;
            }
        } */
    }

  private:
    ///@name Member Variables
    ///@{
    SlaveDataMultiMapType mDataContainer;
    std::string mName;
    ///@}
};

///@name Input/Output funcitons
///@{

///@}

} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
