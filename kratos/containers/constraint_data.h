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
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/node.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"

namespace Kratos
{
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

  private:
    struct key_hash_tuple : public std::unary_function<key_tupple, std::size_t>
    {
        std::size_t operator()(const key_tupple &k) const
        {

            std::size_t seed = 0;
            boost::hash_combine(seed, std::get<0>(k));
            boost::hash_combine(seed, std::get<1>(k));
            boost::hash_combine(seed, std::get<2>(k));
            return seed;
        }
    };

    struct key_equal_tuple : public std::binary_function<key_tupple, key_tupple, bool>
    {
        bool operator()(const key_tupple &v0, const key_tupple &v1) const
        {

            std::size_t seed0 = 0;
            boost::hash_combine(seed0, std::get<0>(v0));
            boost::hash_combine(seed0, std::get<1>(v0));
            boost::hash_combine(seed0, std::get<2>(v0));

            std::size_t seed1 = 0;
            boost::hash_combine(seed1, std::get<0>(v1));
            boost::hash_combine(seed1, std::get<1>(v1));
            boost::hash_combine(seed1, std::get<2>(v1));

            return (seed0 == seed1);
        }
    };

    struct pair_hash
    {
        template <class T1, class T2>
        std::size_t operator()(const std::pair<T1, T2> &p) const
        {
            std::size_t seed0 = 0;
            boost::hash_combine(seed0, p.first);
            boost::hash_combine(seed0, p.second);
            return seed0;
        }
    };

  public:
    typedef std::unordered_map<const key_tupple, double, key_hash_tuple, key_equal_tuple> MasterDofWeightMapType;

    ///@name Life Cycle
    ///@{

    /**
		Creates a MPC data object
		*/
    ConstraintData() : mDofConstraints(), mEquationIdToWeightsMap()
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
        mSlaveEquationIdConstantsMap.clear();
        mEquationIdToWeightsMap.clear();
    }

    /**
		Get the MasterDOFs vector for this slave
		@return MasterDOFs vector for this slave
		*/
    const MasterIdWeightMapType &GetMasterDataForSlave(DofType &SlaveDof)
    {
        return mEquationIdToWeightsMap[SlaveDof.EquationId()];
    }

    /**
		Adds a constraints between the given slave and master with a weight. 		
		*/

    // Takes in a slave dof equationId and a master dof equationId
    void AddConstraint(unsigned int SlaveDofEquationId, unsigned int MasterDofEquationId, double weight, double constant = 0.0)
    {
        mEquationIdToWeightsMap[SlaveDofEquationId].insert(std::pair<unsigned int, double>(MasterDofEquationId, weight));
        mSlaveEquationIdConstantsMap.insert(std::pair<unsigned int, double>(SlaveDofEquationId, constant));
        mSlaveEquationIdConstantsUpdate.insert(std::pair<unsigned int, double>(SlaveDofEquationId, constant));
    }

    // Takes in a slave dof and a master dof
    void AddConstraint(DofType &SlaveDof, DofType &MasterDof, double weight, double constant = 0.0)
    {
        //here we can get the dof since we are sure that such dof exist
        //auto &slave_dof = mp_model_part.Nodes(SlaveNodeId).GetDof(SlaveVariable);
        IndexType MasterNodeId = MasterDof.Id();
        unsigned int MasterVariableKey = (MasterDof).GetVariable().Key();

        unsigned int slaveVariableKey = SlaveDof.GetVariable().Key();

        mDofConstraints[std::make_pair(SlaveDof.Id(), slaveVariableKey)][std::tie(MasterNodeId, MasterVariableKey, constant)] += weight;
    }

    /**
		Get the Total number of MasterDOFs for a given slave dof
		@return Total number of MasterDOFs for a given slave dof
		 */
    unsigned int
    GetNumbeOfMasterDofsForSlave(const DofType &SlaveDof)
    {
        return mDofConstraints[std::make_pair(SlaveDof.Id(), SlaveDof.GetVariable().Key())].size();
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
    }

    ///@name Member Variables
    ///@{
    //this holds the definition of the constraint - can be constructed prior to EquationIds
    std::unordered_map<SlavePairType, MasterDofWeightMapType, pair_hash> mDofConstraints;
    //this stores a much simpler "map of maps" of EquationIds vs EquationId & weight
    // This is to be formulated inside the builder and solver before build() function ideally in initialize solution step
    std::unordered_map<unsigned int,
                       std::unordered_map<unsigned int, double>>
        mEquationIdToWeightsMap;

    std::unordered_map<unsigned int, double> mSlaveEquationIdConstantsMap;
    std::unordered_map<unsigned int, double> mSlaveEquationIdConstantsUpdate;

    bool mActive;
    std::string mName;
    ///@}

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer &rSerializer) const
    {
        rSerializer.save("MpcDataName", mName);
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
        }
    }

    virtual void load(Serializer &rSerializer)
    {
        rSerializer.load("MpcDataName", mName);
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
        }
    }

    ///@}
};

///@name Input/Output funcitons
///@{


///@}

} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
