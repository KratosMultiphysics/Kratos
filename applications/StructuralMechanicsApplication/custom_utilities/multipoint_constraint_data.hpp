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

#if !defined(MULTIPOINT_CONSTRAINT_DATA_H)
#define MULTIPOINT_CONSTRAINT_DATA_H
// System includes
#include <vector>
#include <unordered_map>
#include <iostream>
#include <tuple>
#include <utility>
#include <assert.h>

// project includes
#include <boost/functional/hash.hpp>

namespace Kratos
{
/** \brief Quaternion
	* A simple class that implements the main features of quaternion algebra
	*/
class MpcData
{

  public:
    /// Pointer definition of DataValueContainer
    KRATOS_CLASS_POINTER_DEFINITION(MpcData);

    typedef Dof<double> DofType;
    typedef VariableData VariableDataType;
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>> VariableComponentType;
    typedef unsigned int IndexType;
    typedef std::vector<Dof<double>::Pointer> DofsVectorType;
    typedef std::unordered_map<unsigned int, double> MasterIdWeightMapType;
    typedef std::pair<unsigned int, unsigned int> SlavePairType;
    typedef std::tuple<unsigned int, unsigned int, double> key_tupple;
    typedef Kratos::Variable<double> VariableType;

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

    //friend bool operator == (MpcData &obj1, MpcData &obj2);

    typedef std::unordered_map<const key_tupple, double, key_hash_tuple, key_equal_tuple> MasterDofWeightMapType;
    //typedef std::unordered_map<std::tuple<unsigned int, VariableComponentType, int>, double> ;

    ///@name Life Cycle
    ///@{

    /**
		Creates a MPC data object
		*/
    MpcData() : mDofConstraints(), mEquationIdToWeightsMap()
    {
    }
    /// Destructor.
    virtual ~MpcData(){};

    ///@}

  public:
    ///@name Operators
    ///@{

    ///@}

  public:
    ///@name Access
    ///@{

    /**
		Clears the maps contents
		*/
        void Clear(){
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
        mEquationIdToWeightsMap[SlaveDofEquationId].insert(  std::pair<unsigned int, double>(MasterDofEquationId, weight)  );
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

    // Takes in a slave dof and a list of all the masters associated with it and corresponding weights, partitionIds
    void AddConstraint(DofType &SlaveDof, DofsVectorType MasterDofsVector, std::vector<double> weightsVector, std::vector<double> ConstantVector = std::vector<double>())
    {
        //here we can get the dof since we are sure that such dof exist
        //auto &slave_dof = mp_model_part.Nodes(SlaveNodeId).GetDof(SlaveVariable);
        if (MasterDofsVector.size() != weightsVector.size())
            assert(false);

        unsigned int slaveNodeId = SlaveDof.Id();
        unsigned int slaveVariableKey = SlaveDof.GetVariable().Key();
        unsigned int index = 0;
        for (auto MasterDof : MasterDofsVector)
        {
            IndexType MasterNodeId = (*MasterDof).Id();
            unsigned int MasterVariableKey = (*MasterDof).GetVariable().Key(); // TODO :: Check why do we need a mastervariable ... is a master key not enough ?
            double constant=0.0;
            if (ConstantVector.size() == 0)
                constant = 0;
            else
                constant = ConstantVector[index];

            mDofConstraints[std::make_pair(slaveNodeId, slaveVariableKey)]
                           [std::make_tuple(MasterNodeId, MasterVariableKey, constant)] += weightsVector[index];
            ++index;
        }
    }

    /**
		Get the Total number of MasterDOFs for a given slave dof
		@return Total number of MasterDOFs for a given slave dof
		 */
    unsigned int GetNumbeOfMasterDofsForSlave(const DofType &SlaveDof)
    {
        return mDofConstraints[std::make_pair(SlaveDof.Id(), SlaveDof.GetVariable().Key())].size();
    }

    /**
		Set the name for the current set of constraints. 
		 */
    void SetName(const std::string name)
    {
        mName = name;
    }
    /**
		Get the name for the current set of constraints. 
		 */
    std::string GetName()
    {
        return mName;
    }    

    /**
		Set the activeness for current set of constraints. 
		 */
    void SetActive(const bool isActive)
    {
        mActive = isActive;
    }

    /**
		Returns true if the constraint set is active
		 */
    bool IsActive()
    {
        return mActive;
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
        std::cout << std::endl;
        std::cout << "===============================================================" << std::endl;
        std::cout << "Number of Slave DOFs :: "<< mDofConstraints.size()<<std::endl;
        for(auto i : mDofConstraints)
        {
            std::cout << "Number of Master DOFs :: "<< i.second.size()<<std::endl;
        }

        std::cout << "===============================================================" << std::endl;
        std::cout << std::endl;
    }

    ///@}
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << " MpcData object " << std::endl;
    }

  public:
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MpcData);
    }

    virtual void load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MpcData);
    }

    ///@}
};

///@name Input/Output funcitons
///@{

/*bool operator== (MpcData &obj1, MpcData &obj2) 
{
        return obj1.GetName() == obj2.GetName();
}    */

inline std::istream &operator>>(std::istream &rIStream, MpcData &rThis);

inline std::ostream &operator<<(std::ostream &rOStream, const MpcData &rThis)
{
    return rOStream;
}

///@}

} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
