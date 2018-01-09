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
    typedef Dof<double> DofType;
    struct MasterData
    {
        MasterData(const DofType &iMasterDof, double iWeight=0.0) : masterDof(iMasterDof)
        {
            masterWeight = iWeight;
        }
        bool operator==(const MasterData &iOtherMasterData)
        {
            return (this->masterDof == iOtherMasterData.masterDof);
        }
        unsigned int MasterKey() const {return masterDof.GetVariable().Key();}
        double MasterWeight() const {return masterWeight;}
        unsigned int MasterDofId() const {return masterDof.Id();}
        unsigned int MasterEqId() const {return masterDof.EquationId();}
        double masterWeight;
        const DofType &masterDof;
    };

    std::vector<MasterData> mMasterDataVector;

    const DofType &slaveDof;
    mutable double constant;
    mutable double constantUpdate;

  public:
    KRATOS_CLASS_POINTER_DEFINITION(ConstraintEquation);
    typedef std::vector<MasterData>::const_iterator const_iterator;
    typedef std::vector<MasterData>::iterator iterator;
    friend class ConstraintEquationContainer;

    ConstraintEquation(const DofType &iSlaveDof) : slaveDof(iSlaveDof)
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }

    void SetConstant(double iConstant) { constant = iConstant; }
    void SetConstantUpdate(double iConstantUpdate) { constantUpdate = iConstantUpdate; }
    double Constant() const { return constant; }
    double ConstantUpdate() const { return constantUpdate; }
    unsigned int SlaveDofId() const { return slaveDof.Id(); }
    unsigned int SlaveDofKey() const { return slaveDof.GetVariable().Key(); }
    unsigned int SlaveEquationId() const { return slaveDof.EquationId(); }

    // Add a master to this slave given are the masterDofId, masterDofKey, weight
    void AddMasterData(DofType &iMasterDof, double iWeight)
    {
        MasterData iMasterData(iMasterDof, iWeight); 
        auto it = std::find(mMasterDataVector.begin(), mMasterDataVector.end(), iMasterData);
        if (it != mMasterDataVector.end()) // This master is already present so we add up the weight
        {
            it->masterWeight += iWeight;
        }
        else
        {  // No entry for this master dof ID. So we add the data
            mMasterDataVector.push_back(iMasterData);
        }
    }

    // Get number of masters for this slave
    int NumberOfMasters()
    {
        return mMasterDataVector.size();
    }
    bool operator==(const ConstraintEquation &iOtherConstraintEquation)
    {
        return (this->slaveDof == iOtherConstraintEquation.slaveDof);
    }

    void PrintInfo()
    {
        std::cout << "SlaveDofID :: " << SlaveDofId() << std::endl;
        std::cout << "SlaveDofKey :: " << SlaveDofKey() << std::endl;
        std::cout << "SlaveEquationId :: " << SlaveEquationId() << std::endl;
        std::cout << "Constant :: " << constant << std::endl;
        int index = 0;
        std::cout << "##############################" << std::endl;
        for (auto &master : mMasterDataVector)
        {
            std::cout << index << " Master  ID :: " << master.masterDof << ", weight :: " << master.masterWeight << std::endl;
            index++;
        }
        std::cout << "##############################" << std::endl;
    }

    // To make this class object iterate over all the master data vector.
    const_iterator begin() const { return mMasterDataVector.begin(); }
    const_iterator end() const { return mMasterDataVector.end(); }

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

    typedef std::vector<ConstraintEquation> ConstraintEquationVectorType;
    typedef ConstraintEquationVectorType::const_iterator const_iterator;
    typedef ConstraintEquationVectorType::iterator iterator;

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
    iterator begin() { return mDataContainer.begin(); }
    iterator end() { return mDataContainer.end(); }    

    /**
		Get the Data for this slave
		@return Data vector for this slave
		*/
    const ConstraintEquation &GetConstraintEquation(const DofType &iSlaveDof)
    {
        auto pos = std::find(mDataContainer.begin(), mDataContainer.end(), ConstraintEquation(iSlaveDof));
        return *(pos);
    }

    const ConstraintEquation &GetConstraintEquation(const unsigned int iSlaveEqutionId)
    {
        auto pos = std::find_if(mDataContainer.begin(), mDataContainer.end(), [&iSlaveEqutionId](const ConstraintEquation &obj) { return obj.SlaveEquationId() == iSlaveEqutionId; });
        return *(pos);
    }

    /**
		Get the Total number of MasterDOFs for a given slave dof
		@return Total number of MasterDOFs for a given slave dof
		 */
    unsigned int GetNumbeOfMasterDofsForSlave(const DofType &iSlaveDof)
    {
        auto pos = std::find(mDataContainer.begin(), mDataContainer.end(), iSlaveDof);
        int numMasters = -1;
        if (pos != mDataContainer.end())
            numMasters = (pos)->NumberOfMasters();
        return numMasters;
    }
    unsigned int GetNumbeOfMasterDofsForSlave(unsigned int iSlaveEqutionId)
    {
        // using find_if with a lambda to find the salve with its equation Id
        auto pos = std::find_if(mDataContainer.begin(), mDataContainer.end(), [&iSlaveEqutionId](const ConstraintEquation &obj) { return obj.SlaveEquationId() == iSlaveEqutionId; });
        int numMasters = -1;
        if (pos != mDataContainer.end())
            numMasters = (pos)->NumberOfMasters();
        return numMasters;
    }

    /**
		Get the Data for this slave
		@return Data vector for this slave
		*/
    const ConstraintEquationVectorType &GetData()
    {
        return mDataContainer;
    }

    /**
		Adds a constraints between the given slave and master with a weight. 		
		*/

    // Takes in a slave dof and a master dof
    void AddConstraint(DofType &iSlaveDof, DofType &iMasterDof, double iWeight, double constant = 0.0)
    {
        auto pos = std::find(mDataContainer.begin(), mDataContainer.end(), iSlaveDof);          
        if(pos != mDataContainer.end()){ // Equation already exists 
            pos->AddMasterData(iMasterDof, iWeight);
        } else { // Equation does not exist
            ConstraintEquation newEq = ConstraintEquation(iSlaveDof);
            newEq.AddMasterData(iMasterDof, iWeight);
            mDataContainer.push_back(newEq);
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
    ConstraintEquationVectorType mDataContainer;
    std::string mName;
    ///@}
};

///@name Input/Output funcitons
///@{

///@}

} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
