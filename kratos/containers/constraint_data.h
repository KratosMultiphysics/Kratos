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
        KRATOS_CLASS_POINTER_DEFINITION(MasterData);
        MasterData(DofType const &iMasterDof, double iWeight = 0.0) : masterWeight(iWeight), id(iMasterDof.Id()), key(iMasterDof.GetVariable().Key())
        {
        }
        friend bool operator==(const MasterData::Pointer &iMasterDataOne, const DofType &iMasterDof)
        {
            return ((iMasterDataOne->MasterDofId() == iMasterDof.Id()) && (iMasterDataOne->MasterKey() == iMasterDof.GetVariable().Key()));
        }
        unsigned int MasterKey() { return key; }
        double MasterWeight() { return masterWeight; }
        void UpdateWeight(double iWeight) { masterWeight += iWeight; }
        unsigned int MasterDofId() { return id; }
        void SetEquationId(unsigned int iId) { equationId = iId; }
        unsigned int MasterEqId() { return equationId; }
        void SetMasterEqId(unsigned int iId) { equationId = iId; }

      private:
        double masterWeight;
        const unsigned int id;
        const unsigned int key;
        unsigned int equationId;
    };
    typedef MasterData::Pointer MasterDataPointerType;
    std::vector<MasterDataPointerType> mMasterDataVector;

    const unsigned int id;
    const unsigned int key;
    unsigned int equationId;
    double constant;
    double constantUpdate;

  public:
    KRATOS_CLASS_POINTER_DEFINITION(ConstraintEquation);
    typedef std::vector<MasterDataPointerType>::iterator iterator;
    friend class ConstraintEquationContainer;

    ConstraintEquation(DofType const &iSlaveDof) : id(iSlaveDof.Id()), key(iSlaveDof.GetVariable().Key())
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }

    void SetConstant(double iConstant) { constant = iConstant; }
    void SetConstantUpdate(double iConstantUpdate) { constantUpdate = iConstantUpdate; }
    double Constant() { return constant; }
    double ConstantUpdate() { return constantUpdate; }
    unsigned int SlaveDofId() { return id; }
    unsigned int SlaveDofKey() { return key; }
    unsigned int SlaveEquationId() { return equationId; }
    void SetSlaveEquationId(unsigned int iId) { equationId = iId; }

    // Add a master to this slave given are the masterDofId, masterDofKey, weight
    void AddMasterData(DofType const &iMasterDof, double iWeight)
    {
        MasterDataPointerType masterData = MasterDataPointerType(new MasterData(iMasterDof, iWeight));
        auto it = std::find(mMasterDataVector.begin(), mMasterDataVector.end(), iMasterDof);
        if (it != mMasterDataVector.end()) // This master is already present so we add up the weight
        {
            (*it)->UpdateWeight(iWeight);
        }
        else
        { // No entry for this master dof ID. So we add the data
            mMasterDataVector.push_back(masterData);
        }
    }

    // Get number of masters for this slave
    unsigned int NumberOfMasters()
    {
        return mMasterDataVector.size();
    }
    friend bool operator==(const ConstraintEquation::Pointer &iConstraintDataOne, const DofType &iSlaveDof)
    {
        return ((iConstraintDataOne->SlaveDofId() == iSlaveDof.Id()) && (iConstraintDataOne->SlaveDofKey() == iSlaveDof.GetVariable().Key()));
    }

    void PrintInfo()
    {
        std::cout << "SlaveDofID :: " << SlaveDofId() << std::endl;
        std::cout << "SlaveDofKey :: " << SlaveDofKey() << std::endl;
        std::cout << "SlaveEquationId :: " << SlaveEquationId() << std::endl;
        std::cout << "Constant :: " << Constant() << std::endl;
        int index = 0;
        std::cout << "##############################" << std::endl;
        for (auto &master : mMasterDataVector)
        {
            std::cout << index << " Master  ID :: " << (*master).MasterDofId() << ", weight :: " << (*master).MasterWeight() << std::endl;
            index++;
        }
        std::cout << "##############################" << std::endl;
    }

    // To make this class object iterate over all the master data vector.
    iterator begin() { return mMasterDataVector.begin(); }
    iterator end() { return mMasterDataVector.end(); }

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
    typedef ConstraintEquation::Pointer ConstraintEquationPointerType;
    typedef std::vector<ConstraintEquationPointerType> ConstraintEquationVectorType;

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
    iterator begin() { return mDataContainer.begin(); }
    iterator end() { return mDataContainer.end(); }

    /**
		Get the Data for this slave
		@return Data vector for this slave
		*/
    ConstraintEquation &GetConstraintEquation(DofType &iSlaveDof)
    {
        auto pos = std::find(mDataContainer.begin(), mDataContainer.end(), iSlaveDof);
        return *(*pos);
    }

    ConstraintEquation &GetConstraintEquation(unsigned int iSlaveEqutionId)
    {
        auto pos = std::find_if(mDataContainer.begin(), mDataContainer.end(), [&iSlaveEqutionId](ConstraintEquationPointerType &obj) { return obj->SlaveEquationId() == iSlaveEqutionId; });
        return *(*pos);
    }

    /**
		Get the Total number of MasterDOFs for a given slave dof
		@return Total number of MasterDOFs for a given slave dof
		 */
    int GetNumbeOfMasterDofsForSlave(DofType const &iSlaveDof)
    {
        auto pos = std::find(mDataContainer.begin(), mDataContainer.end(), iSlaveDof);
        int numMasters = -1;
        if (pos != mDataContainer.end())
        {
            numMasters = (*pos)->NumberOfMasters();
        }
        return numMasters;
    }
    int GetNumbeOfMasterDofsForSlave(unsigned int iSlaveEqutionId)
    {
        // using find_if with a lambda to find the salve with its equation Id
        auto pos = std::find_if(mDataContainer.begin(), mDataContainer.end(), [&iSlaveEqutionId](ConstraintEquationPointerType &obj) { return obj->SlaveEquationId() == iSlaveEqutionId; });
        int numMasters = -1;
        if (pos != mDataContainer.end())
        {
            numMasters = (*pos)->NumberOfMasters();
        }
        return numMasters;
    }

    /**
		Get the Data for this slave
		@return Data vector for this slave
		*/
    ConstraintEquationVectorType &GetData()
    {
        return mDataContainer;
    }

    /**
		Adds a constraints between the given slave and master with a weight. 		
		*/

    // Takes in a slave dof and a master dof
    void AddConstraint(DofType const &iSlaveDof, DofType const &iMasterDof, double iWeight, double constant = 0.0)
    {
        auto pos = std::find(mDataContainer.begin(), mDataContainer.end(), iSlaveDof);
        if (pos != mDataContainer.end())
        { // Equation already exists
            (*pos)->AddMasterData(iMasterDof, iWeight);
        }
        else
        { // Equation does not exist
            ConstraintEquationPointerType newEq = ConstraintEquationPointerType(new ConstraintEquation(iSlaveDof));
            newEq->AddMasterData(iMasterDof, iWeight);
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
