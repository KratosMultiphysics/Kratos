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

#if !defined(CONSTRAINT_H)
#define CONSTRAINT_H
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
#include "geometries/geometry.h"
/* //#include "includes/process_info.h"
#include "includes/element.h"
#include "includes/condition.h" */
#include "containers/constraint_data.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"

namespace Kratos
{
/** \brief Constraint
	* A class that implements the interface for different constraints to be applied on a system.
    */

class ProcessInfo;
class Element;
class Condition;

template <class TDenseSpace>
class Constraint
{

  public:
    /// Pointer definition of DataValueContainer
    KRATOS_CLASS_POINTER_DEFINITION(Constraint);

    typedef Dof<double> DofType;
    typedef std::unordered_map<unsigned int, double> MasterIdWeightMapType;
    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;
    typedef std::vector<std::size_t> EquationIdVectorType;
    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;
    typedef typename ConstraintData::VariableDataType VariableDataType;

    ///@name Life Cycle
    ///@{

    /**
	* Creates a Constraint object
	*/
    Constraint(std::string iName="default", bool iIsActive=true): mName(iName), mActive(iIsActive)
    {
    }
    /// Destructor.
    virtual ~Constraint(){};

    ///@}

    ///@name Access
    ///@{

    /**
	* Clears the maps contents
	*/
    void Clear()
    {
        mConstraintData.Clear();
    }

    /**
	* Get the MasterDOFs vector for this slave
	* @return MasterDOFs vector for this slave
	*/
    virtual const SlaveData &GetSlaveData(DofType &SlaveDof)
    {
        return mConstraintData.GetSlaveData(SlaveDof);
    }

    /**
    * Adds a constraints between the given slave and master with a weight. 		
	*/

    // Takes in a slave dof equationId and a master dof equationId
    virtual void AddConstraint(DofType &SlaveDof, DofType &MasterDof, double weight, double constant = 0.0)
    {
        mConstraintData.AddConstraint(SlaveDof, MasterDof, weight, constant);
    }

    virtual void ExecuteBeforeBuilding(NodesContainerType &Nodes)
    {
    }

    virtual void ExecuteAfterBuilding(NodesContainerType &Nodes)
    {
    }

    virtual void ExecuteBeforeSolving(NodesContainerType &Nodes)
    {
    }

    virtual void ExecuteAfterSolving(NodesContainerType &Nodes)
    {
    }

    virtual void FormulateEquationIdRelationMap(NodesContainerType &Nodes)
    {
    }

    virtual void UpdateConstraintEquationsAfterIteration(NodesContainerType &Nodes)
    {
    }

    virtual void Element_ModifyEquationIdsForConstraints(Element& rCurrentElement,
                                                 EquationIdVectorType &EquationId,
                                                 ProcessInfo &CurrentProcessInfo)
    {
    }

    virtual void Condition_ModifyEquationIdsForConstraints(Condition& rCurrentCondition,
                                                   EquationIdVectorType &EquationId,
                                                   ProcessInfo &CurrentProcessInfo)
    {
    }

    virtual void Element_ApplyConstraints(Element& rCurrentElement,
                                                    LocalSystemMatrixType &LHS_Contribution,
                                                    LocalSystemVectorType &RHS_Contribution,
                                                    EquationIdVectorType &EquationId,
                                                    ProcessInfo &CurrentProcessInfo)
    {
    }

    virtual void Condition_ApplyConstraints(Condition& rCurrentCondition,
                                              LocalSystemMatrixType &LHS_Contribution,
                                              LocalSystemVectorType &RHS_Contribution,
                                              EquationIdVectorType &EquationId,
                                              ProcessInfo &CurrentProcessInfo)
    {
    }

    /**
	* Get the Total number of MasterDOFs for a given slave dof
	* @return Total number of MasterDOFs for a given slave dof
	*/
    virtual unsigned int GetNumbeOfMasterDofsForSlave(const DofType &SlaveDof)
    {
        return mConstraintData.GetNumbeOfMasterDofsForSlave(SlaveDof);
    }

    /**
	* Set the name for the current set of constraints. 
	*/
    virtual void SetName(const std::string name)
    {
        mName = name;
    }
    /**
	* Get the name for the current set of constraints. 
	*/
    virtual std::string GetName()
    {
        return mName;
    }

    /**
	* Set the activeness for current set of constraints. 
	*/
    virtual void SetActive(const bool isActive)
    {
        mActive = isActive;
    }

    /**
	* Returns true if the constraint set is active
	*/
    virtual bool IsActive()
    {
        return mActive;
    }
    /**
	* Returns the string containing a detailed description of this object.
	* @return the string with informations
	*/
    virtual void GetInfo() const
    {
    }

    /**
	* Returns the constraint data object of this constraint object object.
	* @return the string with informations
	*/
    virtual ConstraintData &GetData()
    {
        return mConstraintData;
    }

    ///@

    ///@name Static Operations
    ///
    //@{

    ///@}
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << " Constraint base class !" << std::endl;
        mConstraintData.PrintInfo(rOStream);
    }

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer &rSerializer) const
    {
        rSerializer.save("MpcDataName", mName);
        rSerializer.save("isActive", mActive);
        this->mConstraintData.save(rSerializer);
    }

    virtual void load(Serializer &rSerializer)
    {
        rSerializer.load("MpcDataName", mName);
        rSerializer.load("isActive", mActive);
        this->mConstraintData.load(rSerializer);
    }

  private:
    ///@name Member Variables
    ///@{
    bool mActive;
    std::string mName;

    ConstraintData mConstraintData;
    ///@}

    ///@}
};

///@name Input/Output funcitons
///@{

///@}

} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
