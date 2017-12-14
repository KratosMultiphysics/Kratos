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

template <class TSparseSpace,
          class TDenseSpace //= DenseSpace<double>
          >
class Constraint
{

  public:
    /// Pointer definition of DataValueContainer
    KRATOS_CLASS_POINTER_DEFINITION(Constraint);

    typedef Dof<double> DofType;
    typedef Node<3> NodeType;
    typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;
    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;
    typedef typename TSparseSpace::DataType TDataType;
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    ///@name Life Cycle
    ///@{

    /**
	* Creates a Constraint object
	*/
    Constraint(std::string iName = "default", bool iIsActive = true) : mActive(iIsActive), mName(iName)
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
        mConstraintEquationContainer.Clear();
    }

    /**
	* Get the MasterDOFs vector for this slave
	* @return MasterDOFs vector for this slave
	*/
    virtual const ConstraintEquation &GetConstraintEquation(const DofType &SlaveDof)
    {
        return mConstraintEquationContainer.GetConstraintEquation(SlaveDof);
    }

    /**
    * Adds a constraints between the given slave and master with a weight. 		
	*/
    // Takes in a slave dof equationId and a master dof equationId
    virtual void AddConstraint(DofType &SlaveDof, DofType &MasterDof, double weight, double constant = 0.0)
    {
        mConstraintEquationContainer.AddConstraint(SlaveDof, MasterDof, weight, constant);
    }

    /**
     *  Does necessary operations to setup the constraint.
     */
    virtual void SetUp(NodesContainerType &Nodes)
    {
    }

    /**
     *  Does necessary operations on the constraint before the build of master stiffness matrix is commenced
     */
    virtual void ExecuteBeforeBuilding(NodesContainerType &Nodes)
    {
    }
    /**
     *  Does necessary operations on the constraint after the build of master stiffness matrix is commenced
     */
    virtual void ExecuteAfterBuilding(NodesContainerType &Nodes)
    {
    }

    /**
     *  Does necessary operations on global symtem Ax=b on the constraint before they are solved.
     */    
    virtual void ExecuteBeforeSolving(TSystemMatrixType &A,
                                      TSystemVectorType &Dx,
                                      TSystemVectorType &b)
    {
    }

    /**
     *  Does necessary operations on global symtem Ax=b on the constraint after they are solved.
     */    
    virtual void ExecuteAfterSolving(TSystemMatrixType &A,
                                     TSystemVectorType &Dx,
                                     TSystemVectorType &b)
    {
    }

    /**
     *  Does necessary operations on the element freedom table so as to construct the global stiffness matrix.
     *  Mainly to build the sparsity pattern for the global stiffness matrix.
     */        
    virtual void Element_ModifyEquationIdsForConstraints(Element &rCurrentElement,
                                                         EquationIdVectorType &EquationId,
                                                         ProcessInfo &CurrentProcessInfo)
    {
    }

    /**
     *  Does necessary operations on the element freedom table of condition so as to construct the global stiffness matrix.
     *  Mainly to build the sparsity pattern for the global stiffness matrix.
     */            
    virtual void Condition_ModifyEquationIdsForConstraints(Condition &rCurrentCondition,
                                                           EquationIdVectorType &EquationId,
                                                           ProcessInfo &CurrentProcessInfo)
    {
    }

    /**
     *  Does necessary operations on the element stiffness matrix and element rhs to apply the constraints.
     */    
    virtual void Element_ApplyConstraints(Element &rCurrentElement,
                                          LocalSystemMatrixType &LHS_Contribution,
                                          LocalSystemVectorType &RHS_Contribution,
                                          EquationIdVectorType &EquationId,
                                          ProcessInfo &CurrentProcessInfo)
    {
    }

    /**
     *  Does necessary operations on the condition stiffness matrix and element rhs to apply the constraints.
     */        
    virtual void Condition_ApplyConstraints(Condition &rCurrentCondition,
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
        return mConstraintEquationContainer.GetNumbeOfMasterDofsForSlave(SlaveDof);
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
    virtual ConstraintEquationContainer &GetData()
    {
        return mConstraintEquationContainer;
    }

    ///@

    ///@name Static Operations
    ///
    //@{

    ///@}
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << " Constraint base class !" << std::endl;
        mConstraintEquationContainer.PrintInfo(rOStream);
    }

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer &rSerializer) const
    {
        rSerializer.save("ConstraintName", mName);
        rSerializer.save("isActive", mActive);
        this->mConstraintEquationContainer.save(rSerializer);
    }

    virtual void load(Serializer &rSerializer)
    {
        rSerializer.load("ConstraintName", mName);
        rSerializer.load("isActive", mActive);
        this->mConstraintEquationContainer.load(rSerializer);
    }

  private:
    ///@name Member Variables
    ///@{
    bool mActive;
    std::string mName;

    ConstraintEquationContainer mConstraintEquationContainer;
    ///@}

    ///@}
};

///@name Input/Output funcitons
///@{

///@}

} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
