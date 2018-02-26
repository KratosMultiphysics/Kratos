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
     *  Does necessary operations to setup the constraint.
     */
    virtual void SetUp(NodesContainerType &rNodes)
    {
    }

    /**
     *  Does necessary operations on the constraint before the build of master stiffness matrix is commenced
     */
    virtual void ExecuteBeforeBuilding(NodesContainerType &rNodes)
    {
    }
    /**
     *  Does necessary operations on the constraint after the build of master stiffness matrix is commenced
     */
    virtual void ExecuteAfterBuilding(NodesContainerType &rNodes)
    {
    }

    /**
     *  Does necessary operations on global symtem Ax=b on the constraint before they are solved.
     */    
    virtual void ExecuteBeforeSolving(TSystemMatrixType &rA,
                                      TSystemVectorType &rDx,
                                      TSystemVectorType &rb)
    {
    }

    /**
     *  Does necessary operations on global symtem Ax=b on the constraint after they are solved.
     */    
    virtual void ExecuteAfterSolving(TSystemMatrixType &rA,
                                     TSystemVectorType &rDx,
                                     TSystemVectorType &rb)
    {
    }

    /**
     *  Does necessary operations on the element freedom table so as to construct the global stiffness matrix.
     *  Mainly to build the sparsity pattern for the global stiffness matrix.
     */        
    virtual void ModifyEquationIdsForConstraints(Element &rCurrentElement,
                                                         EquationIdVectorType &rEquationId,
                                                         ProcessInfo &rCurrentProcessInfo)
    {
    }

    /**
     *  Does necessary operations on the element freedom table of condition so as to construct the global stiffness matrix.
     *  Mainly to build the sparsity pattern for the global stiffness matrix.
     */            
    virtual void ModifyEquationIdsForConstraints(Condition &rCurrentCondition,
                                                           EquationIdVectorType &rEquationId,
                                                           ProcessInfo &rCurrentProcessInfo)
    {
    }

    /**
     *  Does necessary operations on the element stiffness matrix and element rhs to apply the constraints.
     */    
    virtual void ApplyConstraints(Element &rCurrentElement,
                                          LocalSystemMatrixType &rLHS_Contribution,
                                          LocalSystemVectorType &rRHS_Contribution,
                                          EquationIdVectorType &rEquationId,
                                          ProcessInfo &rCurrentProcessInfo)
    {
    }

    /**
     *  Does necessary operations on the condition stiffness matrix and element rhs to apply the constraints.
     */        
    virtual void ApplyConstraints(Condition &rCurrentCondition,
                                            LocalSystemMatrixType &rLHS_Contribution,
                                            LocalSystemVectorType &rRHS_Contribution,
                                            EquationIdVectorType &rEquationId,
                                            ProcessInfo &rCurrentProcessInfo)
    {
    }

    /**
	* Get the Total number of MasterDOFs for a given slave dof
	* @return Total number of MasterDOFs for a given slave dof
	*/
    virtual unsigned int GetNumberOfMasterDofsForSlave(DofType &rSlaveDof)
    {
        return mConstraintEquationContainer.GetNumberOfMasterDofsForSlave(rSlaveDof);
    }

    /**
	* Set the name for the current set of constraints. 
	*/
    virtual void SetName(const std::string& Name)
    {
        mName = Name;
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
    virtual void SetActive(const bool IsActive)
    {
        mActive = IsActive;
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
