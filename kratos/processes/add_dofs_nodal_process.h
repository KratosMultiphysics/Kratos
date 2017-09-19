//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//


#if !defined(KRATOS_ADD_DOFS_NODAL_PROCESS_H_INCLUDED )
#define  KRATOS_ADD_DOFS_NODAL_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "containers/pointer_vector_set.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
template<class TVariableType>
class AddDofsNodalProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AddDofsNodalProcess
    KRATOS_CLASS_POINTER_DEFINITION(AddDofsNodalProcess);

    typedef Node<3> NodeType;

    typedef Dof<typename TVariableType::Type> DofType;

    typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;

    typedef PointerVectorSet<DofType, IndexedObject> DofsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    AddDofsNodalProcess(TVariableType const& NewVariable, NodesContainerType & rThisNodes, DofsContainerType& rThisDofs)
        : mVariable(NewVariable), mNodes(rThisNodes), mDofs(rThisDofs) {}

    /// Destructor.
    virtual ~AddDofsNodalProcess() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Execute()
    {
        for(typename NodesContainerType::iterator i_node = mNodes.begin() ; i_node != mNodes.end() ; i_node++)
            mDofs.push_back(i_node->pAddDof(mVariable));
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "add "
               << mVariable.Name()
               << " dofs nodal process";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{
    TVariableType const& mVariable;

    NodesContainerType& mNodes;

    DofsContainerType& mDofs;

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    AddDofsNodalProcess& operator=(AddDofsNodalProcess const& rOther);

    /// Copy constructor.
    AddDofsNodalProcess(AddDofsNodalProcess const& rOther);


    ///@}

}; // Class AddDofsNodalProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_ADD_DOFS_NODAL_PROCESS_H_INCLUDED  defined 


