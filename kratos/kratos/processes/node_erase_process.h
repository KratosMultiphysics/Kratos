//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//



#if !defined(KRATOS_NODE_ERASE_PROCESS_INCLUDED )
#define  KRATOS_NODE_ERASE_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes
#include "includes/kratos_flags.h"



// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"


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
//erases the nodes marked as
/** Detail class definition.
*/

class NodeEraseProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NodeEraseProcess
    KRATOS_CLASS_POINTER_DEFINITION(NodeEraseProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NodeEraseProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~NodeEraseProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void Execute()
    {
        KRATOS_TRY;

        mr_model_part.RemoveNodes(TO_ERASE);

        KRATOS_CATCH("")
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
        return "NodeEraseProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NodeEraseProcess";
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


    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mr_model_part;
    PointerVector<Node<3> > mTrashedNodes;


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
    NodeEraseProcess& operator=(NodeEraseProcess const& rOther);

    /// Copy constructor.
    //NodeEraseProcess(NodeEraseProcess const& rOther);


    ///@}

}; // Class NodeEraseProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  NodeEraseProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const NodeEraseProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_NODE_ERASE_PROCESS_INCLUDED  defined 


