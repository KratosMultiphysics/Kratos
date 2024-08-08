//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//



#if !defined(KRATOS_ELIMINATE_ISOLATED_NODES_PROCESS_INCLUDED )
#define  KRATOS_ELIMINATE_ISOLATED_NODES_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
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

class EliminateIsolatedNodesProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EliminateIsolatedNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(EliminateIsolatedNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    EliminateIsolatedNodesProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /// Destructor.
    ~EliminateIsolatedNodesProcess() override
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

    void Execute() override
    {
        KRATOS_TRY;

        ModelPart::NodesContainerType temp_nodes_container;
        temp_nodes_container.reserve(mr_model_part.Nodes().size());

        temp_nodes_container.swap(mr_model_part.Nodes());

        for(ModelPart::NodesContainerType::iterator i_node = temp_nodes_container.begin() ; i_node != temp_nodes_container.end() ; i_node++)
        {
            if( i_node->GetValue(NEIGHBOUR_NODES).size() != 0 )
                (mr_model_part.Nodes()).insert(mr_model_part.Nodes().end(), *(i_node.base()));
        }

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
    std::string Info() const override
    {
        return "EliminateIsolatedNodesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EliminateIsolatedNodesProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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
    PointerVector<Node > mTrashedNodes;


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
    EliminateIsolatedNodesProcess& operator=(EliminateIsolatedNodesProcess const& rOther);

    /// Copy constructor.
    //EliminateIsolatedNodesProcess(EliminateIsolatedNodesProcess const& rOther);


    ///@}

}; // Class EliminateIsolatedNodesProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  EliminateIsolatedNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const EliminateIsolatedNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ELIMINATE_ISOLATED_NODES_PROCESS_INCLUDED  defined


