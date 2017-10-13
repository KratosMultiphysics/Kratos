/*
==============================================================================
KratosPFEMApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-22 17:13:57 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_NODE_AND_ELEMENT_ERASE_PROCESS_INCLUDED )
#define  KRATOS_NODE_AND_ELEMENT_ERASE_PROCESS_INCLUDED



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
#include "includes/kratos_flags.h"

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

class NodeAndElementEraseProcess
        : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NodeAndElementEraseProcess
    KRATOS_CLASS_POINTER_DEFINITION(NodeAndElementEraseProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NodeAndElementEraseProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
        KRATOS_TRY
                KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~NodeAndElementEraseProcess()
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

        // I add this part, needs checking

        ModelPart::ElementsContainerType temp_elements_container;
        temp_elements_container.reserve(mr_model_part.Elements().size());
        temp_elements_container.swap(mr_model_part.Elements());

        ModelPart::NodesContainerType temp_nodes_container;
        temp_nodes_container.reserve(mr_model_part.Nodes().size());
        temp_nodes_container.swap(mr_model_part.Nodes());

        mr_model_part.Elements().clear(); // ??
        mr_model_part.Nodes().clear();




        unsigned int setid = 1;

        for(ModelPart::ElementsContainerType::iterator i_elem = temp_elements_container.begin() ; i_elem != temp_elements_container.end() ; i_elem++)
        {
            if( static_cast<bool>(i_elem->GetGeometry()(0)->Is(TO_ERASE)) == false){

                i_elem->SetId( i_elem->GetGeometry()(0)->Id() );

//                i_elem->SetId( setid ); //****

                (mr_model_part.Elements()).push_back(*(i_elem.base()));

                ModelPart::NodeType::Pointer i_node = (*i_elem).GetGeometry()(0) ;

//                i_node->SetId(setid); //****

                //(*particle_pointer_it)->GetGeometry().pGetPoint(i)
                (mr_model_part.Nodes()).push_back(i_node);

                (i_elem->GetValue(NEIGHBOUR_ELEMENTS)).clear();

                setid++;


            }
        }

        KRATOS_WATCH(setid);




        KRATOS_CATCH("");
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
        return "NodeAndElementEraseProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NodeAndElementEraseProcess";
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
    NodeAndElementEraseProcess& operator=(NodeAndElementEraseProcess const& rOther);

    /// Copy constructor.
    //NodeAndElementEraseProcess(NodeAndElementEraseProcess const& rOther);


    ///@}

}; // Class NodeAndElementEraseProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  NodeAndElementEraseProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const NodeAndElementEraseProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_NODE_AND_ELEMENT_ERASE_PROCESS_INCLUDED   defined

