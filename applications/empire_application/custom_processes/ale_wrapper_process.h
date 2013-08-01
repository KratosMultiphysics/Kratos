//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Pavel Ryzhakov, Daniel Baumgärtner, Johannes Wolf $
//   Date:                $Date: 2013-07-15 17:20:00 $
//   Revision:            $Revision: 1.0 $
//
//  this process defines functions needed for the kratos empire interface

#if !defined(ALE_WRAPPER_PROCESS_INCLUDED )
#define  ALE_WRAPPER_PROCESS_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <algorithm>

// External includes
#include <boost/python.hpp>

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
/** Detail class definition.

*/

class ALEWrapperProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ALEWrapperProcess
    KRATOS_CLASS_POINTER_DEFINITION(ALEWrapperProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ALEWrapperProcess(ModelPart& model_part, ModelPart& interface_part)
        : mr_model_part(model_part), mr_interface_part(interface_part)
    {
    }

    /// Destructor.
    virtual ~ALEWrapperProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ExtractInterface()
    {
        // Add interface nodes to mesh no. 1 of fluid model part
        for(ModelPart::NodeIterator i_node =  mr_model_part.NodesBegin() ;
                                    i_node != mr_model_part.NodesEnd() ;
                                    i_node++)
        {
            if( i_node->FastGetSolutionStepValue(IS_INTERFACE) == 1.0 )
            {
                mr_interface_part.Nodes().push_back( *(i_node.base()) );
            }
        }

        // Generate triangular Conditions from original interface conditions
        for ( ModelPart::ConditionIterator i_condition =  mr_model_part.ConditionsBegin();
                                           i_condition != mr_model_part.ConditionsEnd();
                                           i_condition++)
        {
            if ( ((*i_condition).GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE) == 1.0) &&
                 ((*i_condition).GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE) == 1.0) &&
                 ((*i_condition).GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE) == 1.0))
            {
                mr_interface_part.Conditions().push_back( *(i_condition.base()) );
            }
        }
    }

    void ExtractForcesFromModelPart( boost::python::list& reactions )
    {
        KRATOS_TRY

        for(ModelPart::NodeIterator i_node =  mr_interface_part.NodesBegin() ;
                                    i_node != mr_interface_part.NodesEnd() ;
                                    i_node++)
        {
                double reactionX = i_node->GetSolutionStepValue(REACTION_X);
                double reactionY = i_node->GetSolutionStepValue(REACTION_Y);
                double reactionZ = i_node->GetSolutionStepValue(REACTION_Z);

                reactions.append(reactionX);
                reactions.append(reactionY);
                reactions.append(reactionZ);
        }

        KRATOS_CATCH("")
    }

    void ExtractMeshInfo(boost::python::list& numNodes, boost::python::list& numElems,
                         boost::python::list& nodes, boost::python::list& nodeIDs,
                         boost::python::list& numNodesPerElem, boost::python::list& elems)
    {
        numNodes.append(mr_interface_part.NodesArray().size());
        numElems.append(mr_interface_part.ConditionsArray().size());

        // loop over all fluid nodes
        for ( ModelPart::NodesContainerType::iterator i_fluidNode =  mr_interface_part.NodesBegin();
                                                      i_fluidNode != mr_interface_part.NodesEnd();
                                                      ++i_fluidNode )
        {
            double node_X = i_fluidNode->Coordinates()[0];
            double node_Y = i_fluidNode->Coordinates()[1];
            double node_Z = i_fluidNode->Coordinates()[2];

            // Fill the nodes vector with nodal coordinates
            nodes.append(node_X);
            nodes.append(node_Y);
            nodes.append(node_Z);

            // Fill the nodeIDs vector with the nodal IDs
            nodeIDs.append(i_fluidNode->Id());
        }

        const unsigned int nodesPerElem = 3;

        // loop over all fluid elements
        for( ModelPart::ConditionsContainerType::iterator i_fluidCondition =  mr_interface_part.ConditionsBegin();
                                                          i_fluidCondition != mr_interface_part.ConditionsEnd();
                                                          i_fluidCondition++)
        {
            numNodesPerElem.append(nodesPerElem);

            unsigned int nodeID_1 = i_fluidCondition->GetGeometry()[0].Id();
            unsigned int nodeID_2 = i_fluidCondition->GetGeometry()[1].Id();
            unsigned int nodeID_3 = i_fluidCondition->GetGeometry()[2].Id();

            elems.append(nodeID_1);
            elems.append(nodeID_2);
            elems.append(nodeID_3);
        }
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
        return "ALEWrapperProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ALEWrapperProcess";
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
    ModelPart& mr_interface_part;

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
//		ALEWrapperProcess& operator=(ALEWrapperProcess const& rOther);

    /// Copy constructor.
//		ALEWrapperProcess(ALEWrapperProcess const& rOther);


    ///@}

}; // Class ALEWrapperProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/*
/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ALEWrapperProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ALEWrapperProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
*/
///@}


}  // namespace Kratos.

#endif // ALE_WRAPPER_PROCESS_INCLUDED  defined


