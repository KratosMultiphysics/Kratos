//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Pavel Ryzhakov, Daniel Baumgärtner, Johannes Wolf $
//   Date:                $Date: 2013-07-15 17:20:00 $
//   Revision:            $Revision: 1.0 $
//
//  this process defines functions needed for the kratos empire interface

#if !defined(WRAPPER_PROCESS_INCLUDED )
#define  WRAPPER_PROCESS_INCLUDED

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
#include "geometries/triangle_3d_3.h"

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

class WrapperProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of WrapperProcess
    KRATOS_CLASS_POINTER_DEFINITION(WrapperProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    WrapperProcess(ModelPart& model_part, ModelPart& interface_part)
        : mr_model_part(model_part), mr_interface_part(interface_part)
    {
    }

    /// Destructor.
    virtual ~WrapperProcess()
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
        KRATOS_TRY

        // Add interface nodes (nodes of fluid model part with flag IS_INTERFACE)
        // to interface model part
        for( ModelPart::NodeIterator i_node =  mr_model_part.NodesBegin() ;
                                     i_node != mr_model_part.NodesEnd() ;
                                     i_node++ )
        {
            if( i_node->FastGetSolutionStepValue(IS_INTERFACE) == 1.0 )
            {
                mr_interface_part.Nodes().push_back( *(i_node.base()) );
            }
        }

        // Add interface conditions
        for ( ModelPart::ConditionIterator i_condition =  mr_model_part.ConditionsBegin();
                                           i_condition != mr_model_part.ConditionsEnd();
                                           i_condition++ )
        {
            int size = (*i_condition).GetGeometry().size();

            if(size == 2) // interface conditions are lines (2D fluid mesh)
            {
                if ( ((*i_condition).GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE) == 1.0) &&
                     ((*i_condition).GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE) == 1.0))
                {
                    mr_interface_part.Conditions().push_back( *(i_condition.base()) );
                }
            }
            else if(size == 3) // interface conditions are triangles (3D fluid mesh)
            {
                if ( ((*i_condition).GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE) == 1.0) &&
                     ((*i_condition).GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE) == 1.0) &&
                     ((*i_condition).GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE) == 1.0))
                {
//                    NECESSARY FOR CANTILEVER EXAMPLE
                    double Y0 = (*i_condition).GetGeometry()[0].Y();
                    double Y1 = (*i_condition).GetGeometry()[1].Y();
                    double Y2 = (*i_condition).GetGeometry()[2].Y();

                    if(Y0 > 0 || Y1 > 0 || Y2 > 0 )
                    {
                        mr_interface_part.Conditions().push_back( *(i_condition.base()) );
                    }
                    // NECESSARY FOR TUREK EXAMPLE
//                    double Z0 = (*i_condition).GetGeometry()[0].Z();
//                    double Z1 = (*i_condition).GetGeometry()[1].Z();
//                    double Z2 = (*i_condition).GetGeometry()[2].Z();

//                    if((Z0 < 0 || Z1 < 0 || Z2 < 0) && (Z0 > -0.01 || Z1 > -0.01 || Z2 > -0.01) )
//                    {
//                        mr_interface_part.Conditions().push_back( *(i_condition.base()) );
//                    }
                }
            }

        }

        KRATOS_CATCH("")
    }

    void ExtractPressureFromModelPart( boost::python::list& pressure )
    {
        KRATOS_TRY

        for( ModelPart::NodeIterator i_node =  mr_interface_part.NodesBegin() ;
                                     i_node != mr_interface_part.NodesEnd() ;
                                     i_node++ )
        {
                double p = i_node->GetSolutionStepValue(PRESSURE);

                pressure.append(p);
        }

        KRATOS_CATCH("")
    }

    void ExtractForcesFromModelPart( boost::python::list& forces )
    {
        KRATOS_TRY

        for( ModelPart::NodeIterator i_node =  mr_interface_part.NodesBegin() ;
                                     i_node != mr_interface_part.NodesEnd() ;
                                     i_node++ )
        {
            double r_x = i_node->GetSolutionStepValue(REACTION_X);
            double r_y = i_node->GetSolutionStepValue(REACTION_Y);
            double r_z = i_node->GetSolutionStepValue(REACTION_Z);

            // Negative of reactions = forces to structure
            double f_x = -r_x;
            double f_y = -r_y;
            double f_z = -r_z;


            forces.append(f_x);
            forces.append(f_y);
            forces.append(f_z);
        }

        KRATOS_CATCH("")
    }

    // ##############################################################################
    // Function required for Kratos-Kratos-FSI
    void ExtractDisplacementsFromModelPart( boost::python::list& displacements )
    {
        KRATOS_TRY

        for( ModelPart::NodeIterator i_node =  mr_interface_part.NodesBegin() ;
                                     i_node != mr_interface_part.NodesEnd() ;
                                     i_node++ )
        {
                array_1d<double,3> vect_disp = i_node->GetSolutionStepValue(DISPLACEMENT);

                displacements.append(vect_disp[0]);
                displacements.append(vect_disp[1]);
                displacements.append(vect_disp[2]);
        }

        KRATOS_CATCH("")
    }

    // ##############################################################################

    void CreateEmbeddedStructurePart( boost::python::list nodes , boost::python::list connectivity )
    {
        KRATOS_TRY

        ModelPart& new_model_part = mr_interface_part;

        // ######## ADDING NEW NODES #########

        const unsigned int size_nodes = boost::python::len(nodes) / 3;
        new_model_part.Nodes().reserve(size_nodes);

        for (unsigned int nodesIndex=0; nodesIndex!=size_nodes; ++nodesIndex)
        {
            unsigned int id = nodesIndex*3 + 1;
            boost::python::extract<double> node_X( nodes[id] );
            boost::python::extract<double> node_Y( nodes[id+1] );
            boost::python::extract<double> node_Z( nodes[id+2] );

            Node < 3 >::Pointer pnode = Node < 3 > ::Pointer (new Node < 3 >(id,
                                                                             node_X,
                                                                             node_Y,
                                                                             node_Z));
            new_model_part.Nodes().push_back(pnode);
        }

        // ######## ADDING NEW CONDITIONS #########
        // required information for a new condition: Id, Geometry, Property
        const unsigned int size_conditions = boost::python::len(connectivity) / 3;

        Condition const& rReferenceCondition = KratosComponents<Condition>::Get("Condition3D");
        Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(1);
        new_model_part.Conditions().reserve(size_conditions);

        //#pragma omp parallel for
        for (unsigned int condIndex=0; condIndex!=size_conditions; ++condIndex)
        {
            unsigned int id = condIndex*3 + 1;

            boost::python::extract<unsigned int> node1_ID( connectivity[id] );
            boost::python::extract<unsigned int> node2_ID( connectivity[id+1] );
            boost::python::extract<unsigned int> node3_ID( connectivity[id+2] );

            Triangle3D3<Node < 3 > > triangle(new_model_part.Nodes()(node1_ID),
                                              new_model_part.Nodes()(node2_ID),
                                              new_model_part.Nodes()(node3_ID));
            // Index must take care of the IDs of the other processors --> triangle_id_int + 1 + elements_before + total_existing_elements
            Condition::Pointer p_condition = rReferenceCondition.Create(id,
                                                                        triangle,
                                                                        properties);
            new_model_part.Conditions().push_back(p_condition);
        }

        KRATOS_CATCH("")
    }

    // ##############################################################################

    void ExtractMeshInfo( boost::python::list& numNodes, boost::python::list& numElems,
                          boost::python::list& nodes, boost::python::list& nodeIDs,
                          boost::python::list& numNodesPerElem, boost::python::list& elems )
    {
        KRATOS_TRY

        unsigned int nodesCounter = 0;
        unsigned int elemsCounter = 0;

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

            nodesCounter++;

            // Fill the nodeIDs vector with the nodal IDs
            nodeIDs.append(i_fluidNode->Id());
        }

        // loop over all fluid elements
        for( ModelPart::ConditionsContainerType::iterator i_fluidCondition =  mr_interface_part.ConditionsBegin();
                                                          i_fluidCondition != mr_interface_part.ConditionsEnd();
                                                          i_fluidCondition++)
        {
            unsigned int nodesPerElem = i_fluidCondition->GetGeometry().size();
            numNodesPerElem.append(nodesPerElem);

            if( nodesPerElem == 2 ) // interface conditions are lines (2D fluid mesh)
            {
                unsigned int nodeID_1 = i_fluidCondition->GetGeometry()[0].Id();
                unsigned int nodeID_2 = i_fluidCondition->GetGeometry()[1].Id();

                elems.append(nodeID_1);
                elems.append(nodeID_2);

                elemsCounter++;
            }
            else if( nodesPerElem == 3 ) // interface conditions are triangles (3D fluid mesh)
            {
                unsigned int nodeID_1 = i_fluidCondition->GetGeometry()[0].Id();
                unsigned int nodeID_2 = i_fluidCondition->GetGeometry()[1].Id();
                unsigned int nodeID_3 = i_fluidCondition->GetGeometry()[2].Id();

                elems.append(nodeID_1);
                elems.append(nodeID_2);
                elems.append(nodeID_3);

                elemsCounter++;
            }
        }

        numNodes.append(nodesCounter);
        numElems.append(elemsCounter);

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
        return "WrapperProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "WrapperProcess";
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
//		WrapperProcess& operator=(WrapperProcess const& rOther);

    /// Copy constructor.
//		WrapperProcess(WrapperProcess const& rOther);


    ///@}

}; // Class WrapperProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/*
/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  WrapperProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const WrapperProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
*/
///@}


}  // namespace Kratos.

#endif // WRAPPER_PROCESS_INCLUDED  defined
