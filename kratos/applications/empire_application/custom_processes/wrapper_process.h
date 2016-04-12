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
#include "includes/kratos_flags.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "../structural_application/custom_elements/membrane_element.h"

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

    //definition of node type
    typedef ModelPart::NodeType NodeType;

    typedef NodeType::IndexType IndexType;

    ///definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;

    /// Nodes container. Which is a vector set of nodes with their Id's as key.
    typedef ModelPart::NodesContainerType NodesContainerType;

    /// Nodes map.
    typedef PointerVectorMap<IndexType, NodeType> NodesContainerMapType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    WrapperProcess(ModelPart& model_part, ModelPart& interface_part, int dimension)
        : mr_model_part(model_part),
          mr_interface_part(interface_part),
          m_dimension(dimension)
    {
      KRATOS_TRY
	if (m_dimension != 2 && m_dimension != 3)
	    KRATOS_THROW_ERROR(std::logic_error, "invalid dimension", "");

      mr_interface_part.GetNodalSolutionStepVariablesList() = mr_model_part.GetNodalSolutionStepVariablesList();
      mr_interface_part.SetBufferSize(mr_model_part.GetBufferSize());

      KRATOS_CATCH("")
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

	if (m_dimension == 2)
	  {
	    mr_model_part.Nodes().Sort();

	    IndexType NodeId = mr_model_part.Nodes().back().Id();

	    mDummyNodes2D.clear();

	    // EMPIRE requires a 3D mesh so this 2D adapter creates a set of 
	    // dummy nodes with z-offset in third dimension. 
	    for( ModelPart::NodeIterator i_node =  mr_model_part.NodesBegin() ;
		 i_node != mr_model_part.NodesEnd() ;
		 i_node++ )
	      {
		if ((*i_node).Is(INTERFACE))
		  {
		    double NodeX = (*i_node).X0();
		    double NodeY = (*i_node).Y0();
		    double NodeZ = (*i_node).Z0() + 0.1; // z-offset for dummyNodes2D (must be same on structure side)
		    mr_interface_part.Nodes().push_back( *(i_node.base()) );
		    NodeType::Pointer pNode = mr_interface_part.CreateNewNode(++NodeId,NodeX,NodeY,NodeZ);
		    mDummyNodes2D.push_back( NodesContainerMapType::value_type((*i_node).Id(),pNode) );
		  }
	      }

	    mr_model_part.Conditions().Sort();
	    
	    IndexType ConditionId = mr_model_part.Conditions().back().Id();

	    // The 2D adapter creates quadrilateral elements to map the interface data. 
	    for ( ModelPart::ConditionIterator i_condition =  mr_model_part.ConditionsBegin();
		  i_condition != mr_model_part.ConditionsEnd();
		  i_condition++ )
	      {
		if ((*i_condition).Is(INTERFACE))
		  {
                    typename Node<3>::Pointer pNode1 = i_condition->GetGeometry()(0);
                    typename Node<3>::Pointer pNode2 = i_condition->GetGeometry()(1);
		    // make sure the node was correctly set in the mDummyNodes2D map
		    if (pNode1->IsNot(INTERFACE) || pNode2->IsNot(INTERFACE))
		      KRATOS_THROW_ERROR(std::logic_error, "A condition node is not INTERFACE","");
	    
		    typename Node<3>::Pointer pNode3 = mDummyNodes2D(pNode2->Id());
		    typename Node<3>::Pointer pNode4 = mDummyNodes2D(pNode1->Id());
		    GeometryType::Pointer pGeom(new Quadrilateral3D4<Node<3> >(pNode1,pNode2,pNode3,pNode4));
                    mr_interface_part.Conditions().push_back( Condition::Pointer(new Condition(++ConditionId,pGeom)) );
		  }
	      }
	  }
	else if (m_dimension == 3)
	  {
	    for( ModelPart::NodeIterator i_node =  mr_model_part.NodesBegin() ;
		 i_node != mr_model_part.NodesEnd() ;
		 i_node++ )
	      {
		if ((*i_node).Is(INTERFACE))
		  mr_interface_part.Nodes().push_back( *(i_node.base()) );
	      }
	    
	    for ( ModelPart::ConditionIterator i_condition =  mr_model_part.ConditionsBegin();
		  i_condition != mr_model_part.ConditionsEnd();
		  i_condition++ )
	      {
		if((*i_condition).Is(INTERFACE))
		    mr_interface_part.Conditions().push_back( *(i_condition.base()) );
	      }
	  }

        KRATOS_CATCH("")
    }

    void ExtractPressureFromModelPart( boost::python::list& pressure )
    {
        KRATOS_TRY

	// the dummyNodes2D get copies of data on their parent nodes before mapping.
	if (m_dimension == 2)
	  SynchronizeDummy2DVariable(PRESSURE);

        for( ModelPart::NodeIterator i_node =  mr_interface_part.NodesBegin() ;
                                     i_node != mr_interface_part.NodesEnd() ;
                                     i_node++ )
        {
                double p = i_node->GetSolutionStepValue(PRESSURE);

                pressure.append(p);
        }

        KRATOS_CATCH("")
    }

    void ExtractPressureFromEmbeddedModelPart( boost::python::list& pressure )
    {
        KRATOS_TRY

	if (m_dimension == 2)
	  KRATOS_THROW_ERROR(std::logic_error, "Not implemented for 2D","");

        for( ModelPart::NodeIterator i_node =  mr_interface_part.NodesBegin() ;
                                     i_node != mr_interface_part.NodesEnd() ;
                                     i_node++ )
        {
                double p_pos = i_node->GetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                double p_neg = i_node->GetSolutionStepValue(NEGATIVE_FACE_PRESSURE);

                double p_diff = p_pos - p_neg;
                
                pressure.append(p_diff);
        }

        KRATOS_CATCH("")
    }

    // ##############################################################################

    void ExtractForcesFromModelPart( boost::python::list& forces )
    {
        KRATOS_TRY
	
	  double weight;

	// the dummyNodes2D get copies of data on their parent nodes before mapping.
	if (m_dimension == 2)
	  {
	    SynchronizeDummy2DVariable(REACTION);
	    // weight is 1.0 for 2D-2D coupling and 0.5 * depth for 2D-3D coupling 
	    // with depth the domain length in z-direction
	    weight = 0.5 * 0.1;
	  }
	else
	  weight = 1.0;

        for( ModelPart::NodeIterator i_node =  mr_interface_part.NodesBegin() ;
                                     i_node != mr_interface_part.NodesEnd() ;
                                     i_node++ )
        {
            double r_x = i_node->GetSolutionStepValue(REACTION_X);
            double r_y = i_node->GetSolutionStepValue(REACTION_Y);
            double r_z = i_node->GetSolutionStepValue(REACTION_Z);

            // Negative of reactions = forces to structure
            double f_x = -weight * r_x;
            double f_y = -weight * r_y;
            double f_z = -weight * r_z;

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

	// the dummyNodes2D get copies of data on their parent nodes before mapping.
	if (m_dimension == 2)
	  SynchronizeDummy2DVariable(DISPLACEMENT);

        for( ModelPart::NodeIterator i_node =  mr_interface_part.NodesBegin() ;
                                     i_node != mr_interface_part.NodesEnd() ;
                                     i_node++ )
        {
                array_1d<double,3> disp = i_node->GetSolutionStepValue(DISPLACEMENT,0);

                displacements.append(disp[0]);
                displacements.append(disp[1]);
                displacements.append(disp[2]);
        }

        KRATOS_CATCH("")
    }

    // ##############################################################################

    void CreateEmbeddedInterfacePart( boost::python::list nodeIDs ,
                                      boost::python::list nodeCoordinates,
                                      boost::python::list connectivity,
                                      boost::python::list displacement,
                                      boost::python::list velocity)
    {
        KRATOS_TRY

	if (m_dimension == 2)
	  KRATOS_THROW_ERROR(std::logic_error, "Not implemented for 2D","");
        
        mr_interface_part.Nodes().erase(mr_interface_part.Nodes().begin(), mr_interface_part.Nodes().end());
        mr_interface_part.Elements().erase(mr_interface_part.Elements().begin(), mr_interface_part.Elements().end());

        // Adding new nodes
        const unsigned int size_nodes = boost::python::len(nodeIDs);
        mr_interface_part.Nodes().reserve(size_nodes);

        for (unsigned int nodesIndex=0; nodesIndex!=size_nodes; ++nodesIndex)
        {
            boost::python::extract<double> node_X( nodeCoordinates[3*nodesIndex] );
            boost::python::extract<double> node_Y( nodeCoordinates[3*nodesIndex+1] );
            boost::python::extract<double> node_Z( nodeCoordinates[3*nodesIndex+2] );
            
            boost::python::extract<unsigned int> nodeID( nodeIDs[nodesIndex] );
            
            boost::python::extract<double> vel_X( velocity[3*nodesIndex] );
            boost::python::extract<double> vel_Y( velocity[3*nodesIndex+1] );
            boost::python::extract<double> vel_Z( velocity[3*nodesIndex+2] );
            
            boost::python::extract<double> disp_X( displacement[3*nodesIndex] );
            boost::python::extract<double> disp_Y( displacement[3*nodesIndex+1] );
            boost::python::extract<double> disp_Z( displacement[3*nodesIndex+2] ); 
            
            // ######## ADDING NEW NODE #########
            Node < 3 >::Pointer pnode = mr_interface_part.CreateNewNode(nodeID,node_X,node_Y,node_Z);
            array_1d<double,3> vel;
            vel[0] = vel_X;
            vel[1] = vel_Y;
            vel[2] = vel_Z;
            pnode->GetSolutionStepValue(VELOCITY) = vel;

            array_1d<double,3> disp;
            disp[0] = disp_X;
            disp[1] = disp_Y;
            disp[2] = disp_Z;
            pnode->GetSolutionStepValue(DISPLACEMENT) = disp;
        }

        // Adding new Elements

        // required information for a new element: Id, Geometry
        const unsigned int size_elements = boost::python::len(connectivity) / 3;
        mr_interface_part.Elements().reserve(size_elements);

        // Push back elements
        for (unsigned int elemIndex=0; elemIndex!=size_elements; ++elemIndex)
        {
            boost::python::extract<unsigned int> node1_ID( connectivity[3*elemIndex] );
            boost::python::extract<unsigned int> node2_ID( connectivity[3*elemIndex+1] );
            boost::python::extract<unsigned int> node3_ID( connectivity[3*elemIndex+2] );

            // Create new geometry
            Geometry<Node < 3 > >::Pointer pInterfaceGeom;

            Node < 3 >& point1 = mr_interface_part.Nodes()[node1_ID];
            Node < 3 >& point2 = mr_interface_part.Nodes()[node2_ID];
            Node < 3 >& point3 = mr_interface_part.Nodes()[node3_ID];

            //Node < 3 > ::Pointer pnode1 = Node < 3 > ::Pointer(new Node<3>(node1_ID, point1[0], point1[1], point1[2] ));
            //Node < 3 > ::Pointer pnode2 = Node < 3 > ::Pointer(new Node<3>(node2_ID, point2[0], point2[1], point2[2] ));
            //Node < 3 > ::Pointer pnode3 = Node < 3 > ::Pointer(new Node<3>(node3_ID, point3[0], point3[1], point3[2] ));

            pInterfaceGeom = Geometry<Node < 3 > >::Pointer(new Triangle3D3< Node<3> >(point1,point2,point3));

            // Create and insert new element
            unsigned int elemID = elemIndex + 1;
            Element::Pointer pElement = Element::Pointer(new MembraneElement(
                                                            elemID,
                                                            pInterfaceGeom ) );

            mr_interface_part.Elements().push_back(pElement);
        }

        KRATOS_CATCH("")
    }

    // ##############################################################################

    void ExtractMeshInfo( boost::python::list& numNodes, boost::python::list& numElems,
                          boost::python::list& nodes, boost::python::list& nodeIDs,
                          boost::python::list& numNodesPerElem, boost::python::list& elems)
    {
        KRATOS_TRY

        unsigned int nodesCounter = 0;
        unsigned int elemsCounter = 0;

        // loop over all fluid nodes
        for ( ModelPart::NodesContainerType::iterator i_Node =  mr_interface_part.NodesBegin();
                                                      i_Node != mr_interface_part.NodesEnd();
                                                      ++i_Node )
        {
            double node_X = i_Node->Coordinates()[0];
            double node_Y = i_Node->Coordinates()[1];
            double node_Z = i_Node->Coordinates()[2];

            // Fill the nodes vector with nodal coordinates
            nodes.append(node_X);
            nodes.append(node_Y);
            nodes.append(node_Z);

            nodesCounter++;

            // Fill the nodeIDs vector with the nodal IDs
            nodeIDs.append(i_Node->Id());
        }

        // loop over all fluid elements
        for( ModelPart::ConditionsContainerType::iterator i_Condition =  mr_interface_part.ConditionsBegin();
                                                          i_Condition != mr_interface_part.ConditionsEnd();
                                                          i_Condition++)
        {
            unsigned int nodesPerElem = i_Condition->GetGeometry().size();
            numNodesPerElem.append(nodesPerElem);

            if( nodesPerElem == 2 ) // interface conditions are lines (2D fluid mesh)
            {
                unsigned int nodeID_1 = i_Condition->GetGeometry()[0].Id();
                unsigned int nodeID_2 = i_Condition->GetGeometry()[1].Id();

                elems.append(nodeID_1);
                elems.append(nodeID_2);

                elemsCounter++;
            }
            else if( nodesPerElem == 3 ) // interface conditions are triangles (3D fluid mesh)
            {
                unsigned int nodeID_1 = i_Condition->GetGeometry()[0].Id();
                unsigned int nodeID_2 = i_Condition->GetGeometry()[1].Id();
                unsigned int nodeID_3 = i_Condition->GetGeometry()[2].Id();

                elems.append(nodeID_1);
                elems.append(nodeID_2);
                elems.append(nodeID_3);

                elemsCounter++;
            }
	    else if ( nodesPerElem == 4 ) // interface conditions are quads (2D adaptor with z-offset)
	      {
		unsigned int nodeID_1 = i_Condition->GetGeometry()[0].Id();
                unsigned int nodeID_2 = i_Condition->GetGeometry()[1].Id();
                unsigned int nodeID_3 = i_Condition->GetGeometry()[2].Id();
                unsigned int nodeID_4 = i_Condition->GetGeometry()[3].Id();

                elems.append(nodeID_1);
                elems.append(nodeID_2);
                elems.append(nodeID_3);
                elems.append(nodeID_4);

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
    NodesContainerMapType mDummyNodes2D;
    const int m_dimension;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{
    template< class TDataType >
    void SynchronizeDummy2DVariable(Variable<TDataType>& ThisVariable)
    {
      KRATOS_TRY
	for (NodesContainerMapType::iterator itDummy = mDummyNodes2D.begin(); itDummy != mDummyNodes2D.end(); itDummy++)
	  {
	    NodesContainerType::iterator itNode = mr_interface_part.Nodes().find(itDummy.key());
	    if (itNode == mr_interface_part.Nodes().end())
	      KRATOS_THROW_ERROR(std::logic_error, "Cannot find parent node of dummy node","");
	    itDummy.base()->second->FastGetSolutionStepValue(ThisVariable) = itNode->FastGetSolutionStepValue(ThisVariable);
	  }
      KRATOS_CATCH("")
    }


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
