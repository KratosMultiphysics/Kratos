//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_CIRCLE_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_CIRCLE_BOUNDING_BOX_H_INCLUDED

// External includes

// System includes

// Project includes
#include "custom_bounding/sphere_bounding_box.hpp"
#include "geometries/line_2d_2.h"

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

    This Box represents a 2D wall composed by a circle

    A convexity parameter is given to determine which side of each nose is considered
    the internal or external boundary

    This bounding box is essentially used for rigid wall contact purposes
*/

class CircleBoundingBox
  : public SphereBoundingBox
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CircleBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( CircleBoundingBox );

    //typedef BoundedVector<double, 3>                     PointType;
    typedef array_1d<double, 3>                             PointType;
    typedef ModelPart::NodeType                              NodeType;
    typedef ModelPart::NodesContainerType          NodesContainerType;
    typedef NodesContainerType::Pointer     NodesContainerTypePointer;
    typedef Quaternion<double>                         QuaternionType;
    typedef ModelPart::ElementType                        ElementType;
    typedef ElementType::GeometryType                    GeometryType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CircleBoundingBox() : SphereBoundingBox()
    {
      KRATOS_TRY

      std::cout<< "Calling Rigid Circle Wall BBX empty constructor" <<std::endl;

      KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************

    CircleBoundingBox(Parameters CustomParameters)
      : SphereBoundingBox(CustomParameters)
    {
      KRATOS_TRY

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************


    // General Wall constructor
    CircleBoundingBox(PointType Center,
		      double Radius,
		      PointType Velocity,
		      int Convexity)
      : SphereBoundingBox(Center, Radius, Velocity, Convexity)
    {
      KRATOS_TRY

      KRATOS_CATCH("")
    }



    //**************************************************************************
    //**************************************************************************

    /// Assignment operator.
    CircleBoundingBox& operator=(CircleBoundingBox const& rOther)
    {
      KRATOS_TRY

      SphereBoundingBox::operator=(rOther);
      return *this;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    /// Copy constructor.
    CircleBoundingBox(CircleBoundingBox const& rOther)
      :SphereBoundingBox(rOther)
    {
    }


    //**************************************************************************
    //**************************************************************************

    /// Destructor.
    virtual ~CircleBoundingBox() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //************************************************************************************
    //************************************************************************************

    //Circle
    void CreateBoundingBoxBoundaryMesh(ModelPart& rModelPart, int linear_partitions = 4, int angular_partitions = 4 ) override
    {
      KRATOS_TRY

      unsigned int NodeId = 0;
      if( rModelPart.IsSubModelPart() )
	NodeId = this->GetMaxNodeId( *(rModelPart.GetParentModelPart()) );
      else
	NodeId = this->GetMaxNodeId( rModelPart );

      unsigned int InitialNodeId = NodeId;

      //get boundary model parts ( temporary implementation )
      std::vector<std::string> BoundaryModelPartsName;

      ModelPart* pMainModelPart = &rModelPart;
      if( rModelPart.IsSubModelPart() )
	pMainModelPart = rModelPart.GetParentModelPart();

      for(ModelPart::SubModelPartIterator i_mp= pMainModelPart->SubModelPartsBegin() ; i_mp!=pMainModelPart->SubModelPartsEnd(); i_mp++)
	{
	  if( i_mp->Is(BOUNDARY) || i_mp->Is(ACTIVE) ){
	    for(ModelPart::NodesContainerType::iterator i_node = i_mp->NodesBegin() ; i_node != i_mp->NodesEnd() ; i_node++)
	      {
		if( i_node->Id() == rModelPart.Nodes().front().Id() ){
		  BoundaryModelPartsName.push_back(i_mp->Name());
		  break;
		}
	      }
	  }
	}
      //get boundary model parts ( temporary implementation )

      PointType DirectionX(3);
      DirectionX[0] = 0;
      DirectionX[1] = 0;
      DirectionX[2] = 1;

      PointType DirectionY(3);
      noalias(DirectionY) = ZeroVector(3);
      PointType DirectionZ(3);
      noalias(DirectionZ) = ZeroVector(3);

      this->CalculateOrthonormalBase(DirectionX, DirectionY, DirectionZ);

      PointType BasePoint(3);
      PointType RotationAxis(3);
      PointType RotatedDirectionY(3);

      double alpha = 0;
      QuaternionType Quaternion;

      for(int k=0; k<angular_partitions; k++)
	{
	  alpha = (2.0 * 3.14159262 * k)/(double)angular_partitions;

	  //vector of rotation
	  RotationAxis = DirectionX * alpha;
	  Quaternion   = QuaternionType::FromRotationVector(RotationAxis);

	  RotatedDirectionY = DirectionY;

	  Quaternion.RotateVector3(RotatedDirectionY);

	  //std::cout<<" Rotated "<<RotatedDirectionY<<" alpha "<<alpha<<std::endl;

	  //add the angular_partitions points number along the circle
	  NodeId += 1;

	  std::cout<<" node id "<<NodeId<<std::endl;

	  noalias(BasePoint) = mBox.Center + mBox.Radius * RotatedDirectionY;

	  NodeType::Pointer pNode = this->CreateNode(rModelPart, BasePoint, NodeId);

	  pNode->Set(RIGID,true);

	  rModelPart.AddNode( pNode );

	  //get boundary model parts ( temporary implementation )
	  for(unsigned int j=0; j<BoundaryModelPartsName.size(); j++)
	    (pMainModelPart->GetSubModelPart(BoundaryModelPartsName[j])).AddNode( pNode );
	  //get boundary model parts ( temporary implementation )

	}

      //std::cout<<" Nodes Added "<<NodeId-InitialNodeId<<std::endl;

      this->CreateLinearBoundaryMesh(rModelPart, InitialNodeId);

      KRATOS_CATCH( "" )
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
        return "CircleBoundingBox";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << this->mBox.UpperPoint << " , " << this->mBox.LowerPoint;
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

    //************************************************************************************
    //************************************************************************************

    void CreateLinearBoundaryMesh(ModelPart& rModelPart, const unsigned int& rInitialNodeId)
    {

      KRATOS_TRY

      //add elements to computing model part: (in order to be written)
      ModelPart* pComputingModelPart = NULL;
      if( rModelPart.IsSubModelPart() )
	for(ModelPart::SubModelPartIterator i_mp= rModelPart.GetParentModelPart()->SubModelPartsBegin() ; i_mp!=rModelPart.GetParentModelPart()->SubModelPartsEnd(); ++i_mp)
	  {
	    if( i_mp->Is(ACTIVE) )  //computing_domain
	      pComputingModelPart = &rModelPart.GetParentModelPart()->GetSubModelPart(i_mp->Name());
	  }
      else{
	for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); i_mp++)
	  {
	    if( i_mp->Is(ACTIVE) )  //computing_domain
	      pComputingModelPart = &rModelPart.GetSubModelPart(i_mp->Name());
	  }
      }

      // Create surface of the cylinder/tube with quadrilateral shell conditions
      unsigned int ElementId = 0;
      if( rModelPart.IsSubModelPart() )
	ElementId = this->GetMaxElementId( *(rModelPart.GetParentModelPart()) );
      else
	ElementId = this->GetMaxElementId( rModelPart );

      unsigned int NodeId = 0;
      if( rModelPart.IsSubModelPart() )
	NodeId = this->GetMaxNodeId( *(rModelPart.GetParentModelPart()) );
      else
	NodeId = this->GetMaxNodeId( rModelPart );


      //GEOMETRY:
      GeometryType::Pointer pFace;
      ElementType::Pointer pElement;

      //PROPERTIES:
      int number_of_properties = rModelPart.NumberOfProperties();
      Properties::Pointer pProperties = Kratos::make_shared<Properties>(number_of_properties);

      int counter       = 0;
      unsigned int Id   = rInitialNodeId;

      Vector FaceNodesIds = ZeroVector(2);

      while(Id < NodeId){

	counter += 1;
	ElementId += 1;

	FaceNodesIds[0] = rInitialNodeId + counter ;
	FaceNodesIds[1] = rInitialNodeId + counter + 1;

	//std::cout<<" FaceNodesIds "<<FaceNodesIds<<" element id "<<ElementId<<std::endl;

	GeometryType::PointsArrayType FaceNodes;
	FaceNodes.reserve(2);

	//NOTE: when creating a PointsArrayType
	//important ask for pGetNode, if you ask for GetNode a copy is created
	//if a copy is created a segmentation fault occurs when the node destructor is called

	for(unsigned int j=0; j<2; j++)
	  FaceNodes.push_back(rModelPart.pGetNode(FaceNodesIds[j]));

	pFace    = Kratos::make_shared<Line2D2<NodeType> >(FaceNodes);
        pElement = Kratos::make_shared<Element>(ElementId, pFace, pProperties);

	rModelPart.AddElement(pElement);
	pElement->Set(ACTIVE,false);
	pComputingModelPart->AddElement(pElement);

	Id = rInitialNodeId + counter + 1;

      }


      ElementId += 1;

      //last element
      FaceNodesIds[0] = rInitialNodeId + counter + 1;
      FaceNodesIds[1] = rInitialNodeId + 1;

      GeometryType::PointsArrayType FaceNodes;
      FaceNodes.reserve(2);

      //NOTE: when creating a PointsArrayType
      //important ask for pGetNode, if you ask for GetNode a copy is created
      //if a copy is created a segmentation fault occurs when the node destructor is called

      for(unsigned int j=0; j<2; j++)
	FaceNodes.push_back(rModelPart.pGetNode(FaceNodesIds[j]));

      //std::cout<<" FaceNodesIds "<<FaceNodesIds<<" element id "<<ElementId<<std::endl;

      pFace    = Kratos::make_shared<Line2D2<NodeType> >(FaceNodes);
      pElement = Kratos::make_shared<Element>(ElementId, pFace, pProperties);

      rModelPart.AddElement(pElement);
      pElement->Set(ACTIVE,false);
      pComputingModelPart->AddElement(pElement);


      KRATOS_CATCH( "" )

    }

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


    ///@}


}; // Class CircleBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_CIRCLE_BOUNDING_BOX_H_INCLUDED  defined
