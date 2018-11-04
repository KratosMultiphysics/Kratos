//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_CYLINDER_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_CYLINDER_BOUNDING_BOX_H_INCLUDED

// External includes

// System includes

// Project includes
#include "custom_bounding/spatial_bounding_box.hpp"
#include "geometries/quadrilateral_3d_4.h"

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

    This Box represents a 3D wall composed by a cylinder

    A convexity parameter is given to determine if
    the internal or external space is considered as boundary

    This bounding box is essentially used for rigid wall contact purposes
*/

class CylinderBoundingBox
  : public SpatialBoundingBox
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CylinderBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( CylinderBoundingBox );

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
    CylinderBoundingBox() : SpatialBoundingBox()
    {
      KRATOS_TRY

      std::cout<< "Calling Cylinder Wall BBX empty constructor" <<std::endl;

      KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************

    CylinderBoundingBox(Parameters CustomParameters)
    {

      KRATOS_TRY

      Parameters DefaultParameters( R"(
            {
                "parameters_list":[{
                   "first_center":  [0.0, 0.0, 0.0],
                   "second_center": [0.0, 0.0, 0.0],
                   "radius": 0.0,
                   "convexity": 1
                 }],
                 "velocity": [0.0, 0.0, 0.0]

            }  )" );


      //validate against defaults -- this also ensures no type mismatch
      CustomParameters.ValidateAndAssignDefaults(DefaultParameters);

      if(CustomParameters["parameters_list"].IsArray() == true && CustomParameters["parameters_list"].size() != 1)
        {
	  KRATOS_THROW_ERROR(std::runtime_error,"paramters_list for the Cylinder BBX must contain only one term",CustomParameters.PrettyPrintJsonString());
        }

      mBox.Initialize();

      Parameters BoxParameters = CustomParameters["parameters_list"][0];

      mFirstCenter[0] = BoxParameters["first_center"][0].GetDouble();
      mFirstCenter[1] = BoxParameters["first_center"][1].GetDouble();
      mFirstCenter[2] = BoxParameters["first_center"][2].GetDouble();

      mSecondCenter[0] = BoxParameters["second_center"][0].GetDouble();
      mSecondCenter[1] = BoxParameters["second_center"][1].GetDouble();
      mSecondCenter[2] = BoxParameters["second_center"][2].GetDouble();

      mBox.Center    = 0.5 * (mFirstCenter+mSecondCenter);

      mBox.Radius = BoxParameters["radius"].GetDouble();

      mBox.Velocity[0] = CustomParameters["velocity"][0].GetDouble();
      mBox.Velocity[1] = CustomParameters["velocity"][1].GetDouble();
      mBox.Velocity[2] = CustomParameters["velocity"][2].GetDouble();

      mBox.Convexity = BoxParameters["convexity"].GetInt();

      PointType Side;
      Side[0] = mBox.Radius;
      Side[1] = mBox.Radius;
      Side[2] = mBox.Radius;

      mBox.UpperPoint = mBox.Center + Side;
      mBox.LowerPoint = mSecondCenter - Side;

      mBox.SetInitialValues();

      mRigidBodyCenterSupplied = false;

      KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************


    // General Wall constructor
    CylinderBoundingBox(PointType FirstCenter,
			PointType SecondCenter,
			double Radius,
			PointType Velocity,
			int Convexity)

    {
      KRATOS_TRY

      std::cout<<" [--CYLINDER WALL--] "<<std::endl;

      mFirstCenter   = FirstCenter;
      mSecondCenter  = SecondCenter;

      mBox.Center    = 0.5 * (FirstCenter+SecondCenter);
      mBox.Radius    = Radius;
      mBox.Convexity = Convexity;

      std::cout<<"  [Convexity:"<<mBox.Convexity<<std::endl;
      std::cout<<"  [Radius:"<<mBox.Radius<<std::endl;
      std::cout<<"  [Center1:"<<mFirstCenter<<std::endl;
      std::cout<<"  [Center2:"<<mSecondCenter<<std::endl;
      std::cout<<" [--------] "<<std::endl;

      mBox.Velocity = Velocity;

      mBox.SetInitialValues();

      mRigidBodyCenterSupplied = false;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    /// Assignment operator.
    CylinderBoundingBox& operator=(CylinderBoundingBox const& rOther)
    {
      KRATOS_TRY

      SpatialBoundingBox::operator=(rOther);

      mFirstCenter = rOther.mFirstCenter;
      mSecondCenter = rOther.mSecondCenter;

      return *this;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    /// Copy constructor.
    CylinderBoundingBox(CylinderBoundingBox const& rOther)
      :SpatialBoundingBox(rOther)
      ,mFirstCenter(rOther.mFirstCenter)
      ,mSecondCenter(rOther.mSecondCenter)
    {
    }


    //**************************************************************************
    //**************************************************************************

    /// Destructor.
    virtual ~CylinderBoundingBox() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    //************************************************************************************
    //************************************************************************************

    bool IsInside (const PointType& rPoint, double& rCurrentTime, double Radius = 0) override
    {

      KRATOS_TRY

      bool is_inside = false;

      double CylinderRadius = mBox.Radius;

      //outside
      if( mBox.Convexity == 1)
	CylinderRadius *= 1.25; //increase the bounding box

      //inside
      if( mBox.Convexity == -1)
       	CylinderRadius *= 0.75; //decrease the bounding box

      is_inside = ContactSearch(rPoint, CylinderRadius);

      return is_inside;

      KRATOS_CATCH("")
    }


    //************************************************************************************
    //************************************************************************************

    bool IsInside(BoundingBoxParameters& rValues, const ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY

      bool is_inside = false;

      rValues.SetContactFace(2);

      is_inside = ContactSearch(rValues, rCurrentProcessInfo);

      return is_inside;

      KRATOS_CATCH("")
    }



    //************************************************************************************
    //************************************************************************************

    //Cylinder
    void CreateBoundingBoxBoundaryMesh(ModelPart& rModelPart, int linear_partitions = 4, int angular_partitions = 4 ) override
    {
      KRATOS_TRY

      //std::cout<<" Create Cylinder Mesh "<<std::endl;

      //1.-create generatrix
      PointType Axis = (mSecondCenter - mFirstCenter);
      double AxisLength = norm_2(Axis);

      if( AxisLength )
	Axis/=AxisLength;


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
	    for(ModelPart::NodesContainerType::iterator i_node = i_mp->NodesBegin() ; i_node != i_mp->NodesEnd() ; ++i_node)
	      {
		if( i_node->Id() == rModelPart.Nodes().front().Id() ){
		  BoundaryModelPartsName.push_back(i_mp->Name());
		  break;
		}
	      }
	  }
	}
      //get boundary model parts ( temporary implementation )


      double SingleLength = AxisLength / (double)linear_partitions;

      PointType Point(3);

      // Create Axis generatrix
      // for(int i=0; i<linear_partitions; i++)
      // 	{
      // 	  Point =  mFirstCenter + Axis * ( SingleLength );
      // 	  NodeId += 1;
      // 	  NodeType::Pointer pNode = this->CreateNode(rModelPart, BasePoint, NodeId);

      // 	  pNode->Set(RIGID,true);
      // 	  rModelPart.AddNode( pNode );
      // 	}


      PointType DirectionX(3);
      noalias(DirectionX) = Axis;
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

      for(int i=0; i<linear_partitions; i++)
	{
	  Point =  mFirstCenter + Axis * ( SingleLength ) * i /(double)linear_partitions;

	  for(int k=0; k<angular_partitions; k++)
	    {
	      alpha = (2.0 * 3.14159262 * k)/(double)angular_partitions;

	      //vector of rotation
	      RotationAxis = DirectionX * alpha;
	      Quaternion   = QuaternionType::FromRotationVector(RotationAxis);

	      RotatedDirectionY = DirectionY;

	      Quaternion.RotateVector3(RotatedDirectionY);

	      // std::cout<<" alpha "<<alpha<<"  cos "<<Q(1,1)<<std::endl;
	      //std::cout<<" Rotated "<<RotatedDirectionY<<" alpha "<<alpha<<std::endl;

	      //add the angular_partitions points number along the circular base of the cylinder
	      NodeId += 1;
	      noalias(BasePoint) = Point + mBox.Radius * RotatedDirectionY;

	      //std::cout<<" BasePoint["<<NodeId<<"] "<<BasePoint<<" radius "<<mBox.Radius<<std::endl;

	      NodeType::Pointer pNode = this->CreateNode(rModelPart, BasePoint, NodeId);

	      pNode->Set(RIGID,true);
	      rModelPart.AddNode( pNode );

	      //get boundary model parts ( temporary implementation )
	      for(unsigned int j=0; j<BoundaryModelPartsName.size(); j++)
		(pMainModelPart->GetSubModelPart(BoundaryModelPartsName[j])).AddNode( pNode );
	      //get boundary model parts ( temporary implementation )

	    }
	}

      //std::cout<<" Create Cylinder Mesh NODES "<<std::endl;

      this->CreateQuadrilateralBoundaryMesh(rModelPart, InitialNodeId, angular_partitions);

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
        return "CylinderBoundingBox";
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

    PointType mFirstCenter;
    PointType mSecondCenter;

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


    bool ContactSearch(const PointType& rPoint, const double& rRadius)
    {

      KRATOS_TRY

      //1.-compute point projection
      PointType Axis = (mSecondCenter - mFirstCenter);
      if( norm_2(Axis) )
	Axis/=norm_2(Axis);

      PointType Projection(3);
      Projection = inner_prod( (rPoint-mFirstCenter), Axis) * Axis;

      //2.-compute gap
      double GapNormal = norm_2(rPoint-Projection)-rRadius;

      GapNormal *= mBox.Convexity;

      if(GapNormal<0)
	return true;
      else
	return false;


      KRATOS_CATCH( "" )

   }


    //************************************************************************************
    //************************************************************************************

    //Cylinder (note: box position has been updated previously)
    bool ContactSearch(BoundingBoxParameters& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      const PointType& rPoint = rValues.GetPoint();
      PointType& rNormal      = rValues.GetNormal();
      PointType& rTangent     = rValues.GetTangent();

      double& rGapNormal      = rValues.GetGapNormal();

      rNormal  = ZeroVector(3);
      rTangent = ZeroVector(3);

      //1.-compute point projection
      PointType Axis = (mSecondCenter - mFirstCenter);
      if( norm_2(Axis) )
	Axis/=norm_2(Axis);

      PointType Projection(3);
      Projection = inner_prod( (rPoint-mFirstCenter), Axis) * Axis;


      //2.-compute contact normal
      rNormal = rPoint-Projection;

      if(norm_2(rNormal))
	rNormal /= norm_2(rNormal);

      rNormal *= mBox.Convexity;

      //3.-compute gap
      rGapNormal = norm_2(rPoint-Projection)-mBox.Radius;

      rGapNormal *= mBox.Convexity;

      this->ComputeContactTangent(rValues,rCurrentProcessInfo);

      if(rGapNormal<0)
	return true;
      else
	return false;

      KRATOS_CATCH( "" )
    }


    //************************************************************************************
    //************************************************************************************

    void CreateQuadrilateralBoundaryMesh(ModelPart& rModelPart, const unsigned int& rInitialNodeId, const unsigned int& angular_partitions )
    {

      KRATOS_TRY

      //std::cout<<" Create Cylinder Mesh ELEMENTS "<<std::endl;

      //add elements to computing model part: (in order to be written)
      ModelPart* pComputingModelPart = NULL;
      if( rModelPart.IsSubModelPart() )
	for(ModelPart::SubModelPartIterator i_mp= rModelPart.GetParentModelPart()->SubModelPartsBegin() ; i_mp!=rModelPart.GetParentModelPart()->SubModelPartsEnd(); i_mp++)
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
      GeometryType::Pointer  pFace;
      ElementType::Pointer             pElement;

      //PROPERTIES:
      int number_of_properties = rModelPart.NumberOfProperties();
      Properties::Pointer pProperties = Kratos::make_shared<Properties>(number_of_properties);

      unsigned int local_counter = 1;
      unsigned int counter       = 0;
      unsigned int Id   = rInitialNodeId;

      Vector FaceNodesIds = ZeroVector(4);

      while(Id < NodeId){

	counter += 1;
	ElementId += 1;

	if( local_counter < angular_partitions ){

	  FaceNodesIds[0] = rInitialNodeId + counter ;
	  FaceNodesIds[1] = rInitialNodeId + counter + 1;
	  FaceNodesIds[2] = rInitialNodeId + counter + angular_partitions + 1;
	  FaceNodesIds[3] = rInitialNodeId + counter + angular_partitions;

	  //std::cout<<" ElementA "<<FaceNodesIds<<std::endl;

	  GeometryType::PointsArrayType FaceNodes;
	  FaceNodes.reserve(4);

	  //NOTE: when creating a PointsArrayType
	  //important ask for pGetNode, if you ask for GetNode a copy is created
	  //if a copy is created a segmentation fault occurs when the node destructor is called

	  for(unsigned int j=0; j<4; j++)
	    FaceNodes.push_back(rModelPart.pGetNode(FaceNodesIds[j]));

	  pFace    = Kratos::make_shared<Quadrilateral3D4<NodeType> >(FaceNodes);
	  pElement = Kratos::make_shared<Element>(ElementId, pFace, pProperties);

	  rModelPart.AddElement(pElement);
	  pElement->Set(ACTIVE,false);
	  pComputingModelPart->AddElement(pElement);

	  local_counter++;

	}
	else if( local_counter == angular_partitions ){

	  FaceNodesIds[0] = rInitialNodeId + counter ;
	  FaceNodesIds[1] = rInitialNodeId + counter + 1 - angular_partitions;
	  FaceNodesIds[2] = rInitialNodeId + counter + 1;
	  FaceNodesIds[3] = rInitialNodeId + counter + angular_partitions;

	  //std::cout<<" ElementB "<<FaceNodesIds<<std::endl;

	  GeometryType::PointsArrayType FaceNodes;
	  FaceNodes.reserve(4);

	  for(unsigned int j=0; j<4; j++)
	    FaceNodes.push_back(rModelPart.pGetNode(FaceNodesIds[j]));

	  pFace    = Kratos::make_shared<Quadrilateral3D4<NodeType> >(FaceNodes);
	  pElement = Kratos::make_shared<Element>(ElementId, pFace, pProperties);

	  rModelPart.AddElement(pElement);
	  pElement->Set(ACTIVE,false);
	  pComputingModelPart->AddElement(pElement);

	  local_counter = 1;
	}

	Id = rInitialNodeId + counter + angular_partitions;
      }


      //std::cout<<" Create Cylinder Mesh ELEMENTS CREATED "<<std::endl;

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


}; // Class CylinderBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_CYLINDER_BOUNDING_BOX_H_INCLUDED  defined
