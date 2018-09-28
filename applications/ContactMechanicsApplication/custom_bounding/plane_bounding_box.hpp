//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_PLANE_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_PLANE_BOUNDING_BOX_H_INCLUDED

// External includes

// System includes

// Project includes
#include "custom_bounding/spatial_bounding_box.hpp"
#include "geometries/line_2d_2.h"
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

    This Box represents a 2D wall composed by plane

    A convexity parameter is given to determine which side of each nose is considered
    the internal or external boundary

    This bounding box is essentially used for rigid wall contact purposes
*/

class PlaneBoundingBox
  : public SpatialBoundingBox
{
public:

    //typedef BoundedVector<double, 3>                     PointType;
    typedef array_1d<double, 3>                             PointType;
    typedef ModelPart::NodeType                              NodeType;
    typedef ModelPart::NodesContainerType          NodesContainerType;
    typedef NodesContainerType::Pointer     NodesContainerTypePointer;
    typedef Quaternion<double>                         QuaternionType;
    typedef ModelPart::ElementType                        ElementType;
    typedef ElementType::GeometryType                    GeometryType;

protected:

    typedef struct
    {

      PointType  Point;      // plane point
      PointType  Normal;     // plane normal


    public:

      void Initialize()
      {

	Point.clear();
	Normal.clear();

      }

    } PlaneVariables;


public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PlaneBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( PlaneBoundingBox );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PlaneBoundingBox() : SpatialBoundingBox()
    {
      KRATOS_TRY

      std::cout<< "Calling Rigid Plane Wall BBX empty constructor" <<std::endl;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    PlaneBoundingBox(Parameters CustomParameters)
    {

      KRATOS_TRY

      Parameters DefaultParameters( R"(
            {
                "parameters_list":[{
                   "point": [0.0, 0.0, 0.0],
                   "normal": [0.0, 0.0, 0.0],
                   "convexity": 1
                 }],
                 "velocity": [0.0, 0.0, 0.0],
                 "plane_size": 1.0

            }  )" );


      //validate against defaults -- this also ensures no type mismatch
      CustomParameters.ValidateAndAssignDefaults(DefaultParameters);

      if(CustomParameters["parameters_list"].IsArray() == true && CustomParameters["parameters_list"].size() != 1)
        {
	  KRATOS_THROW_ERROR(std::runtime_error,"paramters_list for the Plane BBX must contain only one term",CustomParameters.PrettyPrintJsonString());
        }

      mBox.Initialize();
      mPlane.Initialize();

      Parameters BoxParameters = CustomParameters["parameters_list"][0];

      mPlane.Point[0] = BoxParameters["point"][0].GetDouble();
      mPlane.Point[1] = BoxParameters["point"][1].GetDouble();
      mPlane.Point[2] = BoxParameters["point"][2].GetDouble();

      mPlane.Normal[0] = BoxParameters["normal"][0].GetDouble();
      mPlane.Normal[1] = BoxParameters["normal"][1].GetDouble();
      mPlane.Normal[2] = BoxParameters["normal"][2].GetDouble();

      mBox.Center = mPlane.Point;

      mBox.Velocity[0] = CustomParameters["velocity"][0].GetDouble();
      mBox.Velocity[1] = CustomParameters["velocity"][1].GetDouble();
      mBox.Velocity[2] = CustomParameters["velocity"][2].GetDouble();

      mBox.Radius = CustomParameters["plane_size"].GetDouble();

      PointType UpperPoint(3);
      PointType LowerPoint(3);

      //calculate upper and lower points
      for(unsigned int i=0; i<3; i++)
	{
	  mBox.UpperPoint[i] = mBox.Radius;
	  mBox.LowerPoint[i] = (-1) * mBox.Radius;
	}

      mBox.UpperPoint += mPlane.Point;
      mBox.LowerPoint += mPlane.Point;

      mBox.Convexity = BoxParameters["convexity"].GetInt();

      mBox.SetInitialValues();

      mRigidBodyCenterSupplied = false;

      KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************

    // General Wall constructor
    PlaneBoundingBox(PointType Point,
		     PointType Normal,
		     PointType Velocity,
		     int Convexity)
    {
      KRATOS_TRY

      std::cout<<" [--PLANE WALL--] "<<std::endl;

      mBox.Center  = Point;
      mPlane.Point = Point;
      mBox.Convexity = Convexity;

      if( norm_2(Normal) )
	mPlane.Normal = Normal/norm_2(Normal);
      else
	std::cout<<" [ERROR: Normal is Zero]"<<std::endl;


      std::cout<<"  [Convexity:"<<mBox.Convexity<<std::endl;
      std::cout<<"  [Point:"<<mPlane.Point<<std::endl;
      std::cout<<"  [Normal:"<<mPlane.Normal<<std::endl;
      std::cout<<" [--------] "<<std::endl;

      mBox.Velocity = Velocity;

      mBox.SetInitialValues();

      mRigidBodyCenterSupplied = false;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************


    /// Assignment operator.
    PlaneBoundingBox& operator=(PlaneBoundingBox const& rOther)
    {
      KRATOS_TRY

      SpatialBoundingBox::operator=(rOther);
      mPlane = rOther.mPlane;
      return *this;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    /// Copy constructor.
    PlaneBoundingBox(PlaneBoundingBox const& rOther)
      :SpatialBoundingBox(rOther)
      ,mPlane(rOther.mPlane)
    {
    }

    //**************************************************************************
    //**************************************************************************

    /// Destructor.
    virtual ~PlaneBoundingBox() {};


    ///@}
    ///@name Operators
    ///@{



    ///@}
    ///@name Operations
    ///@{

    //**************************************************************************
    //**************************************************************************

    void UpdateBoxPosition(const double & rCurrentTime) override
    {

      KRATOS_TRY

      PointType Displacement  =  this->GetBoxDisplacement(rCurrentTime);

      mBox.UpdatePosition(Displacement);

      mPlane.Point = mBox.Center;

      KRATOS_CATCH("")

    }


    //************************************************************************************
    //************************************************************************************


    bool IsInside (const PointType& rPoint, double& rCurrentTime, double Radius = 0) override
    {

      KRATOS_TRY

      bool is_inside = false;

      PointType  Displacement = this->GetBoxDisplacement(rCurrentTime);

      PointType  PlanePoint   = mBox.InitialCenter + Displacement;

      if( mBox.Convexity == 1 )
	PlanePoint += mPlane.Normal * 0.1; //increase the bounding box

      if( mBox.Convexity == -1 )
       	PlanePoint -= mPlane.Normal * 0.1; //decrease the bounding box

      is_inside = ContactSearch(rPoint, PlanePoint);


      return is_inside;

      KRATOS_CATCH("")
    }


    //************************************************************************************
    //************************************************************************************

    bool IsInside(BoundingBoxParameters& rValues, const ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY

      bool is_inside = false;

      rValues.SetContactFace(0);
      rValues.SetRadius(0.0);

      is_inside = ContactSearch(rValues, rCurrentProcessInfo);

      return is_inside;

      KRATOS_CATCH("")
    }


            // *********************************************************************************
            // *********************************************************************************
            virtual void GetParametricDirections(BoundingBoxParameters & rValues, Vector & rT1, Vector & rT2) override
            {
               KRATOS_TRY

               // GetTheNormalOfThePlane
                  PointType Normal(3);
                  noalias(Normal) = mPlane.Normal;
                  PointType T1(3); PointType T2(3);
                  noalias(T1) = ZeroVector(3); noalias(T2) = ZeroVector(3);
                  this->CalculateOrthonormalBase(Normal, T1, T2);

                  for (unsigned int i = 0; i < 3; i++)
                  {
                     rT1(i) = T1(i); rT2(i) = T2(i);
               }

               KRATOS_CATCH("")
            }
    //************************************************************************************
    //************************************************************************************

    //Plane
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

      PointType DirectionX(3);
      noalias(DirectionX) = mPlane.Normal;
      PointType DirectionY(3);
      noalias(DirectionY) = ZeroVector(3);
      PointType DirectionZ(3);
      noalias(DirectionZ) = ZeroVector(3);

      this->CalculateOrthonormalBase(DirectionX, DirectionY, DirectionZ);

      PointType BasePoint(3);
      PointType RotationAxis(3);
      PointType RotatedDirectionY(3);

      //calculate center
      PointType Upper = (mBox.UpperPoint - mPlane.Point);
      Upper -= inner_prod(Upper,mPlane.Normal) * mPlane.Normal;

      PointType Lower = (mBox.LowerPoint - mPlane.Point);
      Lower -= inner_prod(Lower,mPlane.Normal) * mPlane.Normal;

      PointType PlaneCenter = mPlane.Point + 0.5 * (Upper + Lower);
      double PlaneRadius    = norm_2(Upper-Lower);
      PlaneRadius *= 0.5;

      double alpha = 0;
      QuaternionType Quaternion;

      if( rModelPart.GetMesh().WorkingSpaceDimension() == 2 || rModelPart.GetProcessInfo()[SPACE_DIMENSION]==2 )
	angular_partitions = 2;
      else
	angular_partitions = 4;


      for(int k=0; k<angular_partitions; k++)
	{
	  alpha = (2.0 * 3.14159262 * k)/(double)angular_partitions;

	  //vector of rotation
	  RotationAxis = DirectionX * alpha;
	  Quaternion   = QuaternionType::FromRotationVector(RotationAxis);

	  RotatedDirectionY = DirectionY;

	  Quaternion.RotateVector3(RotatedDirectionY);

	  //add the angular_partitions points number along the circle
	  NodeId += 1;

	  noalias(BasePoint) = PlaneCenter + PlaneRadius * RotatedDirectionY;

	  NodeType::Pointer pNode = this->CreateNode(rModelPart, BasePoint, NodeId);

	  pNode->Set(RIGID,true);

	  rModelPart.AddNode( pNode );

	  //get boundary model parts ( temporary implementation )
	  for(unsigned int j=0; j<BoundaryModelPartsName.size(); j++)
	    (pMainModelPart->GetSubModelPart(BoundaryModelPartsName[j])).AddNode( pNode );
	  //get boundary model parts ( temporary implementation )

	}

      //std::cout<<" Nodes Added "<<NodeId-InitialNodeId<<std::endl;
      if( rModelPart.GetMesh().WorkingSpaceDimension() == 2 || rModelPart.GetProcessInfo()[SPACE_DIMENSION]==2 ){
	std::cout<<" CREATE a LINE mesh "<<std::endl;
	this->CreateLinearBoundaryMesh(rModelPart, InitialNodeId);
      }
      else{
	std::cout<<" CREATE a QUADRILATERAL mesh "<<std::endl;
	this->CreateQuadrilateralBoundaryMesh(rModelPart, InitialNodeId);
      }

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
    virtual std::string Info() const override
    {
        return "PlaneBoundingBox";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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

    PlaneVariables mPlane;

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
      GeometryType::Pointer pFace;
      ElementType::Pointer pElement;

      //PROPERTIES:
      int number_of_properties = rModelPart.NumberOfProperties();
      Properties::Pointer pProperties = Kratos::make_shared<Properties>(number_of_properties);

      int counter       = 0;
      unsigned int Id   = rInitialNodeId;

      Vector FaceNodesIds(2);
      noalias( FaceNodesIds ) = ZeroVector(2);

      while(Id < NodeId){

	counter   += 1;
	ElementId += 1;

	FaceNodesIds[0] = rInitialNodeId + counter ;
	FaceNodesIds[1] = rInitialNodeId + counter + 1;

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

	Id = rInitialNodeId + counter + 2;

      }

      KRATOS_CATCH( "" )

    }


    //************************************************************************************
    //************************************************************************************

    void CreateQuadrilateralBoundaryMesh(ModelPart& rModelPart, const unsigned int& rInitialNodeId)
    {

      KRATOS_TRY

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
      GeometryType::Pointer pFace;
      ElementType::Pointer pElement;

      //PROPERTIES:
      int number_of_properties = rModelPart.NumberOfProperties();
      Properties::Pointer pProperties = Kratos::make_shared<Properties>(number_of_properties);

      int counter       = 0;
      unsigned int Id   = rInitialNodeId;

      Vector FaceNodesIds(4);
      noalias( FaceNodesIds ) = ZeroVector(4);

      while(Id < NodeId){

	counter += 1;
	ElementId += 1;

	FaceNodesIds[0] = rInitialNodeId + counter ;
	FaceNodesIds[1] = rInitialNodeId + counter + 1;
	FaceNodesIds[2] = rInitialNodeId + counter + 2;
	FaceNodesIds[3] = rInitialNodeId + counter + 3;

	//std::cout<<" FaceNodesIds "<<FaceNodesIds<<" element id "<<ElementId<<std::endl;

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

	Id = rInitialNodeId + counter + 3;

      }

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
    //************************************************************************************
    //************************************************************************************


    bool ContactSearch(const PointType& rPoint, const PointType& rPlanePoint)
    {

      KRATOS_TRY


      //1.-compute gap
      double GapNormal = inner_prod((rPoint - rPlanePoint), mBox.Convexity*mPlane.Normal);


      if(GapNormal<0)
	return true;
      else
	return false;


      KRATOS_CATCH( "" )

   }

    //************************************************************************************
    //************************************************************************************

    //Plane (note: box position has been updated previously)
    bool ContactSearch(BoundingBoxParameters& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      const PointType& rPoint = rValues.GetPoint();
      PointType& rNormal      = rValues.GetNormal();
      PointType& rTangent     = rValues.GetTangent();

      double& rGapNormal      = rValues.GetGapNormal();

      rNormal  = ZeroVector(3);
      rTangent = ZeroVector(3);

      //1.-compute contact normal
      rNormal = mPlane.Normal;

      rNormal *= mBox.Convexity;

      //2.-compute  normal gap
      rGapNormal = inner_prod((rPoint - mPlane.Point), rNormal);

      this->ComputeContactTangent(rValues,rCurrentProcessInfo);

      if(rGapNormal<0)
	return true;
      else
	return false;

      KRATOS_CATCH( "" )
    }


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


}; // Class PlaneBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_PLANE_BOUNDING_BOX_H_INCLUDED  defined
