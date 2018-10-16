//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_SPHERE_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_SPHERE_BOUNDING_BOX_H_INCLUDED

// External includes

// System includes

// Project includes
#include "custom_bounding/spatial_bounding_box.hpp"
#include "geometries/triangle_3d_3.h"
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

    This Box represents a 3D wall composed by a sphere

    A convexity parameter is given to determine if
    the internal or external space is considered as boundary

    This bounding box is essentially used for rigid wall contact purposes
*/

class SphereBoundingBox
  : public SpatialBoundingBox
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SphereBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( SphereBoundingBox );

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
    SphereBoundingBox() : SpatialBoundingBox()
    {
      KRATOS_TRY

      std::cout<< "Calling Rigid Sphere Wall BBX empty constructor" <<std::endl;

      KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************

    SphereBoundingBox(Parameters CustomParameters)
    {

      KRATOS_TRY

      Parameters DefaultParameters( R"(
            {
                "parameters_list":[{
                   "center": [0.0, 0.0, 0.0],
                   "radius": 0.0,
                   "convexity": 1
                 }],
                 "velocity": [0.0, 0.0, 0.0]

            }  )" );


      //validate against defaults -- this also ensures no type mismatch
      CustomParameters.ValidateAndAssignDefaults(DefaultParameters);

      if(CustomParameters["parameters_list"].IsArray() == true && CustomParameters["parameters_list"].size() != 1)
        {
	  KRATOS_THROW_ERROR(std::runtime_error,"paramters_list for the Sphere BBX must contain only one term",CustomParameters.PrettyPrintJsonString());
        }

      mBox.Initialize();

      Parameters BoxParameters = CustomParameters["parameters_list"][0];

      mBox.Center[0] = BoxParameters["center"][0].GetDouble();
      mBox.Center[1] = BoxParameters["center"][1].GetDouble();
      mBox.Center[2] = BoxParameters["center"][2].GetDouble();

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
      mBox.LowerPoint  = mBox.Center - Side;

      mBox.SetInitialValues();

      mRigidBodyCenterSupplied = false;

      KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************


    // General Wall constructor
    SphereBoundingBox(PointType Center,
		      double Radius,
		      PointType Velocity,
		      int Convexity)

    {
      KRATOS_TRY

      std::cout<<" [--CIRCLE/SPHERE WALL--] "<<std::endl;

      mBox.Center = Center;
      mBox.Radius = Radius;
      mBox.Convexity = Convexity;

      std::cout<<"  [Convexity:"<<mBox.Convexity<<std::endl;
      std::cout<<"  [Radius:"<<mBox.Radius<<std::endl;
      std::cout<<"  [Center:"<<mBox.Center<<std::endl;
      std::cout<<" [--------] "<<std::endl;

      mBox.Velocity = Velocity;

      mBox.SetInitialValues();

      mRigidBodyCenterSupplied = false;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    /// Assignment operator.
    SphereBoundingBox& operator=(SphereBoundingBox const& rOther)
    {
      KRATOS_TRY

      SpatialBoundingBox::operator=(rOther);
      return *this;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    /// Copy constructor.
    SphereBoundingBox(SphereBoundingBox const& rOther)
      :SpatialBoundingBox(rOther)
    {
    }


    //**************************************************************************
    //**************************************************************************

    /// Destructor.
    virtual ~SphereBoundingBox() {};


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

      double SphereRadius = mBox.Radius;

      if( mBox.Convexity == 1)
	SphereRadius *= 1.25; //increase the bounding box

      if( mBox.Convexity == -1)
       	SphereRadius *= 0.75; //decrease the bounding box

      is_inside = ContactSearch(rPoint, SphereRadius);

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

            // *********************************************************************************
            // *********************************************************************************
            virtual void GetParametricDirections(BoundingBoxParameters & rValues, Vector & rT1, Vector & rT2) override
            {
               KRATOS_TRY

               // GetTheNormalOfThePlane
               const PointType  & rPoint = rValues.GetPoint();
               PointType Normal(3);

               Normal = rPoint - mBox.Center;
               Normal /= norm_2(Normal);

               rT1 = Normal;
               if ( fabs(Normal(2))  > 1e-5) {
                  rT1(2) = (-pow(Normal(0),2) - pow(Normal(1),2) )  / Normal(2);
                  rT1 /= norm_2(rT1);
               } else {
                  rT1(0) = 0; rT1(1) = 0; rT1(2) = 1;
               }

               rT2 = ZeroVector(3);
               rT2(0) =  Normal(1) * rT1(2) - Normal(2)*rT1(1);
               rT2(1) = -Normal(0) * rT1(2) + Normal(2)*rT1(0);
               rT2(2) =  Normal(0) * rT1(1) - Normal(1)*rT1(0);

               KRATOS_CATCH("")
            }

    //************************************************************************************
    //************************************************************************************

    //Sphere
    void CreateBoundingBoxBoundaryMesh(ModelPart& rModelPart, int linear_partitions = 4, int angular_partitions = 4 ) override
    {
      KRATOS_TRY

      //std::cout<<" Create Sphere Mesh NODES "<<std::endl;

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

      std::vector<PointType> FacePoints;

      for(int k=0; k<=angular_partitions; k++)
	{
	  alpha = (3.14159262 * k)/(double)angular_partitions;

	  //vector of rotation
	  RotationAxis = DirectionX * alpha;
	  Quaternion   = QuaternionType::FromRotationVector(RotationAxis);

	  RotatedDirectionY = DirectionY;

	  Quaternion.RotateVector3(RotatedDirectionY);

	  //std::cout<<" Rotated "<<RotatedDirectionY<<" alpha "<<alpha<<std::endl;

	  noalias(BasePoint) = mBox.Center + mBox.Radius * RotatedDirectionY;

	  FacePoints.push_back(BasePoint);
	}

      unsigned int size = FacePoints.size();

      for(int k=1; k<=angular_partitions; k++)
	{
	  alpha = (2 * 3.14159262 * k)/((double)angular_partitions+1);

	  //vector of rotation
	  RotationAxis = DirectionY * alpha;
	  Quaternion   = QuaternionType::FromRotationVector(RotationAxis);

	  for(unsigned int j=1; j<size-1; j++)
	    {
	      RotatedDirectionY = FacePoints[j]-mBox.Center;

	      Quaternion.RotateVector3(RotatedDirectionY);

	      //std::cout<<" Rotated "<<RotatedDirectionY<<" alpha "<<alpha<<std::endl;

	      noalias(BasePoint) = mBox.Center + RotatedDirectionY;

	      FacePoints.push_back(BasePoint);
	    }

	}

        //create modelpart nodes
        for(unsigned int i=0; i<FacePoints.size(); i++)
	  {

	    NodeId += 1;

	    //std::cout<<" Node["<<NodeId<<"] "<<FacePoints[i]<<std::endl;

	    NodeType::Pointer pNode = this->CreateNode(rModelPart, FacePoints[i], NodeId);

	    pNode->Set(RIGID,true);

	    rModelPart.AddNode( pNode );

	    //get boundary model parts ( temporary implementation )
	    for(unsigned int j=0; j<BoundaryModelPartsName.size(); j++)
	      (pMainModelPart->GetSubModelPart(BoundaryModelPartsName[j])).AddNode( pNode );
	    //get boundary model parts ( temporary implementation )

	  }


      //std::cout<<" Nodes Added "<<NodeId-InitialNodeId<<std::endl;
      this->CreateSphereBoundaryMesh(rModelPart, InitialNodeId, angular_partitions);

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
        return "SphereBoundingBox";
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
      PointType Projection(3);
      Projection = rRadius * ( (rPoint-mBox.Center)/ norm_2(rPoint-mBox.Center) ) + mBox.Center;


      //2.-compute gap
      double GapNormal = 0;
      if( norm_2(mBox.Center-rPoint) <= rRadius ){
	GapNormal = (-1) * norm_2(rPoint - Projection);
      }
      else{
	GapNormal = norm_2(Projection - rPoint);
      }

      GapNormal *= mBox.Convexity;

      if(GapNormal<0)
	return true;
      else
	return false;


      KRATOS_CATCH( "" )

   }


    //************************************************************************************
    //************************************************************************************

    //Sphere (note: box position has been updated previously)
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
      PointType Projection = mBox.Radius * ( (rPoint-mBox.Center) / norm_2(rPoint-mBox.Center) ) + mBox.Center;

      //2.-compute contact normal
      rNormal = (Projection-mBox.Center)/mBox.Radius;

      rNormal *= mBox.Convexity;

      //3.-compute gap
      if( norm_2(mBox.Center-rPoint) <= mBox.Radius ){
	rGapNormal = (-1) * norm_2(rPoint - Projection);
      }
      else{
	rGapNormal = norm_2(Projection - rPoint);
      }

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

    void CreateSphereBoundaryMesh(ModelPart& rModelPart, const unsigned int& rInitialNodeId, const unsigned int& angular_partitions )
    {

      KRATOS_TRY

      //std::cout<<" Create Sphere Mesh ELEMENTS "<<std::endl;

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


      //GEOMETRY:
      GeometryType::Pointer  pFace;
      ElementType::Pointer   pElement;

      //PROPERTIES:
      int number_of_properties = rModelPart.NumberOfProperties();
      Properties::Pointer pProperties = Kratos::make_shared<Properties>(number_of_properties);

      //Shere numeration matrix:
      matrix<int> Connectivities(angular_partitions+2,angular_partitions+1);
      int counter = 1;
      for(unsigned int i=0; i<=angular_partitions; i++)
	{
	  Connectivities(i,0) = 1;
	  if( i == 0){
	    for(unsigned int j=1; j<=angular_partitions; j++)
	      {
		counter +=1;
		Connectivities(i,j) = counter;
	      }
	  }
	  else{
	    Connectivities(i,angular_partitions) = angular_partitions+1;
	    for(unsigned int j=1; j<angular_partitions; j++)
	      {
		counter +=1;
		Connectivities(i,j) = counter;
	      }
	  }

	}
      for(unsigned int i=0; i<=angular_partitions; i++)
	{
	  Connectivities(angular_partitions+1,i) = Connectivities(0,i);
	}

      //std::cout<<" Connectivities "<<Connectivities<<std::endl;

      Vector FaceNodesIds = ZeroVector(4);

      for( unsigned int i=0; i<=angular_partitions; i++ )
	{
	  ElementId += 1;

	  FaceNodesIds[0] = rInitialNodeId + Connectivities(i,0) ;
	  FaceNodesIds[1] = rInitialNodeId + Connectivities(i,1);
	  FaceNodesIds[2] = rInitialNodeId + Connectivities(i+1,1);

	  GeometryType::PointsArrayType FaceNodes;
	  FaceNodes.reserve(3);

	  //NOTE: when creating a PointsArrayType
	  //important ask for pGetNode, if you ask for GetNode a copy is created
	  //if a copy is created a segmentation fault occurs when the node destructor is called

	  for(unsigned int j=0; j<3; j++)
	    FaceNodes.push_back(rModelPart.pGetNode(FaceNodesIds[j]));

	  //std::cout<<" ElementA "<<FaceNodesIds<<std::endl;

	  pFace    = Kratos::make_shared<Triangle3D3<NodeType> >( FaceNodes );
	  pElement = Kratos::make_shared<Element>(ElementId, pFace, pProperties);

	  rModelPart.AddElement(pElement);
	  pElement->Set(ACTIVE,false);
	  pComputingModelPart->AddElement(pElement);

	}

      counter = 0;
      for( unsigned int i=1; i<angular_partitions-1; i++)
	{
	  for( unsigned int k=0; k<=angular_partitions; k++ )
	    {
	      ElementId += 1;

	      FaceNodesIds[0] = rInitialNodeId + Connectivities(k,i);
	      FaceNodesIds[1] = rInitialNodeId + Connectivities(k,i+1);
	      FaceNodesIds[2] = rInitialNodeId + Connectivities(k+1,i+1);
	      FaceNodesIds[3] = rInitialNodeId + Connectivities(k+1,i);

	      GeometryType::PointsArrayType FaceNodes;
	      FaceNodes.reserve(4);

	      //NOTE: when creating a PointsArrayType
	      //important ask for pGetNode, if you ask for GetNode a copy is created
	      //if a copy is created a segmentation fault occurs when the node destructor is called

	      for(unsigned int j=0; j<4; j++)
		FaceNodes.push_back(rModelPart.pGetNode(FaceNodesIds[j]));

	      //std::cout<<" ElementB "<<FaceNodesIds<<std::endl;

	      pFace    = Kratos::make_shared<Quadrilateral3D4<NodeType> >(FaceNodes);
	      pElement = Kratos::make_shared<Element>(ElementId, pFace, pProperties);

	      rModelPart.AddElement(pElement);
	      pElement->Set(ACTIVE,false);
	      pComputingModelPart->AddElement(pElement);

	    }
	}


      for( unsigned int i=0; i<=angular_partitions; i++ )
	{
	  ElementId += 1;

	  FaceNodesIds[0] = rInitialNodeId + Connectivities(i,angular_partitions-1);
	  FaceNodesIds[1] = rInitialNodeId + Connectivities(i,angular_partitions);
	  FaceNodesIds[2] = rInitialNodeId + Connectivities(i+1,angular_partitions-1);

	  GeometryType::PointsArrayType FaceNodes;
	  FaceNodes.reserve(3);

	  //NOTE: when creating a PointsArrayType
	  //important ask for pGetNode, if you ask for GetNode a copy is created
	  //if a copy is created a segmentation fault occurs when the node destructor is called

	  for(unsigned int j=0; j<3; j++)
	    FaceNodes.push_back(rModelPart.pGetNode(FaceNodesIds[j]));


	  //std::cout<<" ElementC "<<FaceNodesIds<<std::endl;

	  pFace    = Kratos::make_shared<Triangle3D3<NodeType>>(FaceNodes);
	  pElement = Kratos::make_shared<Element>(ElementId, pFace, pProperties);

	  rModelPart.AddElement(pElement);
	  pElement->Set(ACTIVE,false);
	  pComputingModelPart->AddElement(pElement);

	}


      //std::cout<<" Create Sphere Mesh ELEMENTS Created "<<std::endl;

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


}; // Class SphereBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_SPHERE_BOUNDING_BOX_H_INCLUDED  defined
