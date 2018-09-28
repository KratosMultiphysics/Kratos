//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_TUBE_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_TUBE_BOUNDING_BOX_H_INCLUDED

// External includes

// System includes

// Project includes
#include "custom_bounding/spatial_bounding_box.hpp"
#include "geometries/quadrilateral_3d_4.h"
#include "custom_utilities/spline_curve_utilities.hpp"

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

class TubeBoundingBox
  : public SpatialBoundingBox
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TubeBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( TubeBoundingBox );

    //typedef BoundedVector<double, 3>                     PointType;
    typedef array_1d<double, 3>                             PointType;
    typedef Node<3>                                          NodeType;
    typedef ModelPart::NodesContainerType          NodesContainerType;
    typedef NodesContainerType::Pointer     NodesContainerTypePointer;
    typedef ModelPart::ElementsContainerType    ElementsContainerType;
    typedef Quaternion<double>                         QuaternionType;
    typedef ModelPart::ElementType                        ElementType;
    typedef ElementType::GeometryType                    GeometryType;

    //definitions for spatial search
    typedef NodeType::Pointer                         NodePointerType;
    typedef std::vector<NodePointerType>        NodePointerTypeVector;
    typedef std::vector<NodeType>                      NodeTypeVector;
    typedef NodePointerTypeVector::iterator       NodePointerIterator;
    typedef std::vector<double>                        DistanceVector;
    typedef std::vector<double>::iterator            DistanceIterator;
    typedef Bucket<3, NodeType, NodePointerTypeVector, NodePointerType, NodePointerIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> >            KdtreeType; //Kdtree


public:

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TubeBoundingBox() : SpatialBoundingBox()
    {
      KRATOS_TRY

      std::cout<< "Calling Tube BBX empty constructor" <<std::endl;

      KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************

    TubeBoundingBox(ModelPart& rModelPart, Parameters CustomParameters)
    {

      KRATOS_TRY

      Parameters DefaultParameters( R"(
            {
                "parameters_list":[{
                   "radius": 0.0,
                   "convexity": 1
                 }],
                 "velocity": [0.0, 0.0, 0.0]

            }  )" );


      //validate against defaults -- this also ensures no type mismatch
      CustomParameters.ValidateAndAssignDefaults(DefaultParameters);

      if(CustomParameters["parameters_list"].IsArray() == true && CustomParameters["parameters_list"].size() != 1)
        {
	  KRATOS_THROW_ERROR(std::runtime_error,"paramters_list for the Tube BBX must contain only one term",CustomParameters.PrettyPrintJsonString());
        }

      mBox.Initialize();

      Parameters BoxParameters = CustomParameters["parameters_list"][0];

      mBox.Radius = BoxParameters["radius"].GetDouble();

      mBox.Velocity[0] = CustomParameters["velocity"][0].GetDouble();
      mBox.Velocity[1] = CustomParameters["velocity"][1].GetDouble();
      mBox.Velocity[2] = CustomParameters["velocity"][2].GetDouble();

      mBox.Convexity = BoxParameters["convexity"].GetInt();

      mBox.SetInitialValues();

      mRigidBodyCenterSupplied = false;

      this->CreateKnotsList(rModelPart);

      KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************

    TubeBoundingBox(ModelPart& rModelPart, double Radius, int Convexity)
    {

      KRATOS_TRY


      mBox.Initialize();

      mBox.Radius = Radius;

      mBox.Convexity = Convexity;

      mBox.SetInitialValues();

      mRigidBodyCenterSupplied = false;

      this->CreateKnotsList(rModelPart);

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    /// Assignment operator.
    TubeBoundingBox& operator=(TubeBoundingBox const& rOther)
    {
      KRATOS_TRY

      SpatialBoundingBox::operator=(rOther);

      mKnotsList = rOther.mKnotsList;
      mpKnotsKdtree = rOther.mpKnotsKdtree;

      return *this;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    /// Copy constructor.
    TubeBoundingBox(TubeBoundingBox const& rOther)
      :SpatialBoundingBox(rOther)
      ,mKnotsList(rOther.mKnotsList)
      ,mpKnotsKdtree(rOther.mpKnotsKdtree)
    {
    }


    //**************************************************************************
    //**************************************************************************

    /// Destructor.
    virtual ~TubeBoundingBox() {};


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

      double TubeRadius = mBox.Radius-Radius;

      //outside
      if( mBox.Convexity == 1)
	TubeRadius *= 1.25; //increase the bounding box

      //inside
      if( mBox.Convexity == -1)
       	TubeRadius *= 0.75; //decrease the bounding box

      is_inside = ContactSearch(rPoint, TubeRadius);

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

    //Tube
    void CreateBoundingBoxBoundaryMesh(ModelPart& rModelPart, int linear_partitions = 4, int angular_partitions = 4 ) override
    {
      KRATOS_TRY

      //std::cout<<" Create Tube Mesh "<<std::endl;

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


      SplineCurveUtilities::SplineType Spline;
      double PointId = 0;
      int KnotId  = 0;
      int number_of_segments = linear_partitions;

      PointType Point(3);

      // Create Axis generatrix
      // for(int i=1; i<mKnotsList.size()-2; i++)
      // 	{
      // 	  Point[0] = mKnotsList[i]->X();
      // 	  Point[1] = mKnotsList[i]->Y();
      // 	  Point[2] = mKnotsList[i]->Z();

      // 	  NodeId += 1;
      // 	  NodeType::Pointer pNode = this->CreateNode(rModelPart, Point, NodeId);
      // 	  pNode->Set(RIGID,true);
      // 	  rModelPart.AddNode( pNode );

      //     //get boundary model parts ( temporary implementation )
      // 	  for(unsigned int j=0; j<BoundaryModelPartsName.size(); j++)
      // 	    (pMainModelPart->GetSubModelPart(BoundaryModelPartsName[j])).AddNode( pNode );
      //     //get boundary model parts ( temporary implementation )

      //     KnotId = mKnotsList[i]->Id();

      //     mSplineCurveUtilities.SetSpline(Spline, mKnotsList, KnotId);

      //     PointId = KnotId;

      // 	  if( i == mKnotsList.size()-3 ) //last segment write last point only
      // 	    number_of_segments = 1;


      // 	  for(int j=0; j<number_of_segments; j++)
      // 	    {
      // 	      PointId = j/double(linear_partitions);

      // 	      Point = ZeroVector(3);
      // 	      Point = mSplineCurveUtilities.PointOnCurve(Point, Spline, PointId);

      //         NodeId += 1;
      //         NodeType::Pointer pNode = this->CreateNode(rModelPart, Point, NodeId);
      //         pNode->Set(RIGID,true);
      //         rModelPart.AddNode( pNode );

      //         //get boundary model parts ( temporary implementation )
      // 	      for(unsigned int j=0; j<BoundaryModelPartsName.size(); j++)
      // 		(pMainModelPart->GetSubModelPart(BoundaryModelPartsName[j])).AddNode( pNode );
      // 	      //get boundary model parts ( temporary implementation )

      //       }

      // 	}
      // number_of_segments = linear_partitions;
      // InitialNodeId = NodeId;

      PointType DirectionX(3);
noalias(DirectionX) = ZeroVector(3);
      PointType DirectionY(3);
      noalias(DirectionY) = ZeroVector(3);
      PointType DirectionZ(3);
      noalias(DirectionZ) = ZeroVector(3);

      PointType BasePoint(3);
      PointType RotationAxis(3);
      PointType RotatedDirectionY(3);

      double alpha = 0;
      QuaternionType Quaternion;


      for(unsigned int i=1; i<mKnotsList.size()-2; i++)
	{
	  Point[0] = mKnotsList[i]->X();
	  Point[1] = mKnotsList[i]->Y();
	  Point[2] = mKnotsList[i]->Z();

          KnotId = mKnotsList[i]->Id();

          mSplineCurveUtilities.SetSpline(Spline, mKnotsList, KnotId);

          PointId = KnotId;

  	  if( i == mKnotsList.size()-3 ) //last segment write last point only
	    number_of_segments = 1;

	  for(int j=0; j<number_of_segments; j++)
	    {

	      PointId = j/double(linear_partitions);

	      Point = ZeroVector(3);
	      Point = mSplineCurveUtilities.PointOnCurve(Point, Spline, PointId);

              DirectionX = mSplineCurveUtilities.PointOnCurveFirstDerivative(Spline, PointId);

              this->CalculateOrthonormalBase(DirectionX, DirectionY, DirectionZ);


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
	}

      //std::cout<<" Create Tube Mesh NODES "<<std::endl;

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
        return "TubeBoundingBox";
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

    std::vector<NodeType::Pointer> mKnotsList;

    KdtreeType* mpKnotsKdtree;

    SplineCurveUtilities mSplineCurveUtilities;

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

    void CreateKnotsList(ModelPart& rModelPart)
    {
      KRATOS_TRY
      //1.-Create Generatrix
      NodePointerTypeVector GeneratrixPoints;

      //Set generatrix control points for a given set of two noded line conditions
      unsigned int id = 0; //start with 0;

      NodePointerType GeneratrixPoint;

      for(ElementsContainerType::iterator ie = rModelPart.ElementsBegin(); ie!=rModelPart.ElementsEnd(); ++ie)
	{
	  if( ie->GetGeometry().size() > 1){

	    GeneratrixPoint = Kratos::make_shared<NodeType>(id, ie->GetGeometry()[0].X(), ie->GetGeometry()[0].Y(), ie->GetGeometry()[0].Z());
	    GeneratrixPoints.push_back( GeneratrixPoint );

	    //std::cout<<" Point ["<<ie->GetGeometry()[0].X()<<", "<<ie->GetGeometry()[0].Y()<<", "<<ie->GetGeometry()[0].Z()<<"] "<<std::endl;

	    id++;
	  }
	}


      ElementsContainerType::iterator LastElement = rModelPart.ElementsEnd()-1;
      int num_nodes = LastElement->GetGeometry().size()-1;

      GeneratrixPoint = Kratos::make_shared<NodeType>(id,LastElement->GetGeometry()[num_nodes].X(),LastElement->GetGeometry()[num_nodes].Y(),LastElement->GetGeometry()[num_nodes].Z());
      GeneratrixPoints.push_back( GeneratrixPoint );

      std::cout<<"  [DEFINED BY:"<<GeneratrixPoints.size()<<" control points]"<<std::endl;

      //2.-Create Arch Length Parametrized Spline Curve
      int SegmentsNumber = 250;
      mSplineCurveUtilities.CreateParametrizedCurve(GeneratrixPoints, mKnotsList, SegmentsNumber);

      // for(unsigned int i=0; i<mKnotsList.size(); i++)
      // 	std::cout<<" mKnotsListI ["<<i<<"]: "<<mKnotsList[i]->Id()<<std::endl;

      std::cout<<"  [REDEFINED WITH: "<<mKnotsList.size()<<" knots]"<<std::endl;


      this->CreateKnotsKdtree();


      KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void CreateKnotsKdtree()
    {
      KRATOS_TRY


      //creating an auxiliary list for the pre integration points
      unsigned int   bucket_size = 20;

      NodePointerTypeVector KnotsList = mKnotsList;

      unsigned int knots_begin = 1;
      unsigned int knots_end   = 3;

      mpKnotsKdtree = new KdtreeType(KnotsList.begin()+knots_begin,KnotsList.end()-knots_end,bucket_size);


      KRATOS_CATCH("")
    }


    //************************************************************************************
    //************************************************************************************


    bool ContactSearch(const PointType& rPoint, const double& rRadius)
    {

      KRATOS_TRY

      //1.-compute point projection
      PointType Projection(3);
      Projection = mSplineCurveUtilities.CalculatePointProjection( rPoint, *mpKnotsKdtree, mKnotsList, Projection );

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

    //Tube (note: box position has been updated previously)
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
      PointType Projection(3);
      Projection = mSplineCurveUtilities.CalculatePointProjection( rPoint, *mpKnotsKdtree, mKnotsList, Projection );

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

      //std::cout<<" Create Tube Mesh ELEMENTS "<<std::endl;

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


      //std::cout<<" Create Tube Mesh ELEMENTS CREATED "<<std::endl;

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


}; // Class TubeBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_TUBE_BOUNDING_BOX_H_INCLUDED  defined
