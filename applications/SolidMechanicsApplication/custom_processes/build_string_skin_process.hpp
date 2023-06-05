//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               March 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_BUILD_STRING_SKIN_PROCESS_H_INCLUDED)
#define  KRATOS_BUILD_STRING_SKIN_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "geometries/point_3d.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "utilities/beam_math_utilities.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/**  Build Beam or Tube from a generatrix
     set of control points which define a generatrix curve
     set of beam nodes which define the beam generatrix
     radius: define the walls of the tube respect to the generatrix
*/

class BuildStringSkinProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef array_1d<double, 3>                             PointType;
    typedef ModelPart::NodeType                              NodeType;
    typedef ModelPart::NodesContainerType          NodesContainerType;
    typedef ModelPart::ElementsContainerType    ElementsContainerType;
    typedef NodesContainerType::Pointer     NodesContainerPointerType;
    typedef ModelPart::ConditionType                    ConditionType;
    typedef ModelPart::PropertiesType                  PropertiesType;
    typedef ConditionType::GeometryType                  GeometryType;
    typedef Triangle3D3<NodeType>                      Triangle3DType;
    typedef Point3D<NodeType>                             Point3DType;
    typedef Quadrilateral3D4<NodeType>            Quadrilateral3DType;

    typedef BeamMathUtils<double>                   BeamMathUtilsType;
    typedef Quaternion<double>                         QuaternionType;

    typedef GlobalPointersVector<Node >         NodeWeakPtrVectorType;
    typedef GlobalPointersVector<Element>       ElementWeakPtrVectorType;
    typedef GlobalPointersVector<Condition>   ConditionWeakPtrVectorType;
    /// Pointer definition of BuildStringSkinProcess
    KRATOS_CLASS_POINTER_DEFINITION(BuildStringSkinProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    BuildStringSkinProcess(ModelPart& rModelPart,
			   unsigned int sides,
			   double radius
			   ) : Process() , mrModelPart(rModelPart), mSides(sides), mRadius(radius)
    {
        KRATOS_TRY

	mMaxId = GetMaxNodeId(mrModelPart.GetParentModelPart());

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~BuildStringSkinProcess() override {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Execute method is used to execute the BuildStringSkinProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY

	CreateGeneratrix();
	CreateSkinNodes();
	CreateSkinElements();

	TransferSkinToOutput();

	// set nodes to RIGID and ACTIVE
	for (ModelPart::NodeIterator i = mrModelPart.NodesBegin(); i != mrModelPart.NodesEnd(); ++i)
	{
	  (i)->Set(RIGID,true);
	  (i)->Set(ACTIVE,true);
	}

	// set elements to RIGID and ACTIVE
	for (ModelPart::ConditionIterator i = mrModelPart.ConditionsBegin(); i != mrModelPart.ConditionsEnd(); ++i)
	{
	  (i)->Set(RIGID,true);
	  (i)->Set(ACTIVE,true);
	}

        KRATOS_CATCH("")
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    void ExecuteBeforeSolutionLoop() override
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {

        KRATOS_TRY

	SetInActiveFlag(mrModelPart);

	KRATOS_CATCH("")

    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {

        KRATOS_TRY

	MoveSkinNodes();

        KRATOS_CATCH("")

    }


    /// this function will be executed at every time step BEFORE  writing the output
    void ExecuteBeforeOutputStep() override
    {

        KRATOS_TRY

	SetActiveFlag(mrModelPart);

        KRATOS_CATCH("")
    }


    /// this function will be executed at every time step AFTER writing the output
    void ExecuteAfterOutputStep() override
    {
        KRATOS_TRY


        KRATOS_CATCH("")
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    void ExecuteFinalize() override
    {
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
        return "BuildStringSkinProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BuildStringSkinProcess";
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

    /// Copy constructor.
    BuildStringSkinProcess(BuildStringSkinProcess const& rOther);

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

    ModelPart& mrModelPart;

    NodesContainerType mGeneratrixNodes;

    unsigned int mSides;

    unsigned int mMaxId;

    double mRadius;

    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{


    //************************************************************************************
    //************************************************************************************

    void CreateGeneratrix()
    {
      KRATOS_TRY

      //Set generatrix control points for a given set of two noded line conditions
      //unsigned int id = 0; //start with 0;

      //NEIGHBOURS SEARCH
      SearchNeighbours();

      //SEARCH INITIAL NODE OF THE STRING
      Node::Pointer Starter;
      for(auto i_node(mrModelPart.NodesBegin()); i_node != mrModelPart.NodesEnd(); ++i_node)
      {
        NodeWeakPtrVectorType& nNodes = i_node->GetValue(NEIGHBOUR_NODES);
        if( nNodes.size() <= 1 )
          Starter = *i_node.base();
      }

      //SEARCH CONSECUTIVE NODES
      unsigned int previous_id = 0;
      unsigned int current_id  = 0;
      Element::Pointer  CurrentElement;
      for(unsigned int i=0; i<mrModelPart.NumberOfNodes(); i++)
      {

        //std::cout<<" Node ("<<Starter->Id()<<") "<<std::endl;

        current_id = Starter->Id();

        mGeneratrixNodes.push_back( Starter );

        ElementWeakPtrVectorType& nElements = Starter->GetValue(NEIGHBOUR_ELEMENTS);

        for(auto& i_nelem : nElements)
        {
          Element::GeometryType& rGeometry = i_nelem.GetGeometry();

          bool selected = true;
          for(unsigned int j = 0; j < rGeometry.size(); j++)
          {
            if( rGeometry[j].Id() == previous_id )
            {
              selected = false;
            }
          }

          if( selected ){

            for(unsigned int j = 0; j < rGeometry.size(); j++)
            {
              if( rGeometry[j].Id() != current_id )
              {
                previous_id = Starter->Id();
                Starter = rGeometry(j);
              }
            }
          }

        }

      }


      //SET MEAN_RADIUS TO NODES
      /*for(ElementsContainerType::iterator ie = mrModelPart.ElementsBegin(); ie!=mrModelPart.ElementsEnd(); ie++)
	{
	  PointsArrayType& Vertices = ie->GetGeometry().Points();

	  //set radius to nodes
	  PropertiesType& Properties = ie->GetProperties();
	  Radius = Properties[MEAN_RADIUS];

	  Vertices(0)->SetValue(MEAN_RADIUS, Radius);
	}

      ElementsContainerType::iterator LastElement = mrModelPart.ElementsEnd()-1;
      int num_nodes = LastElement->GetGeometry().size()-1;

      PointsArrayType& Vertices = LastElement->GetGeometry().Points();

      //set radius to nodes
      PropertiesType& Properties = LastElement->GetProperties();
      Radius = Properties[MEAN_RADIUS];
      Vertices(num_nodes)->SetValue(MEAN_RADIUS, Radius);*/

      //REMOVE STRING NODES FROM GIVEN SKIN MODELPART AFTER CONSTRUCTING THE GENERATRIX
      SetActiveFlag(mrModelPart,TO_ERASE);

      mrModelPart.RemoveNodes();
      mrModelPart.RemoveElements();

      std::cout<<"  [String_builder] [Defined by "<<mGeneratrixNodes.size()<<" control points]"<<std::endl;

      KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void CreateSkinNodes()
    {
      KRATOS_TRY

      unsigned int node_id = mMaxId;

      PointType Point0     = ZeroVector(3);
      PointType Point      = ZeroVector(3);
      PointType BasePoint  = ZeroVector(3);

      PointType DirectionX = ZeroVector(3);
      PointType DirectionY = ZeroVector(3);
      PointType DirectionZ = ZeroVector(3);

      PointType RotatedDirectionX = ZeroVector(3);

      ModelPart::NodesContainerType::iterator nodes_begin = mGeneratrixNodes.begin();

      double Radius = 0;
      for(unsigned int i=0; i<mGeneratrixNodes.size(); i++)
	{

	  Point0[0] = (nodes_begin+i)->X0();
	  Point0[1] = (nodes_begin+i)->Y0();
	  Point0[2] = (nodes_begin+i)->Z0();

	  Point[0] = (nodes_begin+i)->X();
	  Point[1] = (nodes_begin+i)->Y();
	  Point[2] = (nodes_begin+i)->Z();

	  // rotations
	  double alpha = 0;
	  Matrix Q = ZeroMatrix(3,3); //rotation along the local axis X

	  //vector of beam section rotation
	  array_1d<double,3>& NodeRotation = (nodes_begin+i)->FastGetSolutionStepValue( ROTATION );

	  PointType SectionRotation = ZeroVector(3);

	  for( unsigned int j=0; j<3; j++ )
	    {
	      SectionRotation[j] = NodeRotation[j];
	    }

	  QuaternionType SectionQuaternion;
	  SectionQuaternion = QuaternionType::FromRotationVector(SectionRotation);

	  QuaternionType RotationQuaternion;

	  PointType RotationAxis = ZeroVector(3);

	  PointType DirectionZ1 = ZeroVector(3);
	  PointType DirectionZ2 = ZeroVector(3);

	  PointType DirectionEllipse = ZeroVector(3);
	  double RadiusCorrection = 1;

	  if( i == mGeneratrixNodes.size()-1 ){
	    BasePoint[0] = (nodes_begin+(i-1))->X0();
	    BasePoint[1] = (nodes_begin+(i-1))->Y0();
	    BasePoint[2] = (nodes_begin+(i-1))->Z0();
	    DirectionZ = Point0 - BasePoint;
	  }
	  else if( i == 0 ){
	    BasePoint[0] = (nodes_begin+(i+1))->X0();
	    BasePoint[1] = (nodes_begin+(i+1))->Y0();
	    BasePoint[2] = (nodes_begin+(i+1))->Z0();
	    DirectionZ = BasePoint - Point0;
	  }
	  else{
	    BasePoint[0] = (nodes_begin+(i-1))->X0();
	    BasePoint[1] = (nodes_begin+(i-1))->Y0();
	    BasePoint[2] = (nodes_begin+(i-1))->Z0();

	    DirectionZ1  = (Point0 - BasePoint);
	    DirectionZ1 /= norm_2(DirectionZ1);

	    BasePoint[0] = (nodes_begin+(i+1))->X0();
	    BasePoint[1] = (nodes_begin+(i+1))->Y0();
	    BasePoint[2] = (nodes_begin+(i+1))->Z0();

	    DirectionZ2  = (BasePoint - Point0);
	    DirectionZ2 /= norm_2(DirectionZ2);

	    DirectionZ   = DirectionZ1 + DirectionZ2;

	    if(norm_2(DirectionZ))
	      DirectionZ  /= norm_2(DirectionZ);
	    else
	      DirectionZ = DirectionZ1;

	    DirectionEllipse = DirectionZ1 - DirectionZ2;

	    if(norm_2(DirectionEllipse))
	      DirectionEllipse/=norm_2(DirectionEllipse);

	    RadiusCorrection  = inner_prod(DirectionZ,DirectionZ1);
	    RadiusCorrection += inner_prod(DirectionZ,DirectionZ2);
	    RadiusCorrection *= 0.5;
	    if(RadiusCorrection)
	      RadiusCorrection  = 1.0/RadiusCorrection;
	    else
	      RadiusCorrection  = 1.0;
	  }


	  BeamMathUtilsType::CalculateLocalAxesVectors(DirectionZ,DirectionX,DirectionY);

	  if( DirectionZ[0] == 0 && DirectionZ[1] == 0 && DirectionZ[2] == 1){ //if e3 change the orthornormal base
	    PointType Temp  =  DirectionX;
	    DirectionX = DirectionY;
	    DirectionY = (-1) * Temp;
	  }

	  double EllipsoidalCorrection = 1;

	  for(unsigned int k=0; k<mSides; k++)
	    {
              alpha = (2.0 * Globals::Pi * k)/double(mSides) + 0.25 * Globals::Pi;

	      //vector of rotation
	      RotationAxis = DirectionZ * alpha;

	      RotationQuaternion = QuaternionType::FromRotationVector(RotationAxis);

	      RotatedDirectionX = DirectionX;

	      RotationQuaternion.RotateVector3(RotatedDirectionX);


	      EllipsoidalCorrection = inner_prod(DirectionEllipse,RotatedDirectionX);
	      EllipsoidalCorrection = 1 + ( RadiusCorrection - 1 ) * (EllipsoidalCorrection * EllipsoidalCorrection);

	      //Add four points along the circular base of the tube
	      node_id += 1;

	      //Radius = (nodes_begin+i)->GetValue(MEAN_RADIUS);
	      Radius = mRadius;

	      BasePoint =  Radius * EllipsoidalCorrection * RotatedDirectionX;

	      SectionQuaternion.RotateVector3(BasePoint);

	      BasePoint += Point;

	      // std::cout<<" DirectionZ "<<DirectionZ<<std::endl;
	      // std::cout<<" alpha "<<alpha<<"  cos "<<Q(1,1)<<std::endl;
	      // std::cout<<" Rotated "<<RotatedDirectionX<<" alpha "<<alpha<<std::endl;
	      // std::cout<<" EllipsoidalCorrection "<<EllipsoidalCorrection<<std::endl;
	      // std::cout<<" Radius "<<Radius<<std::endl;
	      // std::cout<<" Base Point ["<<node_id<<"]"<<BasePoint<<std::endl;

	      mrModelPart.AddNode(this->CreateNode(mrModelPart.GetParentModelPart(), BasePoint, node_id));

	    }


	}

      KRATOS_CATCH( "" )

    }

    //************************************************************************************
    //************************************************************************************
    void MoveSkinNodes()
    {

      KRATOS_TRY

      //set new nodes position:

      int number_of_angles = mSides; //number of lines in radius

      PointType Point0     = ZeroVector(3);
      PointType Point      = ZeroVector(3);
      PointType BasePoint  = ZeroVector(3);

      PointType DirectionX = ZeroVector(3);
      PointType DirectionY = ZeroVector(3);
      PointType DirectionZ = ZeroVector(3);

      PointType RotatedDirectionX = ZeroVector(3);

      ModelPart::NodesContainerType::iterator nodes_begin = mGeneratrixNodes.begin();

      ModelPart::NodesContainerType::iterator skin_nodes_begin = mrModelPart.NodesBegin();

      //std::cout<<" Number of Nodes "<<mrModelPart.NumberOfNodes()<<std::endl;

      int counter = 0;
      double Radius = 0;
      for(unsigned int i=0; i<mGeneratrixNodes.size(); i++)
	{

	  Point0[0] = (nodes_begin+i)->X0();
	  Point0[1] = (nodes_begin+i)->Y0();
	  Point0[2] = (nodes_begin+i)->Z0();

	  Point[0] = (nodes_begin+i)->X();
	  Point[1] = (nodes_begin+i)->Y();
	  Point[2] = (nodes_begin+i)->Z();

	  // rotations
	  double alpha = 0;
	  Matrix Q(3,3);
	  noalias(Q) = ZeroMatrix(3,3); //rotation along the local axis X

	  //vector of beam section rotation
	  PointType& NodeRotation = (nodes_begin+i)->FastGetSolutionStepValue( ROTATION );

	  PointType SectionRotation = ZeroVector(3);

	  for( unsigned int j=0; j<3; j++ )
	    {
	      SectionRotation[j] = NodeRotation[j];
	    }

	  QuaternionType SectionQuaternion;
	  SectionQuaternion = QuaternionType::FromRotationVector(SectionRotation);

	  QuaternionType RotationQuaternion;

	  PointType RotationAxis;
	  PointType DirectionZ1;
	  PointType DirectionZ2;

	  PointType DirectionEllipse;
	  double RadiusCorrection = 1;

	  if( i == mGeneratrixNodes.size()-1 ){
	    BasePoint[0] = (nodes_begin+(i-1))->X0();
	    BasePoint[1] = (nodes_begin+(i-1))->Y0();
	    BasePoint[2] = (nodes_begin+(i-1))->Z0();
	    DirectionZ = Point0 - BasePoint;
	  }
	  else if( i == 0 ){
	    BasePoint[0] = (nodes_begin+(i+1))->X0();
	    BasePoint[1] = (nodes_begin+(i+1))->Y0();
	    BasePoint[2] = (nodes_begin+(i+1))->Z0();
	    DirectionZ = BasePoint - Point0;
	  }
	  else{
	    BasePoint[0] = (nodes_begin+(i-1))->X0();
	    BasePoint[1] = (nodes_begin+(i-1))->Y0();
	    BasePoint[2] = (nodes_begin+(i-1))->Z0();

	    DirectionZ1  = (Point0 - BasePoint);
	    DirectionZ1 /= norm_2(DirectionZ1);

	    BasePoint[0] = (nodes_begin+(i+1))->X0();
	    BasePoint[1] = (nodes_begin+(i+1))->Y0();
	    BasePoint[2] = (nodes_begin+(i+1))->Z0();

	    DirectionZ2  = (BasePoint - Point0);
	    DirectionZ2 /= norm_2(DirectionZ2);

	    DirectionZ   = DirectionZ1 + DirectionZ2;
	    DirectionZ  /= norm_2(DirectionZ);

	    DirectionEllipse = DirectionZ1 - DirectionZ2;

	    if(norm_2(DirectionEllipse))
	      DirectionEllipse/=norm_2(DirectionEllipse);

	    RadiusCorrection  = inner_prod(DirectionZ,DirectionZ1);
	    RadiusCorrection += inner_prod(DirectionZ,DirectionZ2);
	    RadiusCorrection *= 0.5;
	    RadiusCorrection  = 1.0/RadiusCorrection;
	  }


	  //std::cout<<" DirectionZ "<<DirectionZ<<std::endl;

	  BeamMathUtilsType::CalculateLocalAxesVectors(DirectionZ,DirectionX,DirectionY);

	  if( DirectionZ[0] == 0 && DirectionZ[1] == 0 && DirectionZ[2] == 1){ //if e3 change the orthornormal base
	    PointType Temp  =  DirectionX;
	    DirectionX = DirectionY;
	    DirectionY = (-1) * Temp;
	  }

	  //std::cout<<" OrthonormalBase "<<DirectionX<<" "<<DirectionY<<std::endl;

	  double EllipsoidalCorrection = 1;

	  for(int k=0; k<number_of_angles; k++)
	    {
	      alpha = (2.0 * Globals::Pi * k)/double(number_of_angles) + 0.25 * Globals::Pi;

	      //vector of rotation
	      RotationAxis = DirectionZ * alpha;

	      RotationQuaternion = QuaternionType::FromRotationVector(RotationAxis);

	      RotatedDirectionX = DirectionX;

	      RotationQuaternion.RotateVector3(RotatedDirectionX);

	      //std::cout<<" alpha "<<alpha<<"  cos "<<Q(1,1)<<std::endl;
	      //std::cout<<" Rotated "<<RotatedDirectionX<<" alpha "<<alpha<<std::endl;

	      EllipsoidalCorrection = inner_prod(DirectionEllipse,RotatedDirectionX);
	      EllipsoidalCorrection = 1 + ( RadiusCorrection - 1 ) * (EllipsoidalCorrection * EllipsoidalCorrection);

	      //Add four points along the circular base of the tube
	      //Radius = (nodes_begin+i)->GetValue(MEAN_RADIUS);
	      Radius = mRadius;

	      BasePoint = Radius * EllipsoidalCorrection * RotatedDirectionX;

	      SectionQuaternion.RotateVector3(BasePoint);

	      //set displacement velocity and acceleration START
	      PointType  RadiusVector = ZeroVector(3);
	      RadiusVector[0] = BasePoint[0];
	      RadiusVector[1] = BasePoint[1];
	      RadiusVector[2] = BasePoint[2];

	      //get coordinates
	      PointType PreviousPosition = ZeroVector(3);
	      PreviousPosition[0] = (skin_nodes_begin+counter)->X0();
	      PreviousPosition[1] = (skin_nodes_begin+counter)->Y0();
	      PreviousPosition[2] = (skin_nodes_begin+counter)->Z0();

	      BasePoint += Point;

	      (skin_nodes_begin+counter)->X() = BasePoint[0];
	      (skin_nodes_begin+counter)->Y() = BasePoint[1];
	      (skin_nodes_begin+counter)->Z() = BasePoint[2];


	      PointType& Displacement = (skin_nodes_begin+counter)->FastGetSolutionStepValue(DISPLACEMENT);

	      for(int j=0; j<3; j++)
		Displacement[j] = BasePoint[j]-PreviousPosition[j];

	      bool DynamicVariables = false;

	      if( DynamicVariables ){

		array_1d<double, 3 >&  Velocity            = (nodes_begin+i)->FastGetSolutionStepValue(VELOCITY);
		array_1d<double, 3 >&  Acceleration        = (nodes_begin+i)->FastGetSolutionStepValue(ACCELERATION);
		array_1d<double, 3 >&  AngularVelocity     = (nodes_begin+i)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
		array_1d<double, 3 >&  AngularAcceleration = (nodes_begin+i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION);

		Matrix SkewSymVariable(3,3);
		noalias(SkewSymVariable) = ZeroMatrix(3,3);
		PointType Variable;
		PointType AngularVariable;

		//********************
		//compute the skewsymmmetric tensor of the angular velocity
		BeamMathUtilsType::VectorToSkewSymmetricTensor(AngularVelocity, SkewSymVariable);

		//compute the contribution of the angular velocity to the velocity v = Wxr
		Variable = prod(SkewSymVariable,RadiusVector);

		(skin_nodes_begin+counter)->FastGetSolutionStepValue(VELOCITY) = Velocity + Variable;

		//********************

		//centripetal acceleration:

		//compute the skewsymmmetric tensor of the angular velocity
		BeamMathUtilsType::VectorToSkewSymmetricTensor(AngularVelocity, SkewSymVariable);

		AngularVariable = prod(SkewSymVariable,Variable); //ac = Wx(Wxr)

		//compute the skewsymmmetric tensor of the angular acceleration
		BeamMathUtilsType::VectorToSkewSymmetricTensor(AngularAcceleration, SkewSymVariable);

		//compute the contribution of the angular velocity to the velocity a = Axr
		Variable = prod(SkewSymVariable,RadiusVector);

		(skin_nodes_begin+counter)->FastGetSolutionStepValue(ACCELERATION) = Acceleration + Variable + AngularVariable;

		//set displacement velocity and acceleration END
	      }

	      counter++;
	    }


	}

      KRATOS_CATCH( "" )

    }

   //************************************************************************************
   //************************************************************************************


    void CreateSkinElements()
    {
      KRATOS_TRY

      //return this->CreateSkinTriangles();
      return this->CreateSkinQuadrilaterals();

      KRATOS_CATCH( "" )

    }


    //************************************************************************************
    //************************************************************************************

    void CreateSkinTriangles()
    {
      KRATOS_TRY

      unsigned int number_of_angles = mSides; //number of lines in radius x 2

      unsigned int wall_nodes_number_id = mMaxId; //used in the creation of the tube surface conditions

      //Triangles:

      // Create surface of the tube with triangular shell conditions
      unsigned int number_of_elements = (mrModelPart.Nodes().back().Id() - wall_nodes_number_id) - (number_of_angles - 1);

      //GEOMETRY:
      GeometryType::Pointer  pFace;
      ConditionType::Pointer pSkinCondition;

      //PROPERTIES:
      //int number_properties = mrModelPart.GetParentModelPart().NumberOfProperties();
      //Properties::Pointer pProperties = mrModelPart.GetParentModelPart().pGetProperties(number_properties-1);
      Properties::Pointer pProperties = mrModelPart.GetParentModelPart().pGetProperties(0);

      //Properties 0 in order to change the Id to 0 and then write tube elements in another layer
      /* ModelPart::PropertiesContainerType::ContainerType& PropertiesArray = mrModelPart.PropertiesArray(); */
      /* PropertiesArray[0] = Kratos::make_shared<PropertiesType>(0); */
      /* PropertiesArray[0]->Data() =PropertiesArray[1]->Data();    */
      /* Properties::Pointer pProperties = PropertiesArray[0]; */

      unsigned int condition_id = GetMaxConditionId(mrModelPart.GetParentModelPart());

      unsigned int counter = 1;

      std::vector<int> FaceNodesIds(3);

      for(unsigned int i=1; i<number_of_elements; i++)
	{
	  condition_id += 1;

	  if( counter < number_of_angles ){

	    //triangle 1
	    FaceNodesIds[0] = wall_nodes_number_id + i ;
	    FaceNodesIds[2] = wall_nodes_number_id + i + number_of_angles;
	    FaceNodesIds[1] = wall_nodes_number_id + i + number_of_angles + 1;

	    GeometryType::PointsArrayType FaceNodes1;
	    FaceNodes1.reserve(3);

	    //NOTE: when creating a PointsArrayType
	    //important ask for pGetNode, if you ask for GetNode a copy is created
	    //if a copy is created a segmentation fault occurs when the node destructor is called

	    for(unsigned int j=0; j<3; j++)
	      FaceNodes1.push_back(mrModelPart.pGetNode(FaceNodesIds[j]));

	    pFace = Kratos::make_shared<Triangle3DType>(FaceNodes1);

	    pSkinCondition = Kratos::make_intrusive<Condition>(condition_id, pFace, pProperties);

	    pSkinCondition->Set(ACTIVE,false);

	    //set to beam tube mesh
	    mrModelPart.AddCondition(pSkinCondition);

	    condition_id += 1;

	    //triangle 2
	    FaceNodesIds[0] = wall_nodes_number_id + i ;
	    FaceNodesIds[2] = wall_nodes_number_id + i + number_of_angles + 1;
	    FaceNodesIds[1] = wall_nodes_number_id + i + 1;

	    GeometryType::PointsArrayType FaceNodes2;
	    FaceNodes2.reserve(3);


	    //NOTE: when creating a PointsArrayType
	    //important ask for pGetNode, if you ask for GetNode a copy is created
	    //if a copy is created a segmentation fault occurs when the node destructor is called

	    for(unsigned int j=0; j<3; j++)
	      FaceNodes2.push_back(mrModelPart.pGetNode(FaceNodesIds[j]));

	    pFace = Kratos::make_shared<Triangle3DType>(FaceNodes2);

	    pSkinCondition = Kratos::make_intrusive<Condition>(condition_id, pFace, pProperties);

	    pSkinCondition->Set(ACTIVE,false);

	    //set to beam tube mesh
	    mrModelPart.AddCondition(pSkinCondition);

	    counter++;

	  }
	  else if( counter == number_of_angles ){

	    //triangle 1
	    FaceNodesIds[0] = wall_nodes_number_id + i ;
	    FaceNodesIds[1] = wall_nodes_number_id + i + 1 - number_of_angles;
	    FaceNodesIds[2] = wall_nodes_number_id + i + 1;

	    GeometryType::PointsArrayType FaceNodes1;
	    FaceNodes1.reserve(3);

	    for(unsigned int j=0; j<3; j++)
	      FaceNodes1.push_back(mrModelPart.pGetNode(FaceNodesIds[j]));

	    pFace = Kratos::make_shared<Triangle3DType>(FaceNodes1);

	    pSkinCondition = Kratos::make_intrusive<Condition>(condition_id, pFace, pProperties);

	    pSkinCondition->Set(ACTIVE,false);

	    //set to beam tube mesh
	    mrModelPart.AddCondition(pSkinCondition);

	    condition_id += 1;

	    //triangle 2
	    FaceNodesIds[0] = wall_nodes_number_id + i ;
	    FaceNodesIds[1] = wall_nodes_number_id + i + 1;
	    FaceNodesIds[2] = wall_nodes_number_id + i + number_of_angles;

	    GeometryType::PointsArrayType FaceNodes2;
	    FaceNodes2.reserve(3);

	    for(unsigned int j=0; j<3; j++)
	      FaceNodes2.push_back(mrModelPart.pGetNode(FaceNodesIds[j]));

	    pFace = Kratos::make_shared<Triangle3DType>(FaceNodes2);

	    pSkinCondition = Kratos::make_intrusive<Condition>(condition_id, pFace, pProperties);

	    pSkinCondition->Set(ACTIVE,false);

	    //set to beam tube mesh
	    mrModelPart.AddCondition(pSkinCondition);

	    counter = 1;
	  }


	}


      KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void CreateSkinQuadrilaterals()
    {
      KRATOS_TRY

      unsigned int number_of_angles = mSides; //number of lines in radius x 2

      unsigned int wall_nodes_number_id = mMaxId; //used in the creation of the tube surface conditions

      //Quadrilaterals:

      // Create surface of the tube with quadrilateral shell conditions
      unsigned int number_of_elements = (mrModelPart.Nodes().back().Id() - wall_nodes_number_id) - (number_of_angles - 1);

      //GEOMETRY:
      GeometryType::Pointer  pFace;
      ConditionType::Pointer pSkinCondition;

      //PROPERTIES:
      //int number_properties = mrModelPart.GetParentModelPart().NumberOfProperties();
      //Properties::Pointer pProperties = mrModelPart.GetParentModelPart().pGetProperties(number_properties-1);
      Properties::Pointer pProperties = mrModelPart.GetParentModelPart().pGetProperties(0);

      //Properties 0 in order to change the Id to 0 and then write tube elements in another layer
      /* ModelPart::PropertiesContainerType::ContainerType& PropertiesArray = mrModelPart.PropertiesArray(); */
      /* PropertiesArray[0] = Kratos::make_shared<PropertiesType>(0); */
      /* PropertiesArray[0]->Data() =PropertiesArray[1]->Data();    */
      /* Properties::Pointer pProperties = PropertiesArray[0]; */

      unsigned int condition_id = GetMaxConditionId(mrModelPart.GetParentModelPart());

      unsigned int counter = 1;

      std::vector<int> FaceNodesIds(4);


      for(unsigned int i=1; i<number_of_elements; i++)
	{
	  condition_id += 1;

	  if( counter < number_of_angles ){

	    FaceNodesIds[0] = wall_nodes_number_id + i ;
	    FaceNodesIds[1] = wall_nodes_number_id + i + number_of_angles;
	    FaceNodesIds[2] = wall_nodes_number_id + i + number_of_angles + 1;
	    FaceNodesIds[3] = wall_nodes_number_id + i + 1;

	    GeometryType::PointsArrayType FaceNodes;
	    FaceNodes.reserve(4);

	    //NOTE: when creating a PointsArrayType
	    //important ask for pGetNode, if you ask for GetNode a copy is created
	    //if a copy is created a segmentation fault occurs when the node destructor is called

	    for(unsigned int j=0; j<4; j++)
	      FaceNodes.push_back(mrModelPart.pGetNode(FaceNodesIds[j]));

	    pFace = Kratos::make_shared<Quadrilateral3DType>(FaceNodes);

	    pSkinCondition = Kratos::make_intrusive<Condition>(condition_id, pFace, pProperties);

	    pSkinCondition->Set(ACTIVE,false);

	    //set to beam tube mesh
	    mrModelPart.AddCondition(pSkinCondition);

	    counter++;

	  }
	  else if( counter == number_of_angles ){

	    FaceNodesIds[0] = wall_nodes_number_id + i ;
	    FaceNodesIds[1] = wall_nodes_number_id + i + 1 - number_of_angles;
	    FaceNodesIds[2] = wall_nodes_number_id + i + 1;
	    FaceNodesIds[3] = wall_nodes_number_id + i + number_of_angles;

	    GeometryType::PointsArrayType FaceNodes;
	    FaceNodes.reserve(4);

	    for(unsigned int j=0; j<4; j++)
	      FaceNodes.push_back(mrModelPart.pGetNode(FaceNodesIds[j]));

	    pFace = Kratos::make_shared<Quadrilateral3DType>(FaceNodes);

	    pSkinCondition = Kratos::make_intrusive<Condition>(condition_id, pFace, pProperties);

	    pSkinCondition->Set(ACTIVE,false);

	    //set to beam tube mesh
	    mrModelPart.AddCondition(pSkinCondition);

	    counter = 1;
	  }


	}

      KRATOS_CATCH( "" )
    }

    //*******************************************************************************************
    //*******************************************************************************************

    static inline unsigned int GetMaxNodeId(ModelPart& rModelPart)
    {
      KRATOS_TRY

      unsigned int max_id = rModelPart.Nodes().back().Id();

      for(ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin(); i_node!= rModelPart.NodesEnd(); i_node++)
	{
	  if(i_node->Id() > max_id)
	    max_id = i_node->Id();
	}

      return max_id;

      KRATOS_CATCH( "" )
    }

    //*******************************************************************************************
    //*******************************************************************************************

    static inline unsigned int GetMaxElementId(ModelPart& rModelPart)
    {
      KRATOS_TRY

      if( rModelPart.NumberOfElements() == 0 )
	return 0;

      unsigned int max_id = rModelPart.Elements().back().Id();

      for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin(); i_elem!= rModelPart.ElementsEnd(); i_elem++)
	{
	  if(i_elem->Id() > max_id)
	    max_id = i_elem->Id();
	}

      return max_id;

      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    static inline unsigned int GetMaxConditionId(ModelPart& rModelPart)
    {
      KRATOS_TRY

      if( rModelPart.NumberOfConditions() == 0 )
	return 0;

      unsigned int max_id = rModelPart.Conditions().back().Id();

      for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); i_cond++)
	{
	  if(i_cond->Id() > max_id)
	    max_id = i_cond->Id();
	}

      return max_id;

      KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    NodeType::Pointer CreateNode (ModelPart& rModelPart, PointType& rPoint, const unsigned int& rNodeId)
    {

      KRATOS_TRY

      NodeType::Pointer Node = rModelPart.CreateNewNode( rNodeId, rPoint[0], rPoint[1], rPoint[2]);

      //generating the dofs
      NodeType::DofsContainerType& reference_dofs = (rModelPart.NodesBegin())->GetDofs();


      for(NodeType::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
      	{
      	  NodeType::DofType& rDof = **iii;
      	  Node->pAddDof( rDof );
      	}

      //set fix dofs:
      NodeType::DofsContainerType& new_dofs = Node->GetDofs();

      for(NodeType::DofsContainerType::iterator iii = new_dofs.begin(); iii != new_dofs.end(); iii++)
      	{
      	  NodeType::DofType& rDof = **iii;
	  rDof.FixDof(); // dofs free
      	}

      //generating step data
      /*unsigned int buffer_size = (rModelPart.NodesBegin())->GetBufferSize();
      unsigned int step_data_size = rModelPart.GetNodalSolutionStepDataSize();
      for(unsigned int step = 0; step<buffer_size; step++)
      	{
      	  double* NodeData = Node->SolutionStepData().Data(step);
      	  double* ReferenceData = (rModelPart.NodesBegin())->SolutionStepData().Data(step);

      	  //copying this data in the position of the vector we are interested in
      	  for(unsigned int j= 0; j<step_data_size; j++)
      	    {
      	      NodeData[j] = ReferenceData[j];
      	    }
      	}
      */

      return Node;

      KRATOS_CATCH( "" )

    }

    //************************************************************************************
    //************************************************************************************

    void CleanElementNeighbours()
    {

      KRATOS_TRY

      NodesContainerType&    rNodes = mrModelPart.Nodes();
      ElementsContainerType& rElements = mrModelPart.Elements();

      //first of all the neighbour nodes and neighbour elements arrays are initialized to the guessed size
      //this cleans the old entries:

      unsigned int AverageNodes = 2;
      unsigned int AverageElements = 2;

      //*************  Erase old node neighbours  *************//
      for(auto& i_node : rNodes)
      {
        NodeWeakPtrVectorType& nNodes = i_node.GetValue(NEIGHBOUR_NODES);
        nNodes.clear();
        nNodes.resize(AverageNodes);

        ElementWeakPtrVectorType& nElements = i_node.GetValue(NEIGHBOUR_ELEMENTS);
        nElements.clear();
        nElements.resize(AverageElements);
      }

      //************* Erase old element neighbours ************//
      for(auto& i_elem : rElements)
        {
	  Element::GeometryType& rGeometry = i_elem.GetGeometry();
	  int size= rGeometry.FacesNumber();

	  ElementWeakPtrVectorType& nElements = i_elem.GetValue(NEIGHBOUR_ELEMENTS);
          nElements.clear();
	  nElements.resize(size);

        }

      KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************
    template<class TDataType> void  AddUniquePointer
    (GlobalPointersVector<TDataType>& v, const typename TDataType::WeakPointer candidate)
    {
      typename GlobalPointersVector< TDataType >::iterator i = v.begin();
      typename GlobalPointersVector< TDataType >::iterator endit = v.end();
      while ( i != endit && (i)->Id() != (candidate)->Id())
      {
        i++;
      }
      if( i == endit )
      {
        v.push_back(candidate);
      }

    }

    //************************************************************************************
    //************************************************************************************

    ElementWeakPtrType CheckForNeighbourElems1D (unsigned int Id_1, ElementWeakPtrVectorType& nElements, ElementsContainerType::iterator i_elem)
    {
      //look for the faces around node Id_1
      for(auto i_nelem(nElements.begin()); i_nelem != nElements.end(); ++i_nelem)
      {
        //look for the nodes of the neighbour faces
        Geometry<Node >& nGeometry = i_nelem->GetGeometry();
        if(nGeometry.LocalSpaceDimension() == 1){
          for(unsigned int node_i = 0; node_i < nGeometry.size(); ++node_i)
          {
            if(nGeometry[node_i].Id() == Id_1)
            {
              if(i_nelem->Id() != i_elem->Id())
              {
                return *i_nelem.base();
              }
            }
          }
        }
      }
      return *i_elem.base();
    }

    //************************************************************************************
    //************************************************************************************

    void SearchNeighbours()
    {
      KRATOS_TRY

      NodesContainerType& rNodes = mrModelPart.Nodes();
      ElementsContainerType& rElements = mrModelPart.Elements();

      //first of all the neighbour nodes and neighbour elements arrays are initialized to the guessed size
      //this cleans the old entries:

      //*****  Erase old node and element neighbours  *********//
      CleanElementNeighbours();


      //*************  Neigbours of nodes  ************//
      //add the neighbour elements to all the nodes in the mesh
      for(auto i_elem(rElements.begin()); i_elem != rElements.end(); ++i_elem)
      {
        Element::GeometryType& rGeometry = i_elem->GetGeometry();
        for(unsigned int i = 0; i < rGeometry.size(); i++)
        {
          rGeometry[i].GetValue(NEIGHBOUR_ELEMENTS).push_back( *i_elem.base() );
        }
      }

      //adding the neighbouring nodes to all nodes in the mesh
      for(auto& i_node : rNodes)
      {
        ElementWeakPtrVectorType& nElements = i_node.GetValue(NEIGHBOUR_ELEMENTS);
        for(auto& i_nelem : nElements)
        {
          Element::GeometryType& rGeometry = i_nelem.GetGeometry();
          for(unsigned int i = 0; i < rGeometry.size(); i++)
          {
            if( rGeometry[i].Id() != i_node.Id() )
            {
              NodeWeakPtrVectorType& nNodes = i_nelem.GetValue(NEIGHBOUR_NODES);
              AddUniquePointer< Node >(nNodes, rGeometry(i));
            }

          }
        }
      }


      //*************  Neigbours of elements  *********//
      //add the neighbour elements to all the elements in the mesh

      //loop over faces
      for(auto i_elem(rElements.begin()); i_elem != rElements.end(); ++i_elem)
      {
        //face nodes
        Geometry<Node >& rGeometry = i_elem->GetGeometry();
        if( rGeometry.FacesNumber() == 2 ){

          ElementWeakPtrVectorType& nElements = i_elem->GetValue(NEIGHBOUR_ELEMENTS);

          //vector of the 2 faces around the given face
          if( nElements.size() != 2 )
            nElements.resize(2);

          //neighb_face is the vector containing pointers to the three faces around ic:

          unsigned int size = rGeometry.size();

          // neighbour element over edge 0 of element ic;
          nElements(0) = CheckForNeighbourElems1D(rGeometry[0].Id(), rGeometry[0].GetValue(NEIGHBOUR_ELEMENTS), i_elem);
          // neighbour element over edge 1 of element ic;
          nElements(1) = CheckForNeighbourElems1D(rGeometry[size-1].Id(), rGeometry[size-1].GetValue(NEIGHBOUR_ELEMENTS), i_elem);

        }
      }


      KRATOS_CATCH( "" )
    }


    void TransferSkinToOutput()
    {
        KRATOS_TRY

	ModelPart& rOutputModelPart = GetOutputModelPart();

	std::vector<std::size_t> NodeIds;
	// set nodes to ACTIVE (to write them in output)
	for (ModelPart::NodeIterator i = mrModelPart.NodesBegin(); i != mrModelPart.NodesEnd(); ++i)
	{
           NodeIds.push_back(i->Id());
	}
	std::vector<std::size_t> ConditionIds;
	// set elements ACTIVE (to write them in output)
	for (ModelPart::ConditionIterator i = mrModelPart.ConditionsBegin(); i != mrModelPart.ConditionsEnd(); ++i)
	{
          ConditionIds.push_back(i->Id());
	}

	rOutputModelPart.AddNodes(NodeIds);
	rOutputModelPart.AddConditions(ConditionIds);

	SetInActiveFlag(rOutputModelPart,TO_ERASE);

        KRATOS_CATCH("")
    }

    ModelPart& GetOutputModelPart()
    {
        KRATOS_TRY

	std::string OutputModelPartName;
	ModelPart& rMainModelPart = mrModelPart.GetParentModelPart();
	for(ModelPart::SubModelPartIterator i_mp= rMainModelPart.SubModelPartsBegin(); i_mp!=rMainModelPart.SubModelPartsEnd(); i_mp++)
	{
	  if( i_mp->Is(ACTIVE) )
	    OutputModelPartName = i_mp->Name();
	}

	return (rMainModelPart.GetSubModelPart(OutputModelPartName));

        KRATOS_CATCH("")
    }


    virtual void SetActiveFlag(ModelPart& rModelPart, Flags id_flag = ACTIVE)
    {
        KRATOS_TRY

	// set nodes to ACTIVE (to write them in output)
	for (ModelPart::NodeIterator i = rModelPart.NodesBegin(); i != rModelPart.NodesEnd(); ++i)
	{
	  (i)->Set(id_flag,true);
	}

	// set elements ACTIVE (to write them in output)
	for (ModelPart::ConditionIterator i = rModelPart.ConditionsBegin(); i != rModelPart.ConditionsEnd(); ++i)
	{
	  (i)->Set(id_flag,true);
	}


        KRATOS_CATCH("")
    }

    virtual void SetInActiveFlag(ModelPart& rModelPart, Flags id_flag = ACTIVE)
    {
        KRATOS_TRY


	// set nodes to NOT ACTIVE (to not consider them in computation)
	for (ModelPart::NodeIterator i = rModelPart.NodesBegin(); i != rModelPart.NodesEnd(); ++i)
	{
	  (i)->Set(id_flag,false);
	}

	// set elements NOT ACTIVE (to not consider them in computation)
	for (ModelPart::ConditionIterator i = rModelPart.ConditionsBegin(); i != rModelPart.ConditionsEnd(); ++i)
	{
	  (i)->Set(id_flag,false);
	}


        KRATOS_CATCH("")
    }

    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    BuildStringSkinProcess& operator=(BuildStringSkinProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class BuildStringSkinProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  BuildStringSkinProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BuildStringSkinProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_BUILD_STRING_SKIN_PROCESS_H_INCLUDED  defined
