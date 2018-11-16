//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_SPATIAL_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_SPATIAL_BOUNDING_BOX_H_INCLUDED


// System includes
#include <limits>

#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/beam_math_utilities.hpp"

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

class KRATOS_API(DELAUNAY_MESHING_APPLICATION) SpatialBoundingBox
{
public:

  //typedef bounded_vector<double, 3>                     PointType;
  typedef array_1d<double, 3>                             PointType;
  typedef ModelPart::NodeType                              NodeType;
  typedef ModelPart::NodesContainerType          NodesContainerType;
  typedef NodesContainerType::Pointer     NodesContainerTypePointer;
  typedef BeamMathUtils<double>                   BeamMathUtilsType;
  typedef Quaternion<double>                         QuaternionType;


  struct BoundingBoxParameters
  {
    KRATOS_CLASS_POINTER_DEFINITION(BoundingBoxParameters);

  private:

    const PointType*  mpPoint;
    const PointType*  mpVelocity;

    const PointType*  mpCurrentDisplacement;
    const PointType*  mpPreviousDisplacement;

    PointType*        mpNormal;
    PointType*        mpTangent;
    PointType*        mpRelativeDisplacement;

    double*           mpGapNormal;
    double*           mpGapTangent;

    double            mRadius;
    int               mContactFace;

  public:

    /**
     * Constructor.
     */
    BoundingBoxParameters ()
    {
      //Initialize pointers to NULL
      mpPoint=NULL;
      mpVelocity=NULL;

      mpCurrentDisplacement=NULL;
      mpPreviousDisplacement=NULL;

      mpNormal=NULL;
      mpTangent=NULL;
      mpRelativeDisplacement=NULL;

      mpGapNormal=NULL;
      mpGapTangent=NULL;

      mRadius      = 0.0;
      mContactFace = 0;

    };

    /**
     * Constructor with Node
     */
    BoundingBoxParameters(const NodeType& rNode)
    {

      mpPoint        = &(rNode.Coordinates());
      mpVelocity     = &(rNode.FastGetSolutionStepValue(VELOCITY));

      mpCurrentDisplacement  = &(rNode.FastGetSolutionStepValue(DISPLACEMENT,0));
      mpPreviousDisplacement = &(rNode.FastGetSolutionStepValue(DISPLACEMENT,1));

      mpNormal=NULL;
      mpTangent=NULL;

      mpRelativeDisplacement=NULL;

      mpGapNormal=NULL;
      mpGapTangent=NULL;

      mRadius      = 0.0;
      mContactFace = 0;

    };

    /**
     * Constructor with Node and Contact Parameters
     */
    BoundingBoxParameters(const NodeType& rNode, double& rGapNormal, double& rGapTangent, PointType& rNormal, PointType& rTangent, PointType& rDisplacement)
    {
      mpPoint        = &(rNode.Coordinates());
      mpVelocity     = &(rNode.FastGetSolutionStepValue(VELOCITY));

      mpCurrentDisplacement  = &(rNode.FastGetSolutionStepValue(DISPLACEMENT,0));
      mpPreviousDisplacement = &(rNode.FastGetSolutionStepValue(DISPLACEMENT,1));

      mpNormal       = &rNormal;
      mpTangent      = &rTangent;

      mpRelativeDisplacement = &rDisplacement;

      mpGapNormal    = &rGapNormal;
      mpGapTangent   = &rGapTangent;

      mRadius        = 0.0;
      mContactFace   = 0;
    };

    //set parameters
    void SetNode(const NodeType& rNode){
      mpPoint        = &(rNode.Coordinates());
      mpVelocity     = &(rNode.FastGetSolutionStepValue(VELOCITY));

      mpCurrentDisplacement  = &(rNode.FastGetSolutionStepValue(DISPLACEMENT,0));
      mpPreviousDisplacement = &(rNode.FastGetSolutionStepValue(DISPLACEMENT,1));
    };

    void SetPoint(const PointType& rPoint) {mpPoint = &rPoint;};
    void SetVelocity(const PointType& rVelocity) {mpVelocity = &rVelocity;};
    void SetCurrentDisplacement(const PointType& rDisplacement) {mpCurrentDisplacement = &rDisplacement;};
    void SetPreviousDisplacement(const PointType& rDisplacement) {mpPreviousDisplacement = &rDisplacement;};

    void SetNormal(PointType& rNormal)   {mpNormal = &rNormal;};
    void SetTangent(PointType& rTangent) {mpTangent = &rTangent;};
    void SetRelativeDisplacement(PointType& rDisplacement) {mpRelativeDisplacement = &rDisplacement;};

    void SetGapNormal(double& rGapNormal)   {mpGapNormal = &rGapNormal;};
    void SetGapTangent(double& rGapTangent) {mpGapTangent = &rGapTangent;};

    void SetRadius(double Radius)          {mRadius = Radius;};
    void SetContactFace(int ContactFace)   {mContactFace = ContactFace;};

    //get parameters
    const PointType& GetPoint()        {return *mpPoint;};
    const PointType& GetVelocity()     {return *mpVelocity;};
    const PointType& GetCurrentDisplacement() {return *mpCurrentDisplacement;};
    const PointType& GetPreviousDisplacement() {return *mpPreviousDisplacement;};

    PointType GetDeltaDisplacement() {return ((*mpCurrentDisplacement)-(*mpPreviousDisplacement));};

    PointType& GetNormal()  {return *mpNormal;};
    PointType& GetTangent() {return *mpTangent;};
    PointType& GetRelativeDisplacement() {return *mpRelativeDisplacement;};


    double& GetGapNormal()  {return *mpGapNormal;};
    double& GetGapTangent() {return *mpGapTangent;};

    double& GetRadius()     {return mRadius;};
    int& GetContactFace()   {return mContactFace;};

  };// struct BoundingBoxParameters end


protected:


  typedef struct
  {

    int        Dimension;          // 2D or 3D
    bool       Axisymmetric;       // true or false
    int        Convexity;          // 1 or -1  if "in" is inside or outside respectively
    double     Radius;             // box radius

    PointType  InitialUpperPoint;  // box highest point
    PointType  InitialLowerPoint;  // box lowest point
    PointType  InitialCenter;      // center current position

    PointType  UpperPoint;         // box highest point
    PointType  LowerPoint;         // box lowest point
    PointType  Center;             // center current position

    PointType  Velocity;           // box velocity
    PointType  AngularVelocity;    // box rotation

    QuaternionType  InitialLocalQuaternion; //initial local axes for the box
    QuaternionType  LocalQuaternion;        //local axes for the box

  public:

    void Initialize()
    {
      Dimension = 3;
      Axisymmetric = false;
      Convexity = 1;
      Radius = 0;

      UpperPoint.clear();
      LowerPoint.clear();
      Center.clear();

      InitialUpperPoint.clear();
      InitialLowerPoint.clear();
      InitialCenter.clear();

      Velocity.clear();
      AngularVelocity.clear();

      Matrix InitialLocalMatrix = IdentityMatrix(3);
      InitialLocalQuaternion = QuaternionType::FromRotationMatrix( InitialLocalMatrix );
      LocalQuaternion        = QuaternionType::FromRotationMatrix( InitialLocalMatrix );

    }

    void SetInitialValues()
    {
      InitialUpperPoint = UpperPoint;
      InitialLowerPoint = LowerPoint;
      InitialCenter     = Center;
    }

    void UpdatePosition( PointType& Displacement )
    {
      UpperPoint = InitialUpperPoint + Displacement;
      LowerPoint = InitialLowerPoint + Displacement;
      Center     = InitialCenter     + Displacement;
    }

    void Print()
    {
      std::cout<<" [--SPATIAL-BOX--] "<<std::endl;
      std::cout<<"  [Center:"<<Center<<std::endl;
      std::cout<<"  [UpperPoint:"<<UpperPoint<<std::endl;
      std::cout<<"  [LowerPoint:"<<LowerPoint<<std::endl;
      std::cout<<" [--------] "<<std::endl;
    }


  } BoundingBoxVariables;



public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( SpatialBoundingBox );



    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpatialBoundingBox()
    {
      KRATOS_TRY

      mBox.Initialize();
      mRigidBodyCenterSupplied = false;
      //std::cout<< " Calling Bounding Box empty constructor" <<std::endl;

      KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************

    SpatialBoundingBox(Parameters CustomParameters)
    {

      KRATOS_TRY

      Parameters DefaultParameters( R"(
            {
                "parameters_list":[{
                    "upper_point": [0.0, 0.0, 0.0],
                    "lower_point": [0.0, 0.0, 0.0],
                    "convexity": 1
                 }],
                 "velocity": [0.0, 0.0, 0.0],
                 "angular_velocity": [0.0, 0.0, 0.0],
                 "local_axes":[ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
            }  )" );


      //validate against defaults -- this also ensures no type mismatch
      CustomParameters.ValidateAndAssignDefaults(DefaultParameters);

      if(CustomParameters["parameters_list"].IsArray() == true && CustomParameters["parameters_list"].size() != 1)
        {
	  KRATOS_THROW_ERROR(std::runtime_error,"paramters_list for the Spatial BBX must contain only one term",CustomParameters.PrettyPrintJsonString());
        }

      mBox.Initialize();

      Parameters BoxParameters = CustomParameters["parameters_list"][0];

      mBox.UpperPoint[0] = BoxParameters["upper_point"][0].GetDouble();
      mBox.UpperPoint[1] = BoxParameters["upper_point"][1].GetDouble();
      mBox.UpperPoint[2] = BoxParameters["upper_point"][2].GetDouble();

      mBox.LowerPoint[0] = BoxParameters["lower_point"][0].GetDouble();
      mBox.LowerPoint[1] = BoxParameters["lower_point"][1].GetDouble();
      mBox.LowerPoint[2] = BoxParameters["lower_point"][2].GetDouble();

      mBox.Center = 0.5 * ( mBox.UpperPoint + mBox.LowerPoint );
      mBox.Radius = 0.5 * norm_2(mBox.UpperPoint - mBox.LowerPoint);

      mBox.Velocity[0] = CustomParameters["velocity"][0].GetDouble();
      mBox.Velocity[1] = CustomParameters["velocity"][1].GetDouble();
      mBox.Velocity[2] = CustomParameters["velocity"][2].GetDouble();

      mBox.AngularVelocity[0] = CustomParameters["angular_velocity"][0].GetDouble();
      mBox.AngularVelocity[1] = CustomParameters["angular_velocity"][1].GetDouble();
      mBox.AngularVelocity[2] = CustomParameters["angular_velocity"][2].GetDouble();

      mBox.Convexity = BoxParameters["convexity"].GetInt();

      Matrix InitialLocalMatrix = IdentityMatrix(3);

      unsigned int size = CustomParameters["local_axes"].size();

      for( unsigned int i=0; i<size; i++ )
	{
	  Parameters LocalAxesRow = CustomParameters["local_axes"][i];

	  InitialLocalMatrix(0,i) = LocalAxesRow[0].GetDouble(); //column disposition
	  InitialLocalMatrix(1,i) = LocalAxesRow[1].GetDouble();
	  InitialLocalMatrix(2,i) = LocalAxesRow[2].GetDouble();
	}

      //set to local frame
      this->MapToLocalFrame(mBox.InitialLocalQuaternion,mBox);

      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, mBox.Velocity);
      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, mBox.AngularVelocity);

      mBox.SetInitialValues();

      mRigidBodyCenterSupplied = false;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    SpatialBoundingBox(const PointType& rLowerPoint, const PointType& rUpperPoint )
    {
      KRATOS_TRY

      mBox.Initialize();
      mBox.UpperPoint = rUpperPoint;
      mBox.LowerPoint = rLowerPoint;

      mBox.Center = 0.5 * ( rUpperPoint + rLowerPoint );
      mBox.Radius = 0.5 * norm_2(rUpperPoint-rLowerPoint);

      //set to local frame
      this->MapToLocalFrame(mBox.InitialLocalQuaternion,mBox);

      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, mBox.Velocity);
      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, mBox.AngularVelocity);

      mBox.SetInitialValues();

      mRigidBodyCenterSupplied = false;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    SpatialBoundingBox(const PointType& rCenter, const double& rRadius)
    {
      KRATOS_TRY

      mBox.Initialize();
      mBox.Center = rCenter;
      mBox.Radius = rRadius;

      PointType Side;
      Side[0] = mBox.Radius;
      Side[1] = mBox.Radius;
      Side[2] = mBox.Radius;

      mBox.UpperPoint = mBox.Center + Side;
      mBox.LowerPoint = mBox.Center - Side;

      //set to local frame
      this->MapToLocalFrame(mBox.InitialLocalQuaternion,mBox);

      mBox.SetInitialValues();

      mRigidBodyCenterSupplied = false;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    SpatialBoundingBox(ModelPart &rModelPart, const double& rRadius, double factor = 0)
    {
      KRATOS_TRY

      double max=std::numeric_limits<double>::max();
      double min=std::numeric_limits<double>::min();


      unsigned int dimension = 2;
      if ( rModelPart.NumberOfElements() > 0)
         dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
      else if ( rModelPart.NumberOfConditions() > 0)
         dimension = rModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
      else
         KRATOS_ERROR << " spatial_bounding_box: the supplied ModelPart does not have elements or conditions " << std::endl;

      PointType Maximum;
      PointType Minimum;

      for(unsigned int i=0; i<3; ++i)
	{
	  Maximum[i] = min;
	  Minimum[i] = max;
	}

      //Get inside point of the subdomains

      for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin(); in!=rModelPart.NodesEnd(); ++in)
	{
	  if(in->Is(BOUNDARY) ){

	    //get maximum
	    if(Maximum[0]<in->X())
	      Maximum[0]=in->X();

	    if(Maximum[1]<in->Y())
	      Maximum[1]=in->Y();

	    if(Maximum[2]<in->Z())
	      Maximum[2]=in->Z();

	    //get minimum
	    if(Minimum[0]>in->X())
	      Minimum[0]=in->X();

	    if(Minimum[1]>in->Y())
	      Minimum[1]=in->Y();

	    if(Minimum[2]>in->Z())
	      Minimum[2]=in->Z();
	  }

	}

      mBox.Initialize();

      mBox.Center = 0.5*(Maximum+Minimum);

      double MaxRadius = min;

      if(Maximum[0]-Minimum[0] > MaxRadius)
	MaxRadius = Maximum[0]-Minimum[0];

      if(Maximum[1]-Minimum[1] > MaxRadius)
	MaxRadius = Maximum[1]-Minimum[1];

      if(Maximum[2]-Minimum[2]>MaxRadius)
	MaxRadius = Maximum[2]-Minimum[2];


      mBox.Radius = rRadius + 0.5*(MaxRadius);

      PointType Side(dimension);
      Side[0] = mBox.Radius + mBox.Radius * factor;
      Side[1] = mBox.Radius + mBox.Radius * factor;
      Side[2] = mBox.Radius + mBox.Radius * factor;

      mBox.UpperPoint = mBox.Center + Side;
      mBox.LowerPoint = mBox.Center - Side;

      //set to local frame
      this->MapToLocalFrame(mBox.InitialLocalQuaternion,mBox);

      mRigidBodyCenterSupplied = false;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************


    /// Assignment operator.
    virtual SpatialBoundingBox& operator=(SpatialBoundingBox const& rOther)
    {
      KRATOS_TRY

      mpRigidBodyCenter = rOther.mpRigidBodyCenter;
      mRigidBodyCenterSupplied = rOther.mRigidBodyCenterSupplied;
      mBox = rOther.mBox;

      return *this;

      KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************

    /// Copy constructor.
    SpatialBoundingBox(SpatialBoundingBox const& rOther)
      :mpRigidBodyCenter(rOther.mpRigidBodyCenter)
      ,mRigidBodyCenterSupplied(rOther.mRigidBodyCenterSupplied)
      ,mBox(rOther.mBox)
    {
    }

    //**************************************************************************
    //**************************************************************************


    /// Destructor.
    virtual ~SpatialBoundingBox() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //**************************************************************************
    //**************************************************************************

    virtual void UpdateBoxPosition(const double & rCurrentTime)
    {

      KRATOS_TRY

      PointType Displacement  =  this->GetBoxDisplacement(rCurrentTime);

      mBox.UpdatePosition(Displacement);

      this->MapToLocalFrame(mBox.LocalQuaternion, mBox);

      KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************


    virtual bool IsInside (const PointType& rPoint, double& rCurrentTime, double Radius = 0)
    {

      KRATOS_TRY

      bool inside = true;

      PointType Displacement  =  this->GetBoxDisplacement(rCurrentTime);

      mBox.UpdatePosition(Displacement);

      this->MapToLocalFrame(mBox.LocalQuaternion, mBox);

      PointType LocalPoint = rPoint;
      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, LocalPoint);

      if(    (mBox.UpperPoint[0]>=LocalPoint[0] && mBox.LowerPoint[0]<=LocalPoint[0])
	  && (mBox.UpperPoint[1]>=LocalPoint[1] && mBox.LowerPoint[1]<=LocalPoint[1])
	  && (mBox.UpperPoint[2]>=LocalPoint[2] && mBox.LowerPoint[2]<=LocalPoint[2]) ){
	inside = true;
      }
      else{
	inside = false;
      }

      QuaternionType LocaQuaternionlConjugate = mBox.LocalQuaternion.conjugate();
      this->MapToLocalFrame(LocaQuaternionlConjugate, mBox);

      return inside;

      KRATOS_CATCH("")

    }


    //**************************************************************************
    //**************************************************************************


    virtual bool IsInside (const PointType& rPoint)
    {

      KRATOS_TRY

      bool inside = true;

      //check if the box is not set
      unsigned int count = 0;
      for(unsigned int i=0; i<mBox.Center.size(); ++i)
	{
	  if( mBox.UpperPoint[i] == mBox.LowerPoint[i] )
	    {
	      count++;
	    }
	}

      if( count == mBox.Center.size() ){
	std::cout<<" IsInside:: warning spatial BOX not set "<<std::endl;
	return true;
      }
      //check if the box is not set


      PointType LocalPoint = rPoint;
      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, LocalPoint);

      // std::cout<<" Local Point "<<LocalPoint<<std::endl;
      // std::cout<<" Upper "<<mBox.UpperPoint<<" Lower "<<mBox.LowerPoint<<std::endl;
      // if(!(mBox.UpperPoint[0]>=LocalPoint[0] && mBox.LowerPoint[0]<=LocalPoint[0]) )
      // 	std::cout<<" first not fit "<<std::endl;
      // if(!(mBox.UpperPoint[1]>=LocalPoint[1] && mBox.LowerPoint[1]<=LocalPoint[1]) )
      // 	std::cout<<" second not fit "<<std::endl;
      // if(!(mBox.UpperPoint[2]>=LocalPoint[2] && mBox.LowerPoint[2]<=LocalPoint[2]) )
      // 	std::cout<<" third not fit "<<std::endl;


      if(    (mBox.UpperPoint[0]>=LocalPoint[0] && mBox.LowerPoint[0]<=LocalPoint[0])
	  && (mBox.UpperPoint[1]>=LocalPoint[1] && mBox.LowerPoint[1]<=LocalPoint[1])
	  && (mBox.UpperPoint[2]>=LocalPoint[2] && mBox.LowerPoint[2]<=LocalPoint[2]) ){
	inside = true;
      }
      else{
	inside = false;
      }


      return inside;

      KRATOS_CATCH("")

    }



    //************************************************************************************
    //************************************************************************************
    virtual bool IsInside(BoundingBoxParameters& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      std::cout<< "Calling empty method" <<std::endl;
      return false;

      KRATOS_CATCH("")
    }


    virtual void GetParametricDirections(BoundingBoxParameters & rValues, Vector & rT1, Vector & rT2)
    {
        KRATOS_TRY

        //std::cout<< "Calling empty method directions" <<std::endl;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Access
    ///@{


    // SET

    void SetUpperPoint(PointType& rUpperPoint)
    {
      mBox.UpperPoint        = rUpperPoint;
      mBox.InitialUpperPoint = rUpperPoint;
    }

    //**************************************************************************
    //**************************************************************************

    void SetLowerPoint(PointType& rLowerPoint)
    {
      mBox.LowerPoint        = rLowerPoint;
      mBox.InitialLowerPoint = rLowerPoint;

    }

    //**************************************************************************
    //**************************************************************************

    void SetVelocity(PointType& rVelocity)
    {
      mBox.Velocity = rVelocity;
      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, mBox.Velocity);
    }

    //**************************************************************************
    //**************************************************************************

    void SetAngularVelocity(PointType& rAngularVelocity)
    {
      mBox.AngularVelocity = rAngularVelocity;
      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, mBox.AngularVelocity);
    }

    //**************************************************************************
    //**************************************************************************

    void SetDimension(int dimension)
    {
        mBox.Dimension = dimension;
    }

    //**************************************************************************
    //**************************************************************************

    void SetAxisymmetric()
    {
        mBox.Axisymmetric = true;
    }

    //**************************************************************************
    //**************************************************************************

    void SetRigidBodyCenter(NodeType::Pointer pCenter)
    {
      mpRigidBodyCenter = pCenter;
      mRigidBodyCenterSupplied = true;
    }

    // GET

    //**************************************************************************
    //**************************************************************************

    virtual double GetRadius()
    {
      return mBox.Radius;
    }

    //**************************************************************************
    //**************************************************************************

    virtual double GetRadius(const PointType& rPoint)
    {
      return mBox.Radius;
    }

    //**************************************************************************
    //**************************************************************************

    virtual PointType GetVelocity()
    {
      PointType Velocity;
      if( mRigidBodyCenterSupplied ){
	array_1d<double, 3>& rVelocity = mpRigidBodyCenter->FastGetSolutionStepValue(VELOCITY);
	for(unsigned int i=0; i<3; ++i)
	  Velocity[i] = rVelocity[i];
      }
      else{
	Velocity = mBox.Velocity;
	BeamMathUtilsType::MapToReferenceLocalFrame(mBox.InitialLocalQuaternion, Velocity);
      }
      return Velocity;
    }

    //**************************************************************************
    //**************************************************************************

    virtual PointType GetCenter()
    {
       BeamMathUtilsType::MapToReferenceLocalFrame(mBox.InitialLocalQuaternion, mBox.Center);
       return mBox.Center;
    }

    //**************************************************************************
    //**************************************************************************

    virtual PointType GetCenter(const PointType& rPoint)
    {
      KRATOS_WARNING("") << "Calling spatial bounding box base class method "<<std::endl;
      BeamMathUtilsType::MapToReferenceLocalFrame(mBox.InitialLocalQuaternion, mBox.Center);
      return mBox.Center;
    }

    //**************************************************************************
    //**************************************************************************

    /// Compute inside holes
    std::vector<PointType > GetHoles(ModelPart &rModelPart)
    {
      //Get inside point of the subdomains
      ModelPart::ConditionsContainerType::iterator condition_begin = rModelPart.ConditionsBegin();
      const unsigned int dimension = condition_begin->GetGeometry().WorkingSpaceDimension();

      std::vector<PointType > Holes;
      PointType Point(dimension);
      for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin(); in!=rModelPart.NodesEnd(); ++in)
	{
	  if(in->IsNot(BOUNDARY) ){
	    Point[0] = in->X();
	    Point[1] = in->Y();

	    if(dimension>2)
	      Point[2] = in->Z();

	    Holes.push_back(Point);
	    break;
	  }
	  else{
	    array_1d<double, 3>& Normal = in->FastGetSolutionStepValue(NORMAL);

	    std::cout<<" Normal "<<Normal<<std::endl;
	    double tolerance = 0.25;

	    Point[0] = in->X() - Normal[0] * tolerance;
	    Point[1] = in->Y() - Normal[1] * tolerance;
	    if(dimension>2)
	      Point[2] = in->Z() - Normal[2] * tolerance;

	    Holes.push_back(Point);
	    break;
	  }
	}


      return Holes;
    }

    //**************************************************************************
    //**************************************************************************

    /// Compute vertices
    void GetVertices(std::vector<PointType >& rVertices, const double& rCurrentTime, const unsigned int& rDimension)
    {

      PointType Displacement = this->GetBoxDisplacement( rCurrentTime );

      PointType Reference = mBox.UpperPoint + Displacement;

      PointType Side = mBox.UpperPoint - mBox.LowerPoint;

      Reference[1] -= Side[1];

      //point 0
      rVertices.push_back(Reference);

      Reference[1] += Side[1];

      //point 1
      rVertices.push_back(Reference);

      Reference[0] -= Side[0];

      //point 2
      rVertices.push_back(Reference);

      Reference[1] -= Side[1];

      //point 3
      rVertices.push_back(Reference);


      if( rDimension == 3 ){

	Reference = mBox.LowerPoint + Displacement;

	Reference[0] += Side[0];

	//point 4
	rVertices.push_back(Reference);

	Reference[0] -= Side[0];

	//point 5
	rVertices.push_back(Reference);

	Reference[1] += Side[1];

	//point 6
	rVertices.push_back(Reference);

	Reference[0] += Side[0];

	//point 7
	rVertices.push_back(Reference);

      }

    }

    //************************************************************************************
    //************************************************************************************

    void GetTriangularFaces(DenseMatrix<unsigned int>& rFaces, const unsigned int& rDimension)
    {
      KRATOS_TRY

      if( rDimension == 2 ){

	if(rFaces.size1() != 4 || rFaces.size2() != 2)
	  rFaces.resize(4,2,false);

	rFaces(0,0)=0;
	rFaces(0,1)=1;

	rFaces(1,0)=1;
	rFaces(1,1)=2;

	rFaces(2,0)=2;
	rFaces(2,1)=3;

	rFaces(3,0)=3;
	rFaces(3,1)=4;

      }
      else if ( rDimension == 3 ){

	if(rFaces.size1() != 12 || rFaces.size2() != 3)
	  rFaces.resize(12,3,false);

	rFaces(0,0)=0;
	rFaces(0,1)=1;
	rFaces(0,2)=3;

	rFaces(1,0)=3;
	rFaces(1,1)=1;
	rFaces(1,2)=2;

	rFaces(2,0)=3;
	rFaces(2,1)=2;
	rFaces(2,2)=6;

	rFaces(3,0)=6;
	rFaces(3,1)=5;
	rFaces(3,2)=3;

	rFaces(4,0)=5;
	rFaces(4,1)=6;
	rFaces(4,2)=7;

	rFaces(5,0)=7;
	rFaces(5,1)=4;
	rFaces(5,2)=5;

	rFaces(6,0)=0;
	rFaces(6,1)=4;
	rFaces(6,2)=7;

	rFaces(7,0)=7;
	rFaces(7,1)=1;
	rFaces(7,2)=0;

	rFaces(8,0)=0;
	rFaces(8,1)=3;
	rFaces(8,2)=5;

	rFaces(9,0)=5;
	rFaces(9,1)=4;
	rFaces(9,2)=0;

	rFaces(10,0)=1;
	rFaces(10,1)=7;
	rFaces(10,2)=6;

	rFaces(11,0)=6;
	rFaces(11,1)=2;
	rFaces(11,2)=1;

      }

      KRATOS_CATCH("")
    }


    void GetQuadrilateralFaces(DenseMatrix<unsigned int>& rFaces, const unsigned int& rDimension)
    {
      KRATOS_TRY

      if( rDimension == 2 ){

	if(rFaces.size1() != 4 || rFaces.size2() != 2)
	  rFaces.resize(4,2,false);

	rFaces(0,0)=0;
	rFaces(0,1)=1;

	rFaces(1,0)=1;
	rFaces(1,1)=2;

	rFaces(2,0)=2;
	rFaces(2,1)=3;

	rFaces(3,0)=3;
	rFaces(3,1)=4;

      }
      else if ( rDimension == 3 ){

	if(rFaces.size1() != 6 || rFaces.size2() != 4)
	  rFaces.resize(6,4,false);

	rFaces(0,0)=0;
	rFaces(0,1)=1;
	rFaces(0,2)=2;
	rFaces(0,3)=3;

	rFaces(1,0)=3;
	rFaces(1,1)=2;
	rFaces(1,2)=6;
	rFaces(1,3)=5;

	rFaces(2,0)=5;
	rFaces(2,1)=6;
	rFaces(2,2)=7;
	rFaces(2,3)=4;

	rFaces(3,0)=4;
	rFaces(3,1)=7;
	rFaces(3,2)=1;
	rFaces(3,3)=0;

	rFaces(4,0)=0;
	rFaces(4,1)=3;
	rFaces(4,2)=5;
	rFaces(4,3)=4;

	rFaces(5,0)=1;
	rFaces(5,1)=7;
	rFaces(5,2)=6;
	rFaces(5,3)=2;

      }

      KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    virtual void CreateBoundingBoxBoundaryMesh(ModelPart& rModelPart, int linear_partitions = 4, int angular_partitions = 4 )
    {
      KRATOS_TRY

	std::cout<< " Calling Spatial Bounding Box empty boundary mesh creation" <<std::endl;

      KRATOS_CATCH("")
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "SpatialBoundingBox";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << mBox.UpperPoint << " , " << mBox.LowerPoint;
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

    NodeType::Pointer   mpRigidBodyCenter;

    bool   mRigidBodyCenterSupplied;

    BoundingBoxVariables mBox;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    //**************************************************************************
    //**************************************************************************

    void MapToLocalFrame(QuaternionType& rQuaternion, BoundingBoxVariables& rBox)
    {
      KRATOS_TRY

      BeamMathUtilsType::MapToCurrentLocalFrame(rQuaternion, rBox.UpperPoint);
      BeamMathUtilsType::MapToCurrentLocalFrame(rQuaternion, rBox.LowerPoint);
      BeamMathUtilsType::MapToCurrentLocalFrame(rQuaternion, rBox.Center);

      KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************

    PointType GetBoxDisplacement(const double& rCurrentTime)
    {

      PointType Displacement = ZeroVector(3);
      PointType Rotation = ZeroVector(3);

      if( mRigidBodyCenterSupplied ){

	array_1d<double, 3 > & CurrentDisplacement = mpRigidBodyCenter->FastGetSolutionStepValue(DISPLACEMENT);
	for( int i=0; i<3; ++i )
	  Displacement[i] = CurrentDisplacement[i];

	if( mpRigidBodyCenter->SolutionStepsDataHas(ROTATION) ){
	  array_1d<double, 3 > & CurrentRotation = mpRigidBodyCenter->FastGetSolutionStepValue(ROTATION);
	  for( int i=0; i<3; ++i )
	    Rotation[i] = CurrentRotation[i];
	}

	//local base rotation
	BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, Rotation);
	mBox.LocalQuaternion = QuaternionType::FromRotationVector(Rotation);

	BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, Displacement);

      }
      else{

	Displacement = mBox.Velocity * rCurrentTime;
	Rotation     = mBox.AngularVelocity * rCurrentTime;

	//local base rotation
	BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, Rotation);
	mBox.LocalQuaternion = QuaternionType::FromRotationVector(Rotation);

      }

      return Displacement;
    }


    //**************************************************************************
    //**************************************************************************

    PointType GetBoxDeltaDisplacement(const double& rCurrentTime, const double& rPreviousTime)
    {

      //Quaternion LocalQuaternion = mBox.LocalQuaternion;

      PointType Displacement = ZeroVector(3);
      //PointType Rotation = ZeroVector(3);

      if( mRigidBodyCenterSupplied ){

	array_1d<double, 3 > & CurrentDisplacement = mpRigidBodyCenter->FastGetSolutionStepValue(DISPLACEMENT,0);
	array_1d<double, 3 > & PreviousDisplacement = mpRigidBodyCenter->FastGetSolutionStepValue(DISPLACEMENT,1);
	for( int i=0; i<3; ++i )
	  Displacement[i] = CurrentDisplacement[i]-PreviousDisplacement[i];

	// TODO: treatment of rotation in displacement calculation
	// if( mpRigidBodyCenter->SolutionStepsDataHas(ROTATION) ){
	//   array_1d<double, 3 > & CurrentRotation = mpRigidBodyCenter->FastGetSolutionStepValue(ROTATION,0);
	//   array_1d<double, 3 > & PreviousRotation = mpRigidBodyCenter->FastGetSolutionStepValue(ROTATION,1);
	//   for( int i=0; i<3; ++i )
	//     Rotation[i] = CurrentRotation[i]-PreviousRotation[i];
	// }

	// //local base rotation
	// BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, Rotation);
	// mBox.LocalQuaternion = QuaternionType::FromRotationVector(Rotation);

	BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, Displacement);

      }
      else{

	Displacement = mBox.Velocity * (rCurrentTime - rPreviousTime);

	// TODO: treatment of rotation in displacement calculation
	// Rotation     = mBox.AngularVelocity * (rCurrentTime - rPreviousTime);
	// //local base rotation
	// BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, Rotation);
	// mBox.LocalQuaternion = QuaternionType::FromRotationVector(Rotation);

      }

      //mBox.LocalQuaternion = LocalQuaternion;

      return Displacement;
    }

    //**************************************************************************
    //**************************************************************************

    void ComputeContactTangent(BoundingBoxParameters& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      PointType& rNormal      = rValues.GetNormal();
      PointType& rTangent     = rValues.GetTangent();
      double& rGapTangent     = rValues.GetGapTangent();

      PointType& rRelativeDisplacement = rValues.GetRelativeDisplacement();

      noalias(rTangent) = ZeroVector(3);

      //1.-compute contact tangent (following relative movement)
      PointType PointDeltaDisplacement = rValues.GetDeltaDisplacement();

      //2.-compute tangent direction
      PointType BoxDeltaDisplacement = this->GetBoxDeltaDisplacement(rCurrentProcessInfo[TIME], rCurrentProcessInfo.GetPreviousTimeStepInfo()[TIME]);

      rRelativeDisplacement = BoxDeltaDisplacement-PointDeltaDisplacement;

      rTangent = (rRelativeDisplacement) - inner_prod(rRelativeDisplacement, rNormal) * rNormal;

      if( !norm_2(rNormal) )
	noalias(rTangent) = ZeroVector(3);

      //3.-compute  normal gap
      rGapTangent = norm_2(rTangent);

      if(norm_2(rTangent))
	rTangent/= norm_2(rTangent);

      //std::cout<<" Normal "<<rNormal<<" Tangent "<<rTangent<<" gapT "<<rGapTangent<<std::endl;

      KRATOS_CATCH( "" )

    }

    //*******************************************************************************************
    //*******************************************************************************************

    static inline unsigned int GetMaxNodeId(ModelPart& rModelPart)
    {
      KRATOS_TRY

      unsigned int max_id = rModelPart.Nodes().back().Id();

      for(ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin(); i_node!= rModelPart.NodesEnd(); ++i_node)
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

      for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin(); i_elem!= rModelPart.ElementsEnd(); ++i_elem)
	{
	  if(i_elem->Id() > max_id)
	    max_id = i_elem->Id();
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


      for(NodeType::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); ++iii)
      	{
      	  NodeType::DofType& rDof = *iii;
      	  Node->pAddDof( rDof );
      	}

      //set fix dofs:
      NodeType::DofsContainerType& new_dofs = Node->GetDofs();

      for(NodeType::DofsContainerType::iterator iii = new_dofs.begin(); iii != new_dofs.end(); ++iii)
      	{
      	  NodeType::DofType& rDof = *iii;
	  rDof.FixDof(); // dofs free
      	}

      //generating step data:
      // unsigned int buffer_size = (rModelPart.NodesBegin())->GetBufferSize();
      // unsigned int step_data_size = rModelPart.GetNodalSolutionStepDataSize();
      // for(unsigned int step = 0; step<buffer_size; ++step)
      // 	{
      // 	  double* NodeData = Node->SolutionStepData().Data(step);
      // 	  double* ReferenceData = (rModelPart.NodesBegin())->SolutionStepData().Data(step);

      // 	  //copying this data in the position of the vector we are interested in
      // 	  for(unsigned int j= 0; j<step_data_size; ++j)
      // 	    {
      // 	      NodeData[j] = ReferenceData[j];
      // 	    }
      // 	}

      return Node;

      KRATOS_CATCH( "" )

    }


    //************************************************************************************
    //************************************************************************************

    void CalculateOrthonormalBase(PointType & rDirectionVectorX, PointType & rDirectionVectorY, PointType & rDirectionVectorZ)
    {
      KRATOS_TRY

      PointType GlobalY = ZeroVector(3);
      GlobalY[1]=1.0;

      PointType GlobalZ = ZeroVector(3);
      GlobalZ[2]=1.0;


      // local x-axis (e1_local) is the beam axis  (in GID is e3_local)
      double VectorNorm = MathUtils<double>::Norm(rDirectionVectorX);
      if( VectorNorm != 0)
	rDirectionVectorX /= VectorNorm;

      // local y-axis (e2_local) (in GID is e1_local)
      rDirectionVectorY = ZeroVector(3);

      double tolerance = 1.0/64.0;
      if(fabs(rDirectionVectorX[0])< tolerance && fabs(rDirectionVectorX[1])< tolerance){
	MathUtils<double>::CrossProduct(rDirectionVectorY, GlobalY, rDirectionVectorX);
      }
      else{
	MathUtils<double>::CrossProduct(rDirectionVectorY, GlobalZ, rDirectionVectorX);
      }


      VectorNorm = MathUtils<double>::Norm(rDirectionVectorY);
      if( VectorNorm != 0)
	rDirectionVectorY /= VectorNorm;

      // local z-axis (e3_local) (in GID is e2_local)
      MathUtils<double>::CrossProduct(rDirectionVectorZ, rDirectionVectorX,rDirectionVectorY);

      VectorNorm = MathUtils<double>::Norm(rDirectionVectorZ);
      if( VectorNorm != 0 )
	rDirectionVectorZ /= VectorNorm;

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


}; // Class SpatialBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_SPATIAL_BOUNDING_BOX_H_INCLUDED  defined
