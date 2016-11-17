//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:          July 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_SPATIAL_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_SPATIAL_BOUNDING_BOX_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>
#include <limits>

#include "includes/kratos_flags.h"
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

class KRATOS_API(PFEM_BASE_APPLICATION) SpatialBoundingBox
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
      
      mpGapNormal=NULL;
      mpGapTangent=NULL;

      mRadius      = 0.0;
      mContactFace = 0;

    };

    /**
     * Constructor with Node and Contact Parameters
     */
    BoundingBoxParameters(const NodeType& rNode, double& rGapNormal, double& rGapTangent, PointType& rNormal, PointType& rTangent)
    {
      mpPoint        = &(rNode.Coordinates());
      mpVelocity     = &(rNode.FastGetSolutionStepValue(VELOCITY));

      mpCurrentDisplacement  = &(rNode.FastGetSolutionStepValue(DISPLACEMENT,0));
      mpPreviousDisplacement = &(rNode.FastGetSolutionStepValue(DISPLACEMENT,1));
      
      mpNormal       = &rNormal;
      mpTangent      = &rTangent;
      
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
      Dimension = 2;
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
      InitialLocalQuaternion  = QuaternionType::FromRotationMatrix( InitialLocalMatrix );
      LocalQuaternion         = QuaternionType::FromRotationMatrix( InitialLocalMatrix );

    }

    void SetInitialValues()
    {
      InitialUpperPoint = Center;
      InitialLowerPoint  = Center;
      InitialCenter    = Center;
    }

    void UpdatePosition( PointType& Displacement )
    {
      UpperPoint = InitialUpperPoint + Displacement;
      LowerPoint = InitialLowerPoint  + Displacement;
      Center     = InitialCenter    + Displacement;
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
      mBox.LowerPoint  = rLowerPoint;

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
      mBox.LowerPoint  = mBox.Center - Side;

      //set to local frame
      this->MapToLocalFrame(mBox.InitialLocalQuaternion,mBox);

      mBox.SetInitialValues();

      mRigidBodyCenterSupplied = false;

      KRATOS_CATCH("")
    }
 
    //**************************************************************************
    //**************************************************************************

    SpatialBoundingBox(ModelPart &rModelPart, const double& rRadius)
    {
      KRATOS_TRY

      double max=std::numeric_limits<double>::max();
      double min=std::numeric_limits<double>::min();
      
      ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();
      const unsigned int dimension = element_begin->GetGeometry().WorkingSpaceDimension();

      PointType Maximum;
      PointType Minimum;

      for(unsigned int i=0; i<3; i++)
	{
	  Maximum[i] = min; 
	  Minimum[i] = max; 
	}

      //Get inside point of the subdomains

      for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin(); in!=rModelPart.NodesEnd(); in++)
	{
	  if(in->Is(BOUNDARY) ){
	  
	    //get maximum
	    if(Maximum[0]<in->X())
	      Maximum[0]=in->X();

	    if(Maximum[1]<in->Y())
	      Maximum[1]=in->Y();

	    if(Maximum[3]<in->Z())
	      Maximum[3]=in->Z();

	    //get minimum
	    if(Minimum[0]>in->X())
	      Minimum[0]=in->X();

	    if(Minimum[1]>in->Y())
	      Minimum[1]=in->Y();
	  
	    if(Minimum[3]>in->Z())
	      Minimum[3]=in->Z();	 
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
      Side[0] = mBox.Radius;
      Side[1] = mBox.Radius;
      Side[2] = mBox.Radius;

      mBox.UpperPoint = mBox.Center + Side;
      mBox.LowerPoint  = mBox.Center - Side;

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


      if(norm_2((mBox.Center-LocalPoint)) > 2 * mBox.Radius)
	inside = false;

      for(unsigned int i=0; i<mBox.Center.size(); i++)
	{
	  if(mBox.UpperPoint[i]<LocalPoint[i]){
	    inside = false;
	    break;
	  }
	}

      for(unsigned int i=0; i<mBox.Center.size(); i++)
	{
	  if(mBox.LowerPoint[i]>LocalPoint[i]){
	    inside = false;
	    break;
	  }
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
      for(unsigned int i=0; i<mBox.Center.size(); i++)
	{
	  if( mBox.UpperPoint[i] == mBox.LowerPoint[i] )
	    {
	      count++;
	    }
	}

      if( count == mBox.Center.size() )
	return true;
      //check if the box is not set


      PointType LocalPoint = rPoint;
      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, LocalPoint);


      if(norm_2((mBox.Center-LocalPoint)) > 2 * mBox.Radius)
	inside = false;

      for(unsigned int i=0; i<mBox.Center.size(); i++)
	{
	  if(mBox.UpperPoint[i]<LocalPoint[i]){
	    inside = false;
	    break;
	  }
	}

      for(unsigned int i=0; i<mBox.Center.size(); i++)
	{
	  if(mBox.LowerPoint[i]>LocalPoint[i]){
	    inside = false;
	    break;
	  }
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

   

    ///@}
    ///@name Access
    ///@{


    // SET

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
	for(unsigned int i=0; i<3; i++)
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
       PointType Center = mBox.Center;
       BeamMathUtilsType::MapToReferenceLocalFrame(mBox.InitialLocalQuaternion, Center);
       return Center;
    }

    //**************************************************************************
    //**************************************************************************

    virtual PointType GetCenter(const PointType& rPoint)
    {
       PointType Center = mBox.Center;
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

      unsigned int start=0;
      unsigned int NumberOfMeshes=rModelPart.NumberOfMeshes();
      if(NumberOfMeshes>1) 
	start=1;

      std::vector<PointType > Holes;
      for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	{
	  PointType Point(dimension);
	  for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin(MeshId); in!=rModelPart.NodesEnd(MeshId); in++)
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

	}
	

      return Holes;
    }

    //**************************************************************************
    //**************************************************************************

    /// Compute vertices
    std::vector<PointType > GetVertices(const double& rCurrentTime, const unsigned int& rDimension)
    {
    
      std::vector<PointType> vertices;

      PointType Displacement = this->GetBoxDisplacement( rCurrentTime );

      PointType Reference = mBox.UpperPoint + Displacement;
      
      PointType Side = mBox.UpperPoint - mBox.LowerPoint;

      //point 1
      vertices.push_back(Reference);
      
      Reference[0] -= Side[0];
      
      //point 2
      vertices.push_back(Reference);
      
      Reference[1] -= Side[1];
      
      //point 3
      vertices.push_back(Reference);
      
      Reference[0] += Side[0];
      
      //point 4
      vertices.push_back(Reference);
      

      if( rDimension == 3 ){ 

	Reference = mBox.LowerPoint + Displacement;
	
	//point 5
	vertices.push_back(Reference);
      
	Reference[0] += Side[0];
      
	//point 6
	vertices.push_back(Reference);
      
	Reference[1] += Side[1];
      
	//point 7
	vertices.push_back(Reference);
      
	Reference[0] -= Side[0];
      
	//point 8
	vertices.push_back(Reference);

      }

      return vertices;
           
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
	for( int i=0; i<3; i++ )
	  Displacement[i] = CurrentDisplacement[i];

	if( mpRigidBodyCenter->SolutionStepsDataHas(ROTATION) ){
	  array_1d<double, 3 > & CurrentRotation = mpRigidBodyCenter->FastGetSolutionStepValue(ROTATION);
	  for( int i=0; i<3; i++ )
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
	for( int i=0; i<3; i++ )
	  Displacement[i] = CurrentDisplacement[i]-PreviousDisplacement[i];

	// TODO: treatment of rotation in displacement calculation
	// if( mpRigidBodyCenter->SolutionStepsDataHas(ROTATION) ){
	//   array_1d<double, 3 > & CurrentRotation = mpRigidBodyCenter->FastGetSolutionStepValue(ROTATION,0);
	//   array_1d<double, 3 > & PreviousRotation = mpRigidBodyCenter->FastGetSolutionStepValue(ROTATION,1);
	//   for( int i=0; i<3; i++ )
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
           
      rTangent = ZeroVector(3);
   
      //1.-compute contact tangent (following relative movement)
      PointType PointDeltaDisplacement = rValues.GetDeltaDisplacement();

      rTangent = PointDeltaDisplacement - inner_prod(PointDeltaDisplacement, rNormal) * rNormal;

      PointType UnitTangent = rTangent;
      if( norm_2(UnitTangent) )
	UnitTangent /= norm_2(UnitTangent);
      
      //std::cout<<" rTangent "<<rTangent<<" PointDelta "<<PointDeltaDisplacement<<std::endl;

      
      //2.-compute tangent direction
      PointType BoxDeltaDisplacement = this->GetBoxDeltaDisplacement(rCurrentProcessInfo[TIME], rCurrentProcessInfo.GetPreviousTimeStepInfo()[TIME]);


      if( norm_2(UnitTangent) != 0 ){

	rTangent = rTangent - inner_prod(BoxDeltaDisplacement, UnitTangent) * UnitTangent;
	
      }
      else{

	rTangent = BoxDeltaDisplacement - inner_prod(BoxDeltaDisplacement, rNormal) * rNormal;
	
      }
      
      //std::cout<<" rTangent "<<rTangent<<" BoxDelta "<<BoxDeltaDisplacement<<" PointDelta "<<PointDeltaDisplacement<<std::endl;
      
      //3.-compute  normal gap
      rGapTangent = norm_2(rTangent);

      
      if(norm_2(rTangent))
	rTangent/= norm_2(rTangent);
	
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


