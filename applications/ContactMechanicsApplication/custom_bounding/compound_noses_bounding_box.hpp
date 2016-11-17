//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_COMPOUND_NOSES_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_COMPOUND_NOSES_BOUNDING_BOX_H_INCLUDED

// External includes

// System includes

// Project includes
#include "custom_bounding/spatial_bounding_box.hpp"


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

    This Box represents a 2D wall composed by a set of two-line sistems called "Noses"
    each pair of lines are tangent to the semi-circle given by a center and a radius
    these semi-circles represent the "noses" or the "tips" of the box wall
    One line is the upper line of the nose defined by a rake angle (respect to the vertical axis)
    The other line is a down line defined by a clearance angle (respect to the horizontal axis)
    
    A convexity parameter is given to determine which side of each nose is considered 
    the internal or external boundary

    This bounding box is essentially used for rigid wall contact purposes
*/

class KRATOS_API(CONTACT_MECHANICS_APPLICATION) CompoundNosesBoundingBox
  : public SpatialBoundingBox
{
public:

  //typedef bounded_vector<double, 3>                     PointType;
  typedef array_1d<double, 3>                             PointType;
  typedef ModelPart::NodeType                              NodeType;
  typedef ModelPart::NodesContainerType          NodesContainerType;
  typedef NodesContainerType::Pointer     NodesContainerTypePointer;
  typedef BeamMathUtils<double>                   BeamMathUtilsType;
  typedef Quaternion<double>                         QuaternionType;

private:

   enum ContactFace{ FreeSurface=0, RakeSurface=1, TipSurface=2, ClearanceSurface=3 };

protected:

   typedef struct
   {
     int    Convexity;      //1 or -1 if "in" is inside or outside respectively
    
     double Radius;         //nose radius
     double RakeAngle;      //top angle,    from vertical axis       --> #RakeAngle      = (90-RakeAngle)
     double ClearanceAngle; //bottom angle, from the horizontal axis --> #ClearanceAngle = (180-ClearanceAngle)
      
     double TangentRakeAngle;      //tan((pi/2)-#RakeAngle)
     double TangentClearanceAngle; //tan((pi/2)-#ClearanceAngle)
    
     PointType  InitialCenter; // center original position
     PointType  Center;        // center current position
    
    
   public:
    
     void Initialize()
     {
       Convexity = 1;
       Radius = 0;
       RakeAngle = 0;
       ClearanceAngle = 0;
      
       TangentRakeAngle = 0;
       TangentClearanceAngle = 0;

       InitialCenter.clear();
       Center.clear();
     }
     
     void SetInitialValues()
     {
       InitialCenter = Center;
     }

     void UpdatePosition( PointType& Displacement )
     {
       Center = InitialCenter + Displacement;
     }

     
   } BoxNoseVariables;


public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CompoundNosesBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( CompoundNosesBoundingBox );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CompoundNosesBoundingBox() : SpatialBoundingBox()
    {
      KRATOS_TRY

      std::cout<< "Calling Rigid Nose Wall BBX empty constructor" <<std::endl;

      KRATOS_CATCH("")     
    }

    //**************************************************************************
    //**************************************************************************

    CompoundNosesBoundingBox(Parameters CustomParameters)
    {

      KRATOS_TRY

      Parameters DefaultParameters( R"(
            {
                "parameters_list": [],
                "velocity": [0.0, 0.0, 0.0],
                "angular_velocity": [0.0, 0.0, 0.0],
                "upper_point": [0.0, 0.0, 0.0],
                "lower_point": [0.0, 0.0, 0.0],
                "local_axes":[ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
            }  )" );


      //items in the parameters_list:
      //{"radius": 0.0,
      // "center": [0.0, 0.0, 0.0],
      // "rake_angle": 0.0,
      // "clearance_angle": 0.0,
      // "convexity": 0}

      //validate against defaults -- this also ensures no type mismatch
      CustomParameters.ValidateAndAssignDefaults(DefaultParameters);       

      unsigned int list_size = CustomParameters["parameters_list"].size();
      
      Vector Radius         (list_size);
      Matrix Centers        (list_size,3);
      Vector RakeAngles     (list_size);
      Vector ClearanceAngles(list_size);
      Vector Convexities    (list_size);

      for( unsigned int i=0; i<list_size; i++ )
	{
	  Parameters NoseParameters = CustomParameters["parameters_list"][i];
	  
	  Radius[i]          = NoseParameters["radius"].GetDouble();

	  Centers(i,0)       = NoseParameters["center"][0].GetDouble();
	  Centers(i,1)       = NoseParameters["center"][1].GetDouble();
	  Centers(i,2)       = NoseParameters["center"][2].GetDouble();

	  RakeAngles[i]      = NoseParameters["rake_angle"].GetDouble();
	  ClearanceAngles[i] = NoseParameters["clearance_angle"].GetDouble();
	  Convexities[i]     = NoseParameters["convexity"].GetInt();
	}

      Vector Velocity(3);
      Velocity[0] = CustomParameters["velocity"][0].GetDouble();
      Velocity[1] = CustomParameters["velocity"][1].GetDouble();
      Velocity[2] = CustomParameters["velocity"][2].GetDouble();

      Vector AngularVelocity(3);
      AngularVelocity[0] = CustomParameters["angular_velocity"][0].GetDouble();
      AngularVelocity[1] = CustomParameters["angular_velocity"][1].GetDouble();
      AngularVelocity[2] = CustomParameters["angular_velocity"][2].GetDouble();

      Vector UpperPoint(3);
      UpperPoint[0] = CustomParameters["upper_point"][0].GetDouble();
      UpperPoint[1] = CustomParameters["upper_point"][1].GetDouble();
      UpperPoint[2] = CustomParameters["upper_point"][2].GetDouble();

      Vector LowerPoint(3);
      LowerPoint[0] = CustomParameters["lower_point"][0].GetDouble();
      LowerPoint[1] = CustomParameters["lower_point"][1].GetDouble();
      LowerPoint[2] = CustomParameters["lower_point"][2].GetDouble();

      Matrix InitialLocalMatrix = IdentityMatrix(3);

      unsigned int size = CustomParameters["local_axes"].size();

      for( unsigned int i=0; i<size; i++ )
	{
	  Parameters LocalAxesRow = CustomParameters["local_axes"][i];

	  InitialLocalMatrix(0,i) = LocalAxesRow[0].GetDouble(); //column disposition
	  InitialLocalMatrix(1,i) = LocalAxesRow[1].GetDouble();
	  InitialLocalMatrix(2,i) = LocalAxesRow[2].GetDouble();
	} 

      this->CreateCompoundNosesBox(Radius,RakeAngles,ClearanceAngles,Centers,Convexities,
				   Velocity, AngularVelocity, UpperPoint, LowerPoint, InitialLocalMatrix);

      KRATOS_CATCH("")

    }

  
    //**************************************************************************
    //**************************************************************************

    // General Wall constructor
    CompoundNosesBoundingBox(Vector Radius,
			     Vector RakeAngles,
			     Vector ClearanceAngles,
			     Matrix Centers,
			     Vector Convexities,
			     Vector Velocity,
			     Vector AngularVelocity,
			     Vector UpperPoint,
			     Vector LowerPoint,
			     Matrix InitialLocalMatrix) //column disposition of the local_axes matrix
    {
    
      KRATOS_TRY

      this->CreateCompoundNosesBox(Radius,RakeAngles,ClearanceAngles,Centers,Convexities,
				   Velocity, AngularVelocity, UpperPoint, LowerPoint, InitialLocalMatrix);

      KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************

    /// Assignment operator.
    CompoundNosesBoundingBox& operator=(CompoundNosesBoundingBox const& rOther)
    {
      SpatialBoundingBox::operator=(rOther);
      mBoxNoses = rOther.mBoxNoses;
      return *this;
    }

    //**************************************************************************
    //**************************************************************************

    /// Copy constructor.
    CompoundNosesBoundingBox(CompoundNosesBoundingBox const& rOther) 
    :SpatialBoundingBox(rOther)
    ,mBoxNoses(rOther.mBoxNoses)
    {
    }

    //**************************************************************************
    //**************************************************************************

    /// Destructor.
    virtual ~CompoundNosesBoundingBox() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
  
    //**************************************************************************
    //**************************************************************************

    double GetRadius(const PointType& rPoint)
    {
      PointType LocalPoint = rPoint;
      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, LocalPoint);
      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.LocalQuaternion, LocalPoint);

      LocalPoint[2] = 0; //2D
	
      return GetNearerRadius(LocalPoint);
    }

    //**************************************************************************
    //**************************************************************************

    PointType GetCenter(const PointType& rPoint)
    {
      PointType LocalPoint = rPoint;
      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, LocalPoint);
      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.LocalQuaternion, LocalPoint);

      LocalPoint[2] = 0; //2D
      
      PointType LocalNearerPoint = GetNearerCenter(LocalPoint);
      BeamMathUtilsType::MapToReferenceLocalFrame(mBox.InitialLocalQuaternion, LocalNearerPoint);
      BeamMathUtilsType::MapToReferenceLocalFrame(mBox.LocalQuaternion, LocalNearerPoint);

      return  LocalNearerPoint;
    }


    //**************************************************************************
    //**************************************************************************

    void UpdateBoxPosition(const double & rCurrentTime)
    {
      
      PointType Displacement  =  this->GetBoxDisplacement(rCurrentTime);

      mBox.UpdatePosition(Displacement);

      this->MapToLocalFrame(mBox.LocalQuaternion, mBox);

      for(unsigned int i=0; i<mBoxNoses.size(); i++)
	{
	  mBoxNoses[i].Center =  mBoxNoses[i].InitialCenter + Displacement;
	  BeamMathUtilsType::MapToCurrentLocalFrame(mBox.LocalQuaternion, mBoxNoses[i].Center);

	  std::cout<<"   BOX NOSE center position "<<mBoxNoses[i].Center<<std::endl;
	}
    }

    //************************************************************************************
    //************************************************************************************
   

    bool IsInside (const PointType& rPoint, double& rCurrentTime, double Radius = 0)
    {
      
      KRATOS_TRY   

      bool is_inside = SpatialBoundingBox::IsInside(rPoint);

      if( is_inside ){

	//std::cout<<" Local Point A "<<rPoint<<std::endl;

	PointType LocalPoint = rPoint;
	BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, LocalPoint);
	BeamMathUtilsType::MapToCurrentLocalFrame(mBox.LocalQuaternion, LocalPoint);  

	LocalPoint[2] = 0; //2D

	//std::cout<<" Local Point B "<<rPoint<<std::endl;

	unsigned int SelectedNose = BoxNoseSearch(LocalPoint);

      	BoxNoseVariables& rWallNose = mBoxNoses[SelectedNose];

	double NoseRadius = rWallNose.Radius;

	if( rWallNose.Convexity == 1)
	  NoseRadius *= 2; //increase the bounding box 

	if( rWallNose.Convexity == -1)
	  NoseRadius *= 0.1; //decrease the bounding box 

	// bool node_in =false;
	// if( rWallNose.Convexity == 1 && (LocalPoint[0]>=95 && LocalPoint[1]>=6.9) )
	//  	node_in = true;

	switch( ContactSearch(LocalPoint, NoseRadius, rWallNose) )
	  {
	  case FreeSurface:      
	    is_inside = false;
	    // ContactFace = 0;
	    // if( node_in )
	    // std::cout<<" Nose "<<SelectedNose<<" [ FreeSurface :"<<rPoint<<"]"<<std::endl;
	    break;
	  case RakeSurface:      
	    is_inside = true;
	    // ContactFace = 1;
	    // if( node_in )
	    // std::cout<<" Nose "<<SelectedNose<<" [ RakeSurface :"<<rPoint<<"]"<<std::endl;
	    break;
	  case TipSurface:       
	    is_inside = true;
	    // ContactFace = 2;
	    // if( node_in )
	    // std::cout<<" Nose "<<SelectedNose<<" [ TipSurface :"<<rPoint<<"]"<<std::endl;
	    break;
	  case ClearanceSurface: 
	    is_inside = true;
	    // ContactFace = 3;
	    // if( node_in )
	    // std::cout<<" Nose "<<SelectedNose<<" [ ClearanceSurface :"<<rPoint<<"]"<<std::endl;
	    break;
	  default:               
	    is_inside = false;
	    break;
	  }
     
      }

      return is_inside;
      
     KRATOS_CATCH("")
 
    }


    //************************************************************************************
    //************************************************************************************
    
    bool IsInside(BoundingBoxParameters& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      const PointType& rPoint = rValues.GetPoint();
      PointType& rNormal      = rValues.GetNormal();     
      double& rGapNormal      = rValues.GetGapNormal();
      
      bool is_inside = SpatialBoundingBox::IsInside(rPoint);

      if( is_inside ){

	PointType LocalPoint = rPoint;
	BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, LocalPoint);
	BeamMathUtilsType::MapToCurrentLocalFrame(mBox.LocalQuaternion, LocalPoint);

	LocalPoint[2] = 0; //2D

	unsigned int SelectedNose = BoxNoseSearch(LocalPoint);
      
	BoxNoseVariables& rWallNose = mBoxNoses[SelectedNose];

	double NoseRadius = rWallNose.Radius;

	//std::cout<<" Convexity ["<<SelectedNose<<" ]: "<< rWallNose.Convexity <<std::endl;

	switch( ContactSearch(LocalPoint, NoseRadius, rWallNose) )
	  {	  
	  case FreeSurface:      
	    is_inside = false;
	    rValues.SetContactFace(0);
	    break;
	  case RakeSurface:      
	    is_inside = CalculateRakeSurface(LocalPoint, NoseRadius, rGapNormal, rNormal, rWallNose);
	    rValues.SetContactFace(1);
	    break;
	  case TipSurface:       
	    is_inside = CalculateTipSurface(LocalPoint, NoseRadius, rGapNormal, rNormal, rWallNose);
	    rValues.SetContactFace(2);
	    break;
	  case ClearanceSurface: 
	    is_inside = CalculateClearanceSurface(LocalPoint, NoseRadius, rGapNormal, rNormal, rWallNose);
	    rValues.SetContactFace(3);
	    break;
	  default:               
	    is_inside = false;
	    break;
	  }

	// if(rWallNose.Convexity == -1  && ContactFace == 3 && rGapNormal<1.0){
	//  	std::cout<<" [ ContactFace: "<<ContactFace<<"; Normal: "<<rNormal<<"; GapNormal: "<< rGapNormal <<" ] "<<rPoint<<std::endl;
	// }
	
	BeamMathUtilsType::MapToReferenceLocalFrame(mBox.InitialLocalQuaternion, rNormal);
	BeamMathUtilsType::MapToReferenceLocalFrame(mBox.LocalQuaternion, rNormal);

	this->ComputeContactTangent(rValues,rCurrentProcessInfo);
      }

      //std::cout<<" ["<<is_inside<<"][ ContactFace: "<<ContactFace<<"; Normal: "<<rNormal<<"; GapNormal: "<< rGapNormal <<" ] "<<rPoint<<std::endl;
      return is_inside;
      
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
        return "CompoundNosesBoundingBox";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    std::vector<BoxNoseVariables> mBoxNoses;

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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    //************************************************************************************
    //************************************************************************************

    PointType GetNearerCenter(const PointType& rPoint)
    {

      unsigned int SelectedNose = BoxNoseSearch(rPoint);
 
      BoxNoseVariables& rWallNose = mBoxNoses[SelectedNose];

      return rWallNose.Center;
    }

    //************************************************************************************
    //************************************************************************************

    double GetNearerRadius(const PointType& rPoint)
    {
      unsigned int SelectedNose = BoxNoseSearch(rPoint);
 
      BoxNoseVariables& rWallNose = mBoxNoses[SelectedNose];

      return rWallNose.Radius;

    }

    //************************************************************************************
    //************************************************************************************

    bool CheckValidAngle( double& rAngle )
    {
      bool valid_angle = true;
      double tol = 1e-3;
      for( int i = 0; i<=360; i+=90 ){
	if( fabs(rAngle) < double(i + tol) && fabs(rAngle) > double(i - tol) ){
	  valid_angle = false;
	  break;
	}
      }
      
      return valid_angle;
    }

    //************************************************************************************
    //************************************************************************************

    // General Wall creator 
    void CreateCompoundNosesBox(Vector& rRadius, Vector& rRakeAngles, Vector& rClearanceAngles, Matrix& rCenters, Vector& rConvexities, 
				Vector& rVelocity, Vector& rAngularVelocity, Vector& rUpperPoint, Vector& rLowerPoint, Matrix& rInitialLocalMatrix)
    {
    
      KRATOS_TRY

      if( rRadius.size() != rRakeAngles.size() || rRakeAngles.size() != rClearanceAngles.size() )
	std::cout<<" Introduced walls are not consistent in sizes "<<std::endl;
      
      double pi = 3.141592654;
      
      std::cout<<" [--NOSE-WALL--] "<<std::endl;
      std::cout<<"  [NOSES:"<<rRadius.size()<<"]"<<std::endl;

      mBox.Initialize();     
      
      for(unsigned int i=0; i<rRadius.size(); i++)
	{
	  BoxNoseVariables WallNose;
	  WallNose.Initialize();

	  WallNose.Convexity      = rConvexities[i];
	  WallNose.Radius         = rRadius[i];

	  // RakeAngle :: Angle given respect to the vertical axis  (is positive line, represents the nose upper part) 
	  // changed to be expressed respect to the horitzontal axis, represents positive (increasing) line
	  WallNose.RakeAngle      = ( 90 - rRakeAngles[i] );        

	  // ClearanceAngle :: Angle given respect to the vertical axis  (is negative line, represents the nose down part) 
	  // changed to represent a negative (decreasing) line
	  WallNose.ClearanceAngle = ( 180 + rClearanceAngles[i] ); 

	  bool valid_angle = CheckValidAngle( WallNose.RakeAngle );

	  //check if the angle is 0 or 180 before performing this operation
	  if( valid_angle ){
	    WallNose.RakeAngle *= pi / 180.0;
	    WallNose.TangentRakeAngle = tan(0.5*pi-WallNose.RakeAngle);
	  }
	  else{
	    WallNose.RakeAngle *= pi / 180.0;
	    WallNose.TangentRakeAngle = 0;
	  }

	  valid_angle = CheckValidAngle( WallNose.ClearanceAngle );

	  //check if the angle is 0 or 180 before performing this operation
	  if( valid_angle ){
	    WallNose.ClearanceAngle *= pi / 180.0;
	    WallNose.TangentClearanceAngle = tan(0.5*pi-WallNose.ClearanceAngle);
	  }
	  else{
	    WallNose.ClearanceAngle *= pi / 180.0;
	    WallNose.TangentClearanceAngle = 0;
	  }
	  
	  for(unsigned int j=0; j<rCenters.size2(); j++)
	    {
	      WallNose.Center[j] = rCenters(i,j);
	    }
	  
	  WallNose.SetInitialValues();

	  //set to local frame
	  BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, WallNose.InitialCenter);
	  BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, WallNose.Center);

	  std::cout<<"  [COMPONENT]["<<i<<"]"<<std::endl;
	  std::cout<<"  [Convexity:"<<WallNose.Convexity<<std::endl;
	  std::cout<<"  [Radius:"<<WallNose.Radius<<std::endl;
	  std::cout<<"  [Center:"<<WallNose.Center<<std::endl;
	  std::cout<<"  [Rake:"<<WallNose.RakeAngle<<std::endl;
	  std::cout<<"  [Clearance:"<<WallNose.ClearanceAngle<<std::endl;
	  std::cout<<"  [TangentRakeAngle:"<<WallNose.TangentRakeAngle<<std::endl;
	  std::cout<<"  [TangentClearanceAngle:"<<WallNose.TangentClearanceAngle<<std::endl;
	  
	  mBoxNoses.push_back(WallNose);

	}

      std::cout<<" [--------] "<<std::endl;
      
      mBox.InitialLocalQuaternion = QuaternionType::FromRotationMatrix( rInitialLocalMatrix );

      mBox.UpperPoint = rUpperPoint;
      mBox.LowerPoint = rLowerPoint;

      mBox.Center = 0.5 * ( rUpperPoint + rLowerPoint );
      mBox.Radius = 0.5 * norm_2(rUpperPoint-rLowerPoint); 

      mBox.Velocity        = rVelocity;
      mBox.AngularVelocity = rAngularVelocity;
      
      //set to local frame
      this->MapToLocalFrame(mBox.InitialLocalQuaternion,mBox);

      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, mBox.Velocity);
      BeamMathUtilsType::MapToCurrentLocalFrame(mBox.InitialLocalQuaternion, mBox.AngularVelocity);

      mBox.SetInitialValues();

      mBox.Print();

      mRigidBodyCenterSupplied = false;

      KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    unsigned int BoxNoseSearch(const PointType& rPoint)
    {

      KRATOS_TRY

      double MinimumDistance       =std::numeric_limits<double>::max();
      double MinimumDistanceRadius =std::numeric_limits<double>::max();

      //std::cout<<" Minimum Distance "<<MinimumDistance<<std::endl;

      unsigned int NumberBoxNoses = mBoxNoses.size();

      unsigned int SelectedNose         = 0;
      unsigned int SelectedNoseDistance = 0;
      unsigned int SelectedNoseRadius   = 0;    

      std::vector<double> NoseDistanceVector (NumberBoxNoses);

      for(unsigned int i=0; i<NumberBoxNoses; i++)
	{

	  //based on distance
	  PointType Distance = (mBoxNoses[i].Center - rPoint);
	  Distance[2] = 0; //2D
	  double NoseDistance = norm_2( Distance );  
	  double NoseDistanceRadius = norm_2( Distance );  

	  if( NoseDistance > 0 ){

	    // //CORRECTION START: add a correction by velocity component for critical points
	    
	    // //based on center movement
	    // PointType VelocityProjection = this->mMovement.Velocity;
	    // double NormVelocity = norm_2(VelocityProjection);
	    // if( NormVelocity > 0 )
	    //   VelocityProjection /= NormVelocity;

	    // Distance += (NoseDistance * 0.01) * VelocityProjection;
	    // NoseDistance = norm_2( Distance );

	    // //CORRECTION END
	    
	    //CORRECTION START: add a correction by radius component for critical points
	    Distance = ( 1.0 - (mBoxNoses[i].Radius)/NoseDistance ) * Distance;
	    NoseDistanceRadius = norm_2( Distance );
	    
	    //CORRECTION END

	  }

	  NoseDistanceVector[i] = NoseDistance;


	  if( NoseDistance < MinimumDistance ){
	    MinimumDistance = NoseDistance;
	    SelectedNoseDistance    = i;
	    //std::cout<<" SelectedNoseDistance: "<<SelectedNoseDistance<<" MinimumDistance "<<MinimumDistance<<std::endl;
	  }

	  if( NoseDistanceRadius < MinimumDistanceRadius ){
	    MinimumDistanceRadius = NoseDistanceRadius;
	    SelectedNoseRadius    = i;
	    //std::cout<<" SelectedNoseDistanceRadiius: "<<SelectedNoseRadius<<" MinimumDistanceRadius "<<MinimumDistanceRadius<<std::endl;
	  }

	  
	}


      if( SelectedNoseDistance ==  SelectedNoseRadius ){

	SelectedNose = SelectedNoseRadius;

      }
      else{
	
	if( mBoxNoses[SelectedNoseDistance].Convexity > mBoxNoses[SelectedNoseRadius].Convexity ){

	  
	  if( NoseDistanceVector[SelectedNoseDistance] > mBoxNoses[SelectedNoseDistance].Radius ){


	    PointType TipPoint =  GetTipPoint( mBoxNoses[SelectedNoseRadius] ); //slave point : convexity -

	    double sign = GetOrientation( rPoint, mBoxNoses[SelectedNoseDistance].Center, mBoxNoses[SelectedNoseRadius].Center, TipPoint );
   
	    if( sign > 0 )
	      SelectedNose = SelectedNoseRadius;
	    else
	      SelectedNose = SelectedNoseDistance;

	    //std::cout<<" SELECTED NOSE A "<<SelectedNose<<" sign "<<sign<<std::endl;

 	    // if ( NoseDistanceVector[SelectedNoseRadius] > mBoxNoses[SelectedNoseRadius].Radius + (mBoxNoses[SelectedNoseDistance].Radius * 0.02) ){
	    //   SelectedNose = SelectedNoseDistance;
	    // }
	    // else{
	    //   SelectedNose = SelectedNoseRadius;
	    // }

	  }
	  else{

	    SelectedNose = SelectedNoseDistance;
	  }


	}
	else if( mBoxNoses[SelectedNoseDistance].Convexity < mBoxNoses[SelectedNoseRadius].Convexity ){


	  if( mBoxNoses[SelectedNoseDistance].Radius < mBoxNoses[SelectedNoseRadius].Radius ){

  
	    PointType TipPoint =  GetTipPoint( mBoxNoses[SelectedNoseRadius] ); //slave point : convexity +

	    double sign = GetOrientation( rPoint, mBoxNoses[SelectedNoseDistance].Center, mBoxNoses[SelectedNoseRadius].Center, TipPoint );
   
	    if( sign > 0 )
	      SelectedNose = SelectedNoseRadius;
	    else
	      SelectedNose = SelectedNoseDistance;

	    //std::cout<<" SELECTED NOSE B "<<SelectedNose<<" sign "<<sign<<std::endl;

	    // if ( NoseDistanceVector[SelectedNoseDistance] < mBoxNoses[SelectedNoseDistance].Radius + (mBoxNoses[SelectedNoseDistance].Radius * 0.01) ){
	    //   SelectedNose = SelectedNoseDistance;
	    // }
	    // else{
	    //   SelectedNose = SelectedNoseRadius;
	    // }


	  }
	  else{

	    SelectedNose = SelectedNoseDistance;

	  }


	}
	else{
	  
	  SelectedNose = SelectedNoseRadius;
	}

      }
	             
      return SelectedNose;

      KRATOS_CATCH("")
      
    }


    //************************************************************************************
    //************************************************************************************

    double GetOrientation(const PointType& rPoint, const PointType& rMasterCenter, const PointType& rSlaveCenter, const PointType& rSlaveTipPoint )
    {
      
      PointType DistanceToPoint  = (rPoint - rMasterCenter);
      DistanceToPoint[2] = 0; //2D
      PointType DistanceToCenter = (rSlaveCenter - rMasterCenter);
      DistanceToCenter[2] = 0; //2D
      PointType DistanceToTip    = (rSlaveTipPoint - rMasterCenter);
      DistanceToTip[2] = 0; //2D

      PointType ReferenceOrientation = MathUtils<double>::CrossProduct( DistanceToTip, DistanceToCenter );
      PointType Orientation = MathUtils<double>::CrossProduct( DistanceToPoint, DistanceToCenter );

      double sign = (Orientation[2] * ReferenceOrientation[2]);

      return sign;
    }


    //************************************************************************************
    //************************************************************************************


    PointType GetTipPoint(const BoxNoseVariables& rWallNose)
    {

      double pi = 3.141592654;

      //-----------
	    
      PointType RakePoint(3);
	    
      RakePoint[0] = rWallNose.Center[0] - rWallNose.Radius * sin(rWallNose.RakeAngle);
      RakePoint[1] = rWallNose.Center[1] + rWallNose.Radius * cos(rWallNose.RakeAngle);
      RakePoint[2] = 0;
      
      PointType ClearancePoint(3);
	    
      ClearancePoint[0] = rWallNose.Center[0] - rWallNose.Radius * sin(rWallNose.ClearanceAngle);
      ClearancePoint[1] = rWallNose.Center[1] + rWallNose.Radius * cos(rWallNose.ClearanceAngle);
      ClearancePoint[2] = 0;
      
      PointType TipPoint(3);
      
      TipPoint    = ( RakePoint - rWallNose.Center ) + ( ClearancePoint - rWallNose.Center );
      TipPoint[2] = 0; //2D
      TipPoint   *= ( rWallNose.Radius/norm_2(TipPoint) );
      
      //open angle to get the correct tip direction
      double OpenAngle = ( rWallNose.ClearanceAngle - rWallNose.RakeAngle ) * ( 180.0 / pi ) - 90;
      
      if( OpenAngle < 90 )
	TipPoint  = rWallNose.Center + TipPoint;
      else
	TipPoint  = rWallNose.Center - TipPoint;
      
      TipPoint[2] = 0; //2D

      return TipPoint;

    }

    //************************************************************************************
    //************************************************************************************


    ContactFace ContactSearch(const PointType& rPoint, const double& rRadius, const BoxNoseVariables& rWallNose)
    {

      KRATOS_TRY

      ContactFace Face = FreeSurface;
           
      double FaceR = CalculateRakeFace( FaceR, rPoint, rRadius, rWallNose );
      double FaceT = CalculateTipFace( FaceT, rPoint, rRadius, rWallNose );
      double FaceC = CalculateClearanceFace( FaceC, rPoint, rRadius, rWallNose ); 
	
      double Face1=0, Face2=0, Face3=0;
      CalculateAuxiliarFaces( Face1, Face2, Face3, rPoint, rRadius, rWallNose );

      
      // bool node_in =false;
      // if( rWallNose.Convexity == 1 && (rPoint[0]>=95 && rPoint[1]>=6.9) )
      //  	node_in = true;

      // if( node_in ){
      // 	std::cout<<" Point : "<<rPoint<<std::endl;
      // 	std::cout<<" [ FaceR: "<<FaceR<<"; FaceT: "<<FaceT<<"; FaceC: "<<FaceC<<" ] "<<std::endl; 
      // 	std::cout<<" [ Face1: "<<Face1<<"; Face2: "<<Face2<<"; Face3: "<<Face3<<" ] "<<std::endl;
      // }

      if(rWallNose.Convexity == 1){

	if(FaceR>0 && Face3<0 && Face1<0){
	  Face = RakeSurface;
	}
	else if(FaceT>0 && Face1>0 && Face2>0){
	  Face = TipSurface;
	}
	else if(FaceC>0 && Face3>0 && Face2<0){
	  Face = ClearanceSurface;
	}
	else{
	  Face = FreeSurface;
	}
	
      }
      else{
	
	if(FaceR<0 && Face3<0 && Face1<0){
	  Face = RakeSurface;
	}
	else if(FaceT<0 && Face1>0 && Face2>0){
	  Face = TipSurface;
	}
	else if(FaceC<0 && Face3>0 && Face2<0){
	  Face = ClearanceSurface;
	}
	else{
	  Face = FreeSurface;
	}


      }

      return Face;

      KRATOS_CATCH( "" )

   }


   //************************************************************************************
   //************************************************************************************


   void PointFaceEvaluation(const PointType& rPoint, const PointType& rLinePoint, const  double& rTangentAngle, PointType& rPointFace)
    {
      KRATOS_TRY

	// if( rTangentAngle != 0 ){
	  
	//   rPointFace[0] = rPoint[0] - ( 1.0/rTangentAngle ) * ( rPoint[1] - rLinePoint[1] ) - rLinePoint[0];
	//   rPointFace[1] = rPoint[1] - (   rTangentAngle   ) * ( rPoint[0] - rLinePoint[0] ) - rLinePoint[1];

	// }
	// else{
	  
      rPointFace[0] = rPoint[0] - rLinePoint[0];
      rPointFace[1] = rPoint[1] - rLinePoint[1];  
	// }	

      rPointFace [2] = 0;
      
      KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************


    double& CalculateRakeFace(double& Face, const PointType& rPoint, const double& rRadius, const BoxNoseVariables& rWallNose)
    {
      KRATOS_TRY
	    

      PointType PointFace(3);

      PointType CenterFace(3);

      PointType RakePoint(3);

      RakePoint[0] = rWallNose.Center[0] - rRadius * sin(rWallNose.RakeAngle);
      RakePoint[1] = rWallNose.Center[1] + rRadius * cos(rWallNose.RakeAngle);
      RakePoint[2] = 0;
 

      // bool node_in =false;
      // if( rWallNose.Convexity == 1 && (rPoint[0]>=95 && rPoint[1]>=6.9) )
      //  	node_in = true;

      // if( node_in ) 
      // 	std::cout<<" RakePoint "<<RakePoint<<" Radius "<<rWallNose.Radius<<" Angle "<<rWallNose.RakeAngle<<" Center "<<rWallNose.Center<<std::endl;


      PointFaceEvaluation(rPoint, RakePoint, rWallNose.TangentRakeAngle, PointFace);
      PointFaceEvaluation(rWallNose.Center, RakePoint, rWallNose.TangentRakeAngle, CenterFace);

      Face = inner_prod(PointFace,CenterFace);

      // PointFace[0] *= CenterFace[0];
      // PointFace[1] *= CenterFace[1];

      // double MaximumDistance=-1;

      // for(unsigned int i=0; i<2; i++)
      // 	{
      // 	  double Distance = fabs( PointFace[i] );
      // 	  if( Distance > MaximumDistance ){
      // 	    MaximumDistance = Distance;
      // 	  }
      // 	}

      // double zero_tol = MaximumDistance * 1e-3;
      

      // //the sign is evaluated (+) accepted and  (-) rejected
      // if( fabs(PointFace[0]) > zero_tol && fabs(PointFace[1]) > zero_tol ){
      // 	Face = PointFace[0] * PointFace[1];
      // }
      // else if( fabs(PointFace[0]) > zero_tol && fabs(PointFace[1]) < zero_tol ){
      // 	Face = PointFace[0];
      // }
      // else if( fabs(PointFace[0]) < zero_tol && fabs(PointFace[1]) > zero_tol ){
      // 	Face = PointFace[1];
      // }
      // else{
      // 	std::cout<<" Critical point reached and accepted on Rake Face: "<<PointFace<<std::endl;
      // 	Face = +1;
      // }

      //Face (+) rPoint is "in" :: (-) rPoint is "out" the rake part of the nose
      return Face;
    
      KRATOS_CATCH( "" )

    }


    //************************************************************************************
    //************************************************************************************

    double& CalculateClearanceFace(double& Face, const PointType& rPoint, const double& rRadius, const BoxNoseVariables& rWallNose)
    {
      KRATOS_TRY

      PointType PointFace(3);

      PointType CenterFace(3);

      PointType ClearancePoint(3);

      ClearancePoint[0] = rWallNose.Center[0] - rRadius * sin(rWallNose.ClearanceAngle);
      ClearancePoint[1] = rWallNose.Center[1] + rRadius * cos(rWallNose.ClearanceAngle);
      ClearancePoint[2] = 0; 

      // bool node_in =false;
      // if( rWallNose.Convexity == 1 && (rPoint[0]>=95 && rPoint[1]>=6.9) )
      //  	node_in = true;

      // if( node_in ) 
      // 	std::cout<<" ClearancePoint "<<ClearancePoint<<" Radius "<<rWallNose.Radius<<" Angle "<<rWallNose.ClearanceAngle<<" Center "<<rWallNose.Center<<std::endl;

      PointFaceEvaluation(rPoint, ClearancePoint, rWallNose.TangentClearanceAngle, PointFace);
      PointFaceEvaluation(rWallNose.Center, ClearancePoint, rWallNose.TangentClearanceAngle, CenterFace);

      Face = inner_prod(PointFace,CenterFace);

      // PointFace[0] *= CenterFace[0];
      // PointFace[1] *= CenterFace[1];

      // double MaximumDistance=-1;

      // for(unsigned int i=0; i<2; i++)
      // 	{
      // 	  double Distance = fabs( PointFace[i] );
      // 	  if( Distance > MaximumDistance ){
      // 	    MaximumDistance = Distance;
      // 	  }
      // 	}

      // double zero_tol = MaximumDistance * 1e-3;
      
      // //the sign is evaluated (+) accepted and  (-) rejected
      // if( fabs(PointFace[0]) > zero_tol && fabs(PointFace[1]) > zero_tol ){
      // 	Face = PointFace[0] * PointFace[1];
      // }
      // else if( fabs(PointFace[0]) > zero_tol && fabs(PointFace[1]) < zero_tol ){
      // 	Face = PointFace[0];
      // }
      // else if( fabs(PointFace[0]) < zero_tol && fabs(PointFace[1]) > zero_tol ){
      // 	Face = PointFace[1];
      // }
      // else{
      // 	std::cout<<" Critical point reached and accepted on Clearance Face :"<<PointFace<<std::endl;
      // 	Face = +1;
      // }
	  

      //Face (+) rPoint is "in" :: (-) rPoint is "out" the clearance part of the nose
      return Face;

    
      KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************


    double& CalculateTipFace(double& Face, const PointType& rPoint, const double& rRadius, const BoxNoseVariables& rWallNose)
    {
      KRATOS_TRY
      
      Face = ( rRadius * rRadius ) - pow( (rPoint[0] - rWallNose.Center[0]), 2 ) - pow( (rPoint[1] - rWallNose.Center[1]), 2 ); 
      
      //Face (+) rPoint is "in" :: (-) rPoint is "out" the nose tip
      return Face;	
    
      KRATOS_CATCH( "" )
    }


    //************************************************************************************
    //************************************************************************************

    void CalculateLineProjection(double& rFace, const PointType& rPoint, const PointType& rReferencePoint, const PointType& rCenterPoint, const PointType& rLinePoint )
    {

      PointType Line(3);

      PointType NormalLine(3);

      //rCenterPoint and rLinePoint define the line:
      Line  = (rLinePoint - rCenterPoint);
      Line *= (1.0/norm_2(Line));
      
      //normal to the line
      NormalLine[0] = -Line[1];
      NormalLine[1] =  Line[0];
      NormalLine[2] =  0;

      //compare the sense of the direction of Line1 and Line2
      PointType Line1 = (rReferencePoint - rCenterPoint);
      PointType Line2 = (rPoint - rCenterPoint);

      //check the Projection with the Line
      double ReferenceDirection = inner_prod(Line1,NormalLine);
      double PointDirection     = inner_prod(Line2,NormalLine);

      //the sign is evaluated (+) accepted and  (-) rejected
      if( PointDirection != 0 ){
	rFace = PointDirection * ReferenceDirection;
      }
      else{
	std::cout<<" Critical point reached and accepted on a Face "<<std::endl;
	std::cout<<" LINE "<<rLinePoint<<" CENTER "<<rCenterPoint<<" REFERENCE "<<rReferencePoint<<" POINT "<<rPoint<<std::endl;
	rFace = +1;
      }

    }


    //************************************************************************************
    //************************************************************************************


    void CalculateAuxiliarFaces(double& rFace1, double& rFace2,  double& rFace3, const PointType& rPoint, const double& rRadius, const BoxNoseVariables& rWallNose)
    {
      KRATOS_TRY

      double pi = 3.141592654;

      //-----------

      PointType RakePoint(3);

      RakePoint[0] = rWallNose.Center[0] - rRadius * sin(rWallNose.RakeAngle);
      RakePoint[1] = rWallNose.Center[1] + rRadius * cos(rWallNose.RakeAngle);
      RakePoint[2] = 0;

      PointType ClearancePoint(3);

      ClearancePoint[0] = rWallNose.Center[0] - rRadius * sin(rWallNose.ClearanceAngle);
      ClearancePoint[1] = rWallNose.Center[1] + rRadius * cos(rWallNose.ClearanceAngle);
      ClearancePoint[2] = 0;

      PointType TipPoint(3);

      TipPoint  = ( RakePoint - rWallNose.Center ) + ( ClearancePoint - rWallNose.Center );
      TipPoint *= ( rRadius/norm_2(TipPoint) );
      
      //open angle to get the correct tip direction
      double OpenAngle = ( rWallNose.ClearanceAngle - rWallNose.RakeAngle ) * ( 180.0 / pi ) - 90;

      if( OpenAngle < 90 )
	TipPoint  = rWallNose.Center + TipPoint;
      else
	TipPoint  = rWallNose.Center - TipPoint;

      // std::cout<<" Center "<<rWallNose.Center<<" Radius "<<rRadius<<std::endl;
      // std::cout<<" RakePoint "<<RakePoint<<" x "<<-sin(rWallNose.RakeAngle)<<" y "<<cos(rWallNose.RakeAngle)<<std::endl;
      // std::cout<<" ClearancePoint "<<ClearancePoint<<" x "<<-sin(rWallNose.ClearanceAngle)<<" y "<<cos(rWallNose.ClearanceAngle)<<std::endl;
      // std::cout<<" TipPoint "<<TipPoint<<std::endl;
      

      //-----------
      
 
      //Center to RakePoint line  (TipPoint:= NoseCenter)  (rFace1)
      
      CalculateLineProjection(rFace1, rPoint, TipPoint, rWallNose.Center, RakePoint);

      //Face1 (+) rPoint is "in" the same part as the TipPoint :: (-) rPoint is "out" of the nose tip
      

      //-----------


      //Center to ClearancePoint line (TipPoint:= NoseCenter)  (rFace2)

      CalculateLineProjection(rFace2, rPoint, TipPoint, rWallNose.Center, ClearancePoint);

      //Face2 (+) rPoint is "in" the same part as the TipPoint :: (-) rPoint is "out" of the same part of the ClearancePoint



      //-----------

      //Center to TipPoint line (ClearancePoint:= NoseCenter) (rFace3)

      CalculateLineProjection(rFace3, rPoint, ClearancePoint, rWallNose.Center, TipPoint);

      //Face3 (+) rPoint is "in" the same part as the ClearancePoint :: (-) rPoint is "in" the same part of the RakePoint


      KRATOS_CATCH( "" )
    }


    //************************************************************************************
    //************************************************************************************


    bool CalculateRakeSurface(const PointType& rPoint, const double & rRadius, double& rGapNormal, PointType& rNormal, const BoxNoseVariables& rWallNose)
    {
      KRATOS_TRY
     
      rNormal  = ZeroVector(3);
 
      //1.-compute contact normal
      rNormal[0] = -sin(rWallNose.RakeAngle);
      rNormal[1] =  cos(rWallNose.RakeAngle);
      rNormal[2] = 0;

      rNormal   *= rWallNose.Convexity; 

      //2.-compute point projection
      PointType RakePoint(3);

      RakePoint[0] = rWallNose.Center[0] - rRadius * sin(rWallNose.RakeAngle);
      RakePoint[1] = rWallNose.Center[1] + rRadius * cos(rWallNose.RakeAngle);
      RakePoint[2] = 0;

      //3.-compute gap
      rGapNormal = inner_prod((rPoint - RakePoint), rNormal);

       
      if(rGapNormal<0)
	return true;
      else
	return false;

      KRATOS_CATCH( "" )
    }
  
    //************************************************************************************
    //************************************************************************************

    bool CalculateTipSurface(const PointType& rPoint, const double& rRadius, double& rGapNormal, PointType& rNormal, const BoxNoseVariables& rWallNose)
    {
      KRATOS_TRY

      rNormal  = ZeroVector(3);

      //1.-compute point projection
      PointType Projection(3);
      Projection = rRadius * ( (rPoint-rWallNose.Center)/ norm_2(rPoint-rWallNose.Center) ) + rWallNose.Center;
      
      //2.-compute contact normal
      rNormal = (Projection-rWallNose.Center)/rRadius;

      rNormal   *= rWallNose.Convexity; 

      //3.-compute gap
      if( norm_2(rWallNose.Center-rPoint) <= rRadius ){
	rGapNormal = (-1) * norm_2(rPoint - Projection);
      }
      else{
	rGapNormal = norm_2(Projection - rPoint);
      }
           
      rGapNormal *= rWallNose.Convexity;

      if(rGapNormal<0)
	return true;
      else
	return false;

    
      KRATOS_CATCH( "" )
	}


    //************************************************************************************
    //************************************************************************************


    bool CalculateClearanceSurface(const PointType& rPoint, const double& rRadius, double& rGapNormal, PointType& rNormal, const BoxNoseVariables& rWallNose)
    {
      KRATOS_TRY

      rNormal  = ZeroVector(3);
      
      //1.-compute contact normal
      rNormal[0] = -sin(rWallNose.ClearanceAngle);
      rNormal[1] =  cos(rWallNose.ClearanceAngle);
      rNormal[2] = 0;

      rNormal   *= rWallNose.Convexity; 

      //2.-compute point projection
      PointType ClearancePoint(3);

      ClearancePoint[0] = rWallNose.Center[0] - rRadius * sin(rWallNose.ClearanceAngle);
      ClearancePoint[1] = rWallNose.Center[1] + rRadius * cos(rWallNose.ClearanceAngle);
      ClearancePoint[2] = 0;
 
      //3.-compute gap
      rGapNormal = inner_prod((rPoint - ClearancePoint), rNormal);

      if(rGapNormal<0)
	return true;
      else
	return false;
     
      KRATOS_CATCH( "" )
     }


    //************************************************************************************
    //************************************************************************************

    static inline double inner_prod(const PointType& a, const PointType& b)
    {
        double temp =a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
        return temp;
    }

    //************************************************************************************
    //************************************************************************************

    static inline double norm_2(const PointType& a)
    {
        double temp = pow(a[0],2) + pow(a[1],2) + pow(a[2],2);
        temp = sqrt(temp);
        return temp;
    }

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


}; // Class CompoundNosesBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_COMPOUND_NOSES_BOUNDING_BOX_H_INCLUDED  defined 


