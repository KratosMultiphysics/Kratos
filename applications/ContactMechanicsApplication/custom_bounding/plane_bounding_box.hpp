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

class KRATOS_API(CONTACT_MECHANICS_APPLICATION) PlaneBoundingBox
  : public SpatialBoundingBox
{
public:

    //typedef bounded_vector<double, 3>                     PointType;
    typedef array_1d<double, 3>                             PointType;
    typedef ModelPart::NodeType                              NodeType;
    typedef ModelPart::NodesContainerType          NodesContainerType;
    typedef NodesContainerType::Pointer     NodesContainerTypePointer;


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
                 "velocity": [0.0, 0.0, 0.0]
                
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

    void UpdateBoxPosition(const double & rCurrentTime)
    {

      KRATOS_TRY

      PointType Displacement  =  this->GetBoxDisplacement(rCurrentTime);

      mBox.UpdatePosition(Displacement);

      mPlane.Point = mBox.Center;

      KRATOS_CATCH("")
      
    }


    //************************************************************************************
    //************************************************************************************
   

    bool IsInside (const PointType& rPoint, double& rCurrentTime, double Radius = 0)
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
    
    bool IsInside(BoundingBoxParameters& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      bool is_inside = false;

      rValues.SetContactFace(0);
      rValues.SetRadius(0.0);
      
      is_inside = ContactSearch(rValues, rCurrentProcessInfo);

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
        return "PlaneBoundingBox";
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

    PlaneVariables mPlane;

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


