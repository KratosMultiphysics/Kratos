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

class KRATOS_API(CONTACT_MECHANICS_APPLICATION) SphereBoundingBox
  : public SpatialBoundingBox
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SphereBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( SphereBoundingBox );

    //typedef bounded_vector<double, 3>                     PointType;
    typedef array_1d<double, 3>                             PointType;
    typedef ModelPart::NodeType                              NodeType;
    typedef ModelPart::NodesContainerType          NodesContainerType;
    typedef NodesContainerType::Pointer     NodesContainerTypePointer;

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

      std::cout<<" [--CIRCLE WALL--] "<<std::endl;
      
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
   
    bool IsInside (const PointType& rPoint, double& rCurrentTime, double Radius = 0)
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
    
    bool IsInside(BoundingBoxParameters& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      bool is_inside = false;

      rValues.SetContactFace(2);

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
        return "SphereBoundingBox";
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


