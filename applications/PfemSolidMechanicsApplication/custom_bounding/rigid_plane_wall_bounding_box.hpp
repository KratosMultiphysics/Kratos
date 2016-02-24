//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_RIGID_PLANE_WALL_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_RIGID_PLANE_WALL_BOUNDING_BOX_H_INCLUDED

// External includes

// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/model_part.h"

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

class RigidPlaneWallBoundingBox
  : public SpatialBoundingBox
{
protected:

   typedef struct
   {

     TPointType  Point;      // plane point
     TPointType  Normal;     // plane normal

     
   public:
     
     void clear()
     {

       Point.clear();
       Normal.clear();

      }


   } PlaneVariables;


public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RigidPlaneWallBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( RigidPlaneWallBoundingBox );

    typedef Vector TPointType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RigidPlaneWallBoundingBox() : SpatialBoundingBox()
    {
        std::cout<< "Calling Rigid Plane Wall BBX empty constructor" <<std::endl;
    }

    // General Wall constructor
    RigidPlaneWallBoundingBox( int Label,
			       int Convexity,
			       TPointType Point,
			       TPointType Normal,
			       TPointType Velocity,
			       TPointType AngularVelocity,
			       TPointType RotationCenter)
    {
                  
      std::cout<<" [--PLANE WALL--]["<<Label<<"]"<<std::endl;
      
      mBox.OriginalCenter = Point;
      mBox.Center    = Point;
      mBox.Convexity = Convexity;
      
      mPlane.Point  = Point;

      if( norm_2(Normal) )
	mPlane.Normal = Normal/norm_2(Normal);
      else
	std::cout<<" [ERROR: Normal is Zero]"<<std::endl;


      std::cout<<"  [Convexity:"<<mBox.Convexity<<std::endl;
      std::cout<<"  [Point:"<<mPlane.Point<<std::endl;
      std::cout<<"  [Normal:"<<mPlane.Normal<<std::endl;
      std::cout<<" [--------] "<<std::endl;
      
      this->mMovement.Label                   = Label;
      this->mMovement.Velocity                = Velocity;
      this->mMovement.AngularVelocity         = AngularVelocity;
      this->mMovement.OriginalRotationCenter  = RotationCenter;
      this->mMovement.RotationCenter          = RotationCenter;

    }


    /// Assignment operator.
    RigidPlaneWallBoundingBox& operator=(RigidPlaneWallBoundingBox const& rOther)
    {
      SpatialBoundingBox::operator=(rOther);
      return *this;
    }


    /// Copy constructor.
    RigidPlaneWallBoundingBox(RigidPlaneWallBoundingBox const& rOther) 
    :SpatialBoundingBox(rOther)
    {
    }

    /// Destructor.
    virtual ~RigidPlaneWallBoundingBox() {};


    ///@}
    ///@name Operators
    ///@{



    ///@}
    ///@name Operations
    ///@{

    double GetRadius()
    {
        return mBox.Radius;
    }

    void SetRadius(double& rRadius)
    {
      //used to set a comparisson radius
    }  

    virtual void UpdatePosition(double & rTime)
    {
      
      SpatialBoundingBox::UpdatePosition(rTime);
      
      mPlane.Point =  mBox.OriginalCenter + this->mMovement.Velocity * rTime;

    }


    //************************************************************************************
    //************************************************************************************
   

    bool IsInside (const TPointType& rPoint, double& rCurrentTime, int & ContactFace, double Radius = 0)
    {
      
      ContactFace = 0; //FreeSurface
      
      return IsInside(rPoint,rCurrentTime,Radius);
      
    } 

    //************************************************************************************
    //************************************************************************************
   

    bool IsInside (const TPointType& rPoint, double& rCurrentTime, double Radius = 0)
    {
      
      bool is_inside = false;
      
      TPointType  PlanePoint = mPlane.Point;
      
      if( mBox.Convexity == 1 )
	PlanePoint += mPlane.Normal * 0.1; //increase the bounding box 

      if( mBox.Convexity == -1 )
       	PlanePoint -= mPlane.Normal * 0.1; //decrease the bounding box 

      is_inside = ContactSearch(rPoint, PlanePoint);


      return is_inside;
      
    } 


   //************************************************************************************
    //************************************************************************************
   
    bool IsInside(const TPointType& rPoint, double& rGapNormal, double& rGapTangent, TPointType& rNormal, TPointType& rTangent, int& ContactFace, double Radius = 0)
    {

      ContactFace = 0; //FreeSurface
      
      return IsInside(rPoint,rGapNormal,rGapTangent,rNormal,rTangent,Radius);
      
    }


    //************************************************************************************
    //************************************************************************************
   
    bool IsInside(const TPointType& rPoint, double& rGapNormal, double& rGapTangent, TPointType& rNormal, TPointType& rTangent, double Radius = 0)
    {
      bool is_inside = false;

      rGapNormal  = 0;
      rGapTangent = 0;
      rNormal.clear();
      rTangent.clear();

      is_inside = ContactSearch(rPoint,rGapNormal,rGapTangent,rNormal,rTangent);

      return is_inside;
      
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
        return "RigidPlaneWallBoundingBox";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << this->mBox.HighPoint << " , " << this->mBox.LowPoint;
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


    bool ContactSearch(const TPointType& rPoint, const TPointType& rPlanePoint)
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

    //Plane
    bool ContactSearch(const TPointType& rPoint, double& rGapNormal, double& rGapTangent, TPointType& rNormal, TPointType& rTangent)
    {
      KRATOS_TRY
     
      rNormal  = ZeroVector(3);
      rTangent = ZeroVector(3);
 
      //1.-compute contact normal
      rNormal = mPlane.Normal;

      rNormal *= mBox.Convexity; 

      //2.-compute gap
      rGapNormal = inner_prod((rPoint - mPlane.Point), rNormal);


      //3.-compute contact tangent
      rTangent = (rPoint - mPlane.Point) - rGapNormal * rNormal;
      
      if(norm_2(rTangent))
	rTangent/= norm_2(rTangent);

       
      if(rGapNormal<0)
	return true;
      else
	return false;

      KRATOS_CATCH( "" )
    }



    //************************************************************************************
    //************************************************************************************

    static inline double inner_prod(const TPointType& a, const TPointType& b)
    {
        double temp =a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
        return temp;
    }

    //************************************************************************************
    //************************************************************************************

    static inline double norm_2(const TPointType& a)
    {
        double temp = pow(a[0],2) + pow(a[1],2) + pow(a[2],2);
        temp = sqrt(temp);
        return temp;
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


}; // Class RigidPlaneWallBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TPointType, class TPointerType>
inline std::istream& operator >> (std::istream& rIStream,
                                  RigidPlaneWallBoundingBox& rThis);

/// output stream function
template<class TPointType, class TPointerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RigidPlaneWallBoundingBox& rThis)
{
    // rThis.PrintInfo(rOStream);
    // rOStream << std::endl;
    // rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_RIGID_PLANE_WALL_BOUNDING_BOX_H_INCLUDED  defined 


