//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_RIGID_CIRCLE_WALL_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_RIGID_CIRCLE_WALL_BOUNDING_BOX_H_INCLUDED

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

    This Box represents a 2D wall composed by a circle
    
    A convexity parameter is given to determine which side of each nose is considered 
    the internal or external boundary

    This bounding box is essentially used for rigid wall contact purposes
*/

class RigidCircleWallBoundingBox
  : public SpatialBoundingBox
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RigidCircleWallBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( RigidCircleWallBoundingBox );

    typedef Vector TPointType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RigidCircleWallBoundingBox() : SpatialBoundingBox()
    {
        std::cout<< "Calling Rigid Circle Wall BBX empty constructor" <<std::endl;
    }

    // General Wall constructor
    RigidCircleWallBoundingBox( int Label,
				int Convexity,
				double Radius,
				TPointType Center,
				TPointType Velocity,
				TPointType AngularVelocity,
				TPointType RotationCenter)
    {
            
      
      std::cout<<" [--CIRCLE WALL--]["<<Label<<"]"<<std::endl;
      
      mBox.OriginalCenter = Center;
      mBox.Center = Center;
      mBox.Radius = Radius;
      mBox.Convexity = Convexity;

      std::cout<<"  [Convexity:"<<mBox.Convexity<<std::endl;
      std::cout<<"  [Radius:"<<mBox.Radius<<std::endl;
      std::cout<<"  [Center:"<<mBox.Center<<std::endl;
      std::cout<<" [--------] "<<std::endl;
      
      this->mMovement.Label                   = Label;
      this->mMovement.Velocity                = Velocity;
      this->mMovement.AngularVelocity         = AngularVelocity;
      this->mMovement.OriginalRotationCenter  = RotationCenter;
      this->mMovement.RotationCenter          = RotationCenter;

    }


    /// Assignment operator.
    RigidCircleWallBoundingBox& operator=(RigidCircleWallBoundingBox const& rOther)
    {
      SpatialBoundingBox::operator=(rOther);
      return *this;
    }


    /// Copy constructor.
    RigidCircleWallBoundingBox(RigidCircleWallBoundingBox const& rOther) 
    :SpatialBoundingBox(rOther)
    {
    }

    /// Destructor.
    virtual ~RigidCircleWallBoundingBox() {};


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
      
    }


    //************************************************************************************
    //************************************************************************************
   

    bool IsInside (const TPointType& rPoint, double& rCurrentTime, int & ContactFace, double Radius = 0)
    {
      
      ContactFace = 2; //TipSurface
      
      return IsInside(rPoint,rCurrentTime,Radius);
      
    } 

    //************************************************************************************
    //************************************************************************************
   

    bool IsInside (const TPointType& rPoint, double& rCurrentTime, double Radius = 0)
    {
      
      bool is_inside = false;
      
      double CircleRadius = mBox.Radius;

      if( mBox.Convexity == 1)
	CircleRadius *= 1.25; //increase the bounding box 

      if( mBox.Convexity == -1)
       	CircleRadius *= 0.75; //decrease the bounding box 

      is_inside = ContactSearch(rPoint, CircleRadius);


      return is_inside;
      
    } 


    //************************************************************************************
    //************************************************************************************
   
    bool IsInside(const TPointType& rPoint, double& rGapNormal, double& rGapTangent, TPointType& rNormal, TPointType& rTangent, int& ContactFace, double Radius = 0)
    {

      ContactFace = 2; //TipSurface
      
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

      double CircleRadius = mBox.Radius;

      is_inside = ContactSearch(rPoint,CircleRadius,rGapNormal,rGapTangent,rNormal,rTangent);

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
        return "RigidCircleWallBoundingBox";
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


    bool ContactSearch(const TPointType& rPoint, const double& rRadius)
    {

      KRATOS_TRY

      //1.-compute point projection
      TPointType Projection(3);
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

    //Circle
    bool ContactSearch(const TPointType& rPoint, const double& rRadius, double& rGapNormal, double& rGapTangent, TPointType& rNormal, TPointType& rTangent)
    {
      KRATOS_TRY

      rNormal  = ZeroVector(3);
      rTangent = ZeroVector(3);

      //1.-compute point projection
      TPointType Projection(3);
      Projection = rRadius * ( (rPoint-mBox.Center)/ norm_2(rPoint-mBox.Center) ) + mBox.Center;
      
      //2.-compute contact normal
      rNormal = (Projection-mBox.Center)/rRadius;

      rNormal   *= mBox.Convexity; 

      rTangent[0] =  rNormal[1];
      rTangent[1] = -rNormal[0];
      rTangent[2] =  0;

      //3.-compute gap
      if( norm_2(mBox.Center-rPoint) <= rRadius ){
	rGapNormal = (-1) * norm_2(rPoint - Projection);
      }
      else{
	rGapNormal = norm_2(Projection - rPoint);
      }
           
      rGapNormal *= mBox.Convexity;

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


}; // Class RigidCircleWallBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TPointType, class TPointerType>
inline std::istream& operator >> (std::istream& rIStream,
                                  RigidCircleWallBoundingBox& rThis);

/// output stream function
template<class TPointType, class TPointerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RigidCircleWallBoundingBox& rThis)
{
    // rThis.PrintInfo(rOStream);
    // rOStream << std::endl;
    // rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_RIGID_CIRCLE_WALL_BOUNDING_BOX_H_INCLUDED  defined 


