//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_RIGID_TOOL_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_RIGID_TOOL_BOUNDING_BOX_H_INCLUDED

// External includes

// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/model_part.h"

#include "custom_modelers/spatial_bounding_box.hpp"


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

class RigidToolBoundingBox
  : public SpatialBoundingBox
{
private:

   enum ContactFace{ FreeSurface=0, RakeSurface=1, TipSurface=2, ClearanceSurface=3 };

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RigidToolBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( RigidToolBoundingBox );

    typedef Vector TPointType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RigidToolBoundingBox() : SpatialBoundingBox()
    {
        std::cout<< "Calling empty constructor" <<std::endl;
    }

    RigidToolBoundingBox( double Radius,
			  double RakeAngle,
			  double ClearanceAngle,
			  TPointType  Center,
			  TPointType  Velocity) :  SpatialBoundingBox(Radius, RakeAngle, ClearanceAngle, Center, Velocity)
    {
      
      double pi = 3.141592654;

      this->mBox.Radius         = Radius;
      this->mBox.RakeAngle      = RakeAngle * pi / 180;
      this->mBox.ClearanceAngle = ClearanceAngle * pi / 180;

      this->mBox.m_factor       = tan(0.5*pi-this->mBox.RakeAngle);
      this->mBox.n_factor       = tan(this->mBox.ClearanceAngle);

      this->mBox.OriginalCenter = Center;
      this->mBox.Center         = Center;
      this->mBox.Velocity       = Velocity;

      std::cout<<" [TOOL:                        ] "<<std::endl;
      std::cout<<" [Radius:"<<this->mBox.Radius<<"            ] "<<std::endl;
      std::cout<<" [Center:"<<this->mBox.Center<<"  ] "<<std::endl;
      std::cout<<" [Velocity:"<<this->mBox.Velocity<<"      ] "<<std::endl;
      std::cout<<" [Rake:"<<this->mBox.RakeAngle<<"               ] "<<std::endl;
      std::cout<<" [Clearance:"<<this->mBox.ClearanceAngle<<"          ] "<<std::endl;
      std::cout<<" [m_factor:"<<this->mBox.m_factor<<"          ] "<<std::endl;
      std::cout<<" [n_factor:"<<this->mBox.n_factor<<"          ] "<<std::endl;

    }

    /// Assignment operator.
    RigidToolBoundingBox& operator=(RigidToolBoundingBox const& rOther)
    {
      SpatialBoundingBox::operator=(rOther);
      return *this;
    }


    /// Copy constructor.
    RigidToolBoundingBox(RigidToolBoundingBox const& rOther) 
    :SpatialBoundingBox(rOther)
    {
    }

    /// Destructor.
    virtual ~RigidToolBoundingBox() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    bool IsInside (const TPointType& rPoint, double& rCurrentTime)
    {
      bool is_inside = false;

      this->mBox.Radius *= 2; //increase the bounding box 

      switch( ContactSearch(rPoint) )
	{
	  
	case FreeSurface:      
	  is_inside = false;
	  break;
	case RakeSurface:      
	  is_inside = true;
	  break;
	case TipSurface:       
	  is_inside = true;
	  break;
	case ClearanceSurface: 
	  is_inside = true;
	  break;
	default:               
	  is_inside = false;
	  break;
	}

      this->mBox.Radius *= 0.5; //restore the bounding box

      return is_inside;
      
    } 

    //************************************************************************************
    //************************************************************************************
   
    bool IsInside(const TPointType& rPoint, double& rGapNormal, double& rGapTangent, TPointType& rNormal, TPointType& rTangent, int ContactFace = 0)
    {
      bool is_inside = false;

      rGapNormal  = 0;
      rGapTangent = 0;
      rNormal.clear();
      rTangent.clear();

      switch( ContactSearch(rPoint) )
	{	  
	case FreeSurface:      
	  is_inside = false;
	  break;
	case RakeSurface:      
	  is_inside = CalculateRakeSurface(rPoint, rGapNormal, rGapTangent, rNormal, rTangent);
	  break;
	case TipSurface:       
	  is_inside = CalculateTipSurface(rPoint, rGapNormal, rGapTangent, rNormal, rTangent);
	  break;
	case ClearanceSurface: 
	  is_inside = CalculateClearanceSurface(rPoint, rGapNormal, rGapTangent, rNormal, rTangent);
	  break;
	default:               
	  is_inside = false;
	  break;
	}

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
        return "RigidToolBoundingBox";
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

    ContactFace ContactSearch(const TPointType& rPoint)
    {

      KRATOS_TRY

      ContactFace Face = FreeSurface;
           
      double FaceR = CalculateRakeFace( FaceR, rPoint );
      double FaceT = CalculateTipFace( FaceT, rPoint );
      double FaceC = CalculateClearanceFace( FaceC, rPoint );      
	
      double Face1=0,Face2=0,Face3=0;
      CalculateAuxiliarFaces( Face1, Face2, Face3, rPoint );

      //The nodes in the wall tip, are marked as TO_SPLIT 
      //in order to be susceptible to refine
      //rPoint.Reset(TO_SPLIT);


      if(this->mBox.RakeAngle>0){
	if(FaceR<=0 && Face3>=0 && Face1>=0){
	  Face = RakeSurface;
	}
	else if(FaceT<=0 && Face3<0 && Face2<=0){
	  Face = TipSurface;
	  //It must be set to be able to refine boundaries later on REFINE
	  //rPoint.Set(TO_SPLIT);
	}
	else if(FaceC>=0 && Face2>=0 && Face1<0){
	  Face = ClearanceSurface;
	}
	else{
	  Face = FreeSurface;
	}
      }
      else if(this->mBox.RakeAngle==0){

	if(FaceR<=0 && Face3>=0 && Face1>=0){
	  Face = RakeSurface;
	}
	else if(FaceT<=0 && Face3<=0 && Face2<=0){
	  Face = TipSurface;
	  //It must be set to be able to refine boundaries later on REFINE
	  //rPoint.Set(TO_SPLIT);
	}
	else if(FaceC>=0 && Face2>=0 && Face1<=0){
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


    double& CalculateRakeFace(double& Face, const TPointType& rPoint)
    {
      KRATOS_TRY
	    
      Face = rPoint[1] - this->mBox.m_factor * ( rPoint[0] + this->mBox.Radius * cos(this->mBox.RakeAngle) - this->mBox.Center[0]) - this->mBox.Center[1] - this->mBox.Radius * sin(this->mBox.RakeAngle); 
 
      if(this->mBox.RakeAngle == 0)
	Face = rPoint[0] - this->mBox.Center[0] - this->mBox.Radius; 


      return Face;
    
      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************


    double& CalculateTipFace(double& Face, const TPointType& rPoint)
    {
      KRATOS_TRY
      
      Face = pow((rPoint[0] - this->mBox.Center[0]),2) + pow((rPoint[1] - this->mBox.Center[1]),2) - this->mBox.Radius * this->mBox.Radius; 
      
      return Face;	
    
      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************


    double& CalculateClearanceFace(double& Face, const TPointType& rPoint)
    {
      KRATOS_TRY
      
      Face = rPoint[1] - this->mBox.n_factor * ( rPoint[0] - this->mBox.Center[0] - this->mBox.Radius * sin(this->mBox.ClearanceAngle)) - this->mBox.Center[1]  +  this->mBox.Radius * cos(this->mBox.ClearanceAngle); 
 
      return Face;
 

    
      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************


    void CalculateAuxiliarFaces(double& rFace1, double& rFace2,  double& rFace3, const TPointType& rPoint)
    {
      KRATOS_TRY

	double pi = 3.141592654;

	rFace1 = rPoint[1] - tan(pi*0.25) * ( rPoint[0] - this->mBox.Center[0] ) - this->mBox.Center[1];

	rFace2 = rPoint[0]  + this->mBox.n_factor * ( rPoint[1] - this->mBox.Center[1] ) - this->mBox.Center[0];

	rFace3 = rPoint[1] + tan(this->mBox.RakeAngle) * ( rPoint[0] - this->mBox.Center[0] ) - this->mBox.Center[1];

      if(this->mBox.RakeAngle==0)
      	rFace3 = rPoint[1] - this->mBox.Center[1];

      KRATOS_CATCH( "" )
	}


    //************************************************************************************
    //************************************************************************************


    bool CalculateRakeSurface(const TPointType& rPoint, double& rGapNormal, double& rGapTangent, TPointType& rNormal, TPointType& rTangent)
    {
      KRATOS_TRY
      
      //1.-compute contact normal
      rNormal[0] = -cos(this->mBox.RakeAngle);
      rNormal[1] =  sin(this->mBox.RakeAngle);
      rNormal[2] = 0;

      //2.-compute point projection
      TPointType Projection;

      Projection[0] = rPoint[0] + this->mBox.m_factor * ( rPoint[1] - this->mBox.m_factor * ( this->mBox.Radius * cos(this->mBox.RakeAngle) - this->mBox.Center[0] ) - this->mBox.Center[1] ) - this->mBox.Radius * sin(this->mBox.RakeAngle);

      Projection[1] = this->mBox.m_factor * ( this->mBox.m_factor * rPoint[1] + rPoint[0] + this->mBox.Radius * cos(this->mBox.RakeAngle) - this->mBox.Center[0] ) + this->mBox.Center[1] + this->mBox.Radius * sin(this->mBox.RakeAngle);
      
      Projection[2] = 0;
      
      Projection /= (1+this->mBox.m_factor * this->mBox.m_factor);


      if(this->mBox.RakeAngle == 0){
	Projection[0] = this->mBox.Center[0] - this->mBox.Radius; 
	Projection[1] = rPoint[1]; 
      }

      //3.-compute gap

      rGapNormal = inner_prod((rPoint - Projection), rNormal);
      

      if(rGapNormal<0)
	return true;
      else
	return false;


      KRATOS_CATCH( "" )
	}

    //************************************************************************************
    //************************************************************************************

    bool CalculateTipSurface(const TPointType& rPoint, double& rGapNormal, double& rGapTangent, TPointType& rNormal, TPointType& rTangent)
    {
      KRATOS_TRY

 
      //1.-compute point projection
      TPointType Projection;
      Projection = this->mBox.Radius * ( (rPoint-this->mBox.Center)/ norm_2(rPoint-this->mBox.Center) ) + this->mBox.Center;

      
      //2.-compute contact normal
      rNormal = (Projection-this->mBox.Center)/this->mBox.Radius;


      //3.-compute gap

      if( norm_2(this->mBox.Center-rPoint) <= this->mBox.Radius ){
	rGapNormal = (-1) * norm_2(rPoint - Projection);
      }
      else{
	rGapNormal = norm_2(Projection - rPoint);
      }
      

   
      if(rGapNormal<0)
	return true;
      else
	return false;

    
      KRATOS_CATCH( "" )
	}


    //************************************************************************************
    //************************************************************************************


    bool CalculateClearanceSurface(const TPointType& rPoint, double& rGapNormal, double& rGapTangent, TPointType& rNormal, TPointType& rTangent)
    {
      KRATOS_TRY

      //1.-compute contact normal
      rNormal[0] =  sin(this->mBox.ClearanceAngle);
      rNormal[1] = -cos(this->mBox.ClearanceAngle);
      rNormal[2] = 0;

      //2.-compute point projection
      TPointType Projection;

      Projection[0] = rPoint[0] + this->mBox.n_factor * ( rPoint[1] + this->mBox.n_factor * ( this->mBox.Center[0] + this->mBox.Radius * sin(this->mBox.ClearanceAngle) ) - this->mBox.Center[1] + this->mBox.Radius * cos(this->mBox.ClearanceAngle) );

      Projection[1] = this->mBox.n_factor * ( this->mBox.n_factor * rPoint[1] + rPoint[0] + this->mBox.Center[0] - this->mBox.Radius * sin(this->mBox.ClearanceAngle) ) + this->mBox.Center[1] - this->mBox.Radius * cos(this->mBox.ClearanceAngle);
      
      Projection[2] = 0;
      
      Projection /= (1+this->mBox.n_factor * this->mBox.n_factor);


      //3.-compute gap

      rGapNormal = inner_prod((rPoint - Projection), rNormal);
      
      if(rGapNormal<0)
	return true;
      else
	return false;
     
      KRATOS_CATCH( "" )
     }


    //************************************************************************************
    //************************************************************************************

   static inline double inner_prod(const array_1d<double, 3>& a, const array_1d<double, 3>& b)
    {
        double temp =a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
        return temp;
    }

    //************************************************************************************
    //************************************************************************************

    static inline double norm_2(const array_1d<double, 3>& a)
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


}; // Class RigidToolBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TPointType, class TPointerType>
inline std::istream& operator >> (std::istream& rIStream,
                                  RigidToolBoundingBox& rThis);

/// output stream function
template<class TPointType, class TPointerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RigidToolBoundingBox& rThis)
{
    // rThis.PrintInfo(rOStream);
    // rOStream << std::endl;
    // rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_RIGID_TOOL_BOUNDING_BOX_H_INCLUDED  defined 


