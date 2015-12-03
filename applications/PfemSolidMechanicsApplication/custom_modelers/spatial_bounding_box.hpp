//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
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
#include "includes/model_part.h"


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

class SpatialBoundingBox
{
public:

  typedef Vector                                         TPointType;
  typedef ModelPart::NodesContainerType          NodesContainerType;
  typedef NodesContainerType::Pointer     NodesContainerTypePointer;
 
protected:

  typedef struct
  {
    int    Dimension;           //2D or 3D
    bool   Axisymmetric;        //true or false
    int    Convexity;           //1 or -1  if "in" is inside or outside respectively   
    double Radius;              // box radius
    
    TPointType  HighPoint;      // box highest point
    TPointType  LowPoint;       // box lowest point

    TPointType  OriginalCenter; // center original position
    TPointType  Center;         // center current position
    
  public:
    
    void clear()
    {
      Dimension = 2;
      Axisymmetric = false;
      Convexity = 1;
      Radius = 0;

      HighPoint.clear();
      LowPoint.clear();

      OriginalCenter.clear();
      Center.clear();
    }


  } BoundingBoxVariables;



  typedef struct
  {

    int Label;

    TPointType  Velocity;         // velocity, relative to rotation center and box center (the same)
    TPointType  AngularVelocity;  // angular velocity, relative to the Rotation Center 

    TPointType  OriginalRotationCenter;   // original position for the rotation center
    TPointType  RotationCenter;           // rotation center


  public:
    
    void clear()
    {
      Label = 0;
      Velocity.clear();
      AngularVelocity.clear();
      RotationCenter.clear();
      OriginalRotationCenter.clear();
    }

  } BoxMovementVariables;


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
      mBox.clear();
      //std::cout<< " Calling Bounding Box empty constructor" <<std::endl;
    }

    SpatialBoundingBox(const TPointType& rPoint)
    {
      mBox.clear();
      mBox.HighPoint = rPoint;
      mBox.LowPoint  = rPoint;
    }

    SpatialBoundingBox(const TPointType& rLowPoint, const TPointType& rHighPoint )
    {

      mBox.clear();
      mBox.HighPoint = rHighPoint;
      mBox.LowPoint  = rLowPoint;

      mBox.Center = 0.5 *( rHighPoint + rLowPoint );

      mBox.Radius = 0.5 * norm_2(rHighPoint-rLowPoint);   
    }

    SpatialBoundingBox(const TPointType& rCenter, const double& rRadius)
    {

      mBox.clear();     
      mBox.Center = rCenter;
      mBox.Radius = rRadius;

      TPointType Side;
      Side[0] = mBox.Radius;
      Side[1] = mBox.Radius;

      if(rCenter.size()>2)
	Side[2] = mBox.Radius;

      mBox.HighPoint = rCenter + Side;
      mBox.LowPoint  = rCenter - Side;
    }

    SpatialBoundingBox(const TPointType& rCenter, const double& rRadius, const TPointType& rVelocity)
    {
      mBox.clear();     
      mBox.Center   = rCenter;
      mBox.Radius   = rRadius;
      mMovement.Velocity = rVelocity;

      TPointType Side(rCenter.size());
      Side[0] = 1.8 * mBox.Radius;
      Side[1] = mBox.Radius;

      if(rCenter.size()>2)
	Side[2] = 1.8 * mBox.Radius;

      mBox.HighPoint = rCenter + Side;
      mBox.LowPoint  = rCenter - Side;

    }


    SpatialBoundingBox(  int Label,
			 int Convexity,
			 double Radius,
			 NodesContainerTypePointer GeneratrixPoints)
    {
      std::cout<<" Calling a Base Class Constructor: RIGID TUBE Bounding Box must be called "<<std::endl;
    }


    SpatialBoundingBox( Vector Convexities,
			Vector Radius,
			Vector RakeAngles,
			Vector ClearanceAngles,
			Matrix Centers,
			TPointType Velocity,
			TPointType AngularVelocity,
			TPointType RotationCenter)
    {
      std::cout<<" Calling a Base Class Constructor: RIGID WALL Bounding Box must be called "<<std::endl;
    }


   SpatialBoundingBox(ModelPart &rModelPart,const double& rRadius)
    {
      
      double max=std::numeric_limits<double>::max();
      double min=std::numeric_limits<double>::min();
      
      ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();
      const unsigned int dimension = element_begin->GetGeometry().WorkingSpaceDimension();

      TPointType Maximum(dimension, min);
      TPointType Minimum(dimension, max);

      //Get inside point of the subdomains

      for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin(); in!=rModelPart.NodesEnd(); in++){
	if(in->Is(BOUNDARY) ){
	  
	  //get maximum
	  if(Maximum[0]<in->X())
	     Maximum[0]=in->X();

	  if(Maximum[1]<in->Y())
	     Maximum[1]=in->Y();

	  //get minimum
	  if(Minimum[0]>in->X())
	     Minimum[0]=in->X();

	  if(Minimum[1]>in->Y())
	     Minimum[1]=in->Y();

	  if(dimension>2){

	    if(Maximum[3]<in->Z())
	      Maximum[3]=in->Z();
	    
	    if(Minimum[3]>in->Z())
	      Minimum[3]=in->Z();
	  }
	     
	} 
      }
    
      mBox.Center = 0.5*(Maximum+Minimum);

      double MaxRadius = Maximum[0]-Minimum[0];
      if(Maximum[1]-Minimum[1]>MaxRadius)
	MaxRadius = Maximum[1]-Minimum[1];
      
      if(dimension>2){
	if(Maximum[2]-Minimum[2]>MaxRadius)
	  MaxRadius = Maximum[2]-Minimum[2];
      }
	  
      mBox.Radius = rRadius + 0.5*(MaxRadius);

      mMovement.Velocity = ZeroVector(3);

      TPointType Side(dimension);
      Side[0] = mBox.Radius;
      Side[1] = mBox.Radius;

      if(dimension>2)
	Side[2] = mBox.Radius;

      mBox.HighPoint = mBox.Center + Side;
      mBox.LowPoint  = mBox.Center - Side;

    }


    /// Assignment operator.
    virtual SpatialBoundingBox& operator=(SpatialBoundingBox const& rOther)
    {
        mBox = rOther.mBox;
        return *this;
    }


    /// Copy constructor.
    SpatialBoundingBox(SpatialBoundingBox const& rOther) 
    :mBox(rOther.mBox)
    {
    }


    /// Destructor.
    virtual ~SpatialBoundingBox() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual bool IsInside (const TPointType& rPoint, double& rCurrentTime, double Radius = 0)
    {
      bool inside = true;

      TPointType Reference = mBox.Center + mMovement.Velocity * rCurrentTime;
      
      if(norm_2((Reference-rPoint)) > 2 * mBox.Radius)
	inside = false;

      Reference = mBox.HighPoint + mMovement.Velocity * rCurrentTime;

      for(unsigned int i=0; i<mBox.Center.size(); i++)
	{
	  if(Reference[i]<rPoint[i]){
	    inside = false;
	    break;
	  }
	}

      Reference = mBox.LowPoint + mMovement.Velocity * rCurrentTime;

      for(unsigned int i=0; i<mBox.Center.size(); i++)
	{
	  if(Reference[i]>rPoint[i]){
	    inside = false;
	    break;
	  }
	}

           
      return inside;
    }


    //************************************************************************************
    //************************************************************************************

    virtual bool IsInside (const TPointType& rPoint, double& rCurrentTime, int& ContactFace, double Radius = 0)
    {
      ContactFace = 0;

      return IsInside(rPoint,rCurrentTime,Radius);
    }


    //************************************************************************************
    //************************************************************************************

    virtual bool IsInside(const TPointType& rPoint, double& rGapNormal, double& rGapTangent, TPointType& rNormal, TPointType& rTangent, double Radius = 0)
    {
      std::cout<< "Calling empty method" <<std::endl;
      return false;
    }


    //************************************************************************************
    //************************************************************************************
    virtual bool IsInside(const TPointType& rPoint, double& rGapNormal, double& rGapTangent, TPointType& rNormal, TPointType& rTangent, int& ContactFace, double Radius = 0)
    {
      std::cout<< "Calling empty method" <<std::endl;
      ContactFace = 0;
      return false;
    }

    ///@}
    ///@name Access
    ///@{

    void Set(const TPointType& rLowPoint, const TPointType& rHighPoint)
    {
        mBox.LowPoint  = rLowPoint;
        mBox.HighPoint = rHighPoint;
    }

    TPointType const& HighPoint() const
    {
        return mBox.HighPoint;
    }


    TPointType& HighPoint()
    {
        return mBox.HighPoint;
    }

    TPointType const& LowPoint() const
    {
        return mBox.LowPoint;
    }


    TPointType& LowPoint()
    {
        return mBox.LowPoint;
    }

    const TPointType& Center()
    {
        return mBox.Center;
    }


    TPointType& OriginalCenter()
    {
        return mBox.OriginalCenter;
    }

    const double& Radius()
    {
        return mBox.Radius;
    }

    virtual TPointType GetCenter()
    {
        return mBox.Center;
    }

    virtual TPointType GetCenter(const TPointType& rPoint)
    {
       return this->GetCenter();
    }


    virtual double GetRadius()
    {
        return mBox.Radius;
    }

    virtual double GetRadius(const TPointType& rPoint)
    {
       return this->GetRadius();
    }

    virtual void SetRadius(double rRadius)
    {
        mBox.Radius = rRadius;
    }

    TPointType& Velocity()
    {
        return mMovement.Velocity;
    }


    int GetMovementLabel()
    {
        return mMovement.Label;
    }

    void SetMovementLabel(int &rLabel)
    {
        mMovement.Label = rLabel;
    }


    void SetDimension(int dimension)
    {
        mBox.Dimension = dimension;
    }

    int GetDimension()
    {
        return mBox.Dimension;
    }


    void SetAxisymmetric()
    {
        mBox.Axisymmetric = true;
    }

    bool& Axisymmetric()
    {
        return mBox.Axisymmetric;
    }

  
    virtual void UpdatePosition(double & rTime)
    {
      
      mBox.Center = mBox.OriginalCenter + mMovement.Velocity * rTime;
      
    }

    /// Compute inside holes
    std::vector<TPointType > GetHoles(ModelPart &rModelPart)
    {
      //Get inside point of the subdomains
      ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();
      const unsigned int dimension = element_begin->GetGeometry().WorkingSpaceDimension();

      unsigned int start=0;
      unsigned int NumberOfMeshes=rModelPart.NumberOfMeshes();
      if(NumberOfMeshes>1) 
	start=1;

      std::vector<TPointType > Holes;
      for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	{
	  TPointType Point(dimension);
	  for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin(MeshId); in!=rModelPart.NodesEnd(MeshId); in++){
	    if(in->IsNot(BOUNDARY) ){
	      Point[0] = in->X();	
	      Point[1] = in->Y();

	      if(dimension>2)
		Point[2] = in->Z();		

	      Holes.push_back(Point);
	      break;
	    }

	  }
	}

      return Holes;
    }


    /// Compute vertices
    std::vector<TPointType > GetVertices(double& rCurrentTime)
    {
    
      std::vector<TPointType> vertices;

      TPointType Reference = mBox.HighPoint + mMovement.Velocity * rCurrentTime;
      
      double Side =  2.0 * (mBox.HighPoint[0] - mBox.Center[0]);

      //point 1
      vertices.push_back(Reference);
      
      Reference[0] -= Side;
      
      //point 2
      vertices.push_back(Reference);
      
      Reference[1] -= Side;
      
      //point 3
      vertices.push_back(Reference);
      
      Reference[0] += Side;
      
      //point 4
      vertices.push_back(Reference);
      

      if( mBox.Center.size() > 2){

	Reference = mBox.LowPoint + mMovement.Velocity * rCurrentTime;
	
	//point 5
	vertices.push_back(Reference);
      
	Reference[0] += Side;
      
	//point 6
	vertices.push_back(Reference);
      
	Reference[1] += Side;
      
	//point 7
	vertices.push_back(Reference);
      
	Reference[0] -= Side;
      
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
        rOStream << mBox.HighPoint << " , " << mBox.LowPoint;
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

    BoundingBoxVariables mBox;

    BoxMovementVariables mMovement;

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


/// input stream function
template<class TPointType, class TPointerType>
inline std::istream& operator >> (std::istream& rIStream,
                                  SpatialBoundingBox& rThis);

/// output stream function
template<class TPointType, class TPointerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SpatialBoundingBox& rThis)
{
    // rThis.PrintInfo(rOStream);
    // rOStream << std::endl;
    // rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SPATIAL_BOUNDING_BOX_H_INCLUDED  defined 


