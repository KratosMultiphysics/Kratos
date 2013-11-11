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

  typedef Vector              TPointType;
  
  
protected:

  typedef struct
  {
    double Radius;         //tool tip radius
    double RakeAngle;      //top angle,    from vertical axis
    double ClearanceAngle; //bottom angle, from the horizontal axis
      
    double m_factor;       //cotan(RakeAngle) == 1/tan(RakeAngle) == tan((pi/2)-RakeAngle)
    double n_factor;       //tan(ClearanceAngle)
    
    TPointType  HighPoint;      // box highest point
    TPointType  LowPoint;       // box lowest point

    TPointType  OriginalCenter; // center original position
    TPointType  Center;         // center current position
    TPointType  Velocity;       // velocity, expressed on the center velocity
    
  public:
    
    void clear()
    {
      Radius = 0;
      RakeAngle = 0;
      ClearanceAngle = 0;
      
      m_factor = 0;
      n_factor = 0;

      HighPoint.clear();
      LowPoint.clear();

      OriginalCenter.clear();
      Center.clear();
      Velocity.clear();
    }


  } BoundingBoxVariables;



public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION(SpatialBoundingBox);



    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpatialBoundingBox()
    {
      mBox.clear();
      std::cout<< "Calling empty constructor" <<std::endl;
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
      mBox.Velocity = rVelocity;

      TPointType Side(rCenter.size());
      Side[0] = 1.8 * mBox.Radius;
      Side[1] = mBox.Radius;

      if(rCenter.size()>2)
	Side[2] = 1.8 * mBox.Radius;

      mBox.HighPoint = rCenter + Side;
      mBox.LowPoint  = rCenter - Side;

    }


    SpatialBoundingBox( double Radius,
			double RakeAngle,
			double ClearanceAngle,
			TPointType  Center,
			TPointType  Velocity)
    {
      mBox.clear();
      mBox.Radius         = Radius;
      mBox.RakeAngle      = RakeAngle;
      mBox.ClearanceAngle = ClearanceAngle; 
      mBox.Center         = Center;
      mBox.Velocity       = Velocity;
    }


   SpatialBoundingBox(ModelPart &rModelPart,const double& rRadius)
    {
      
      double max=9.999999999999999e300;
      //double min=1e-300;
      
      ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();
      const unsigned int dimension = element_begin->GetGeometry().WorkingSpaceDimension();

      TPointType Maximum(dimension,-max);
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
      mBox.Radius = rRadius + 0.5*(Maximum[0]-Minimum[0]);

      mBox.Velocity = ZeroVector(3);

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

    virtual bool IsInside (const TPointType& rPoint, double& rCurrentTime)
    {
      bool inside = true;

      TPointType Reference = mBox.Center + mBox.Velocity * rCurrentTime;
      
      if(norm_2((Reference-rPoint)) > 2 * mBox.Radius)
	inside = false;

      Reference = mBox.HighPoint + mBox.Velocity * rCurrentTime;

      for(unsigned int i=0; i<mBox.Center.size(); i++)
	{
	  if(Reference[i]<rPoint[i]){
	    inside = false;
	    break;
	  }
	}

      Reference = mBox.LowPoint + mBox.Velocity * rCurrentTime;

      for(unsigned int i=0; i<mBox.Center.size(); i++)
	{
	  if(Reference[i]>rPoint[i]){
	    inside = false;
	    break;
	  }
	}

           
      return inside;
    }



    virtual bool IsInside(const TPointType& rPoint, double& rGapNormal, double& rGapTangent, TPointType& rNormal, TPointType& rTangent, int ContactFace = 0)
    {
      std::cout<< "Calling empty method" <<std::endl;
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

    TPointType& Center()
    {
        return mBox.Center;
    }

    TPointType& OriginalCenter()
    {
        return mBox.OriginalCenter;
    }

    double  Radius()
    {
        return mBox.Radius;
    }

    TPointType& Velocity()
    {
        return mBox.Velocity;
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

      TPointType Reference = mBox.HighPoint + mBox.Velocity * rCurrentTime;
      
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

	Reference = mBox.LowPoint + mBox.Velocity * rCurrentTime;
	
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


