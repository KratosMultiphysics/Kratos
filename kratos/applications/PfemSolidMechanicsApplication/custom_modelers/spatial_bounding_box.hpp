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
#include "pfem_solid_mechanics_application.h"

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
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION(SpatialBoundingBox);

    typedef Vector TPointType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpatialBoundingBox() : mHighPoint(), mLowPoint()
    {
        std::cout<< "Calling empty constructor" <<std::endl;
    }

    SpatialBoundingBox(const TPointType& Point) :  mHighPoint(Point), mLowPoint(Point)
    {
    }

    SpatialBoundingBox(const TPointType& LowPoint, const TPointType& HighPoint ) :  mHighPoint(HighPoint), mLowPoint(LowPoint)
    {
      mCenter = 0.5 *( HighPoint + LowPoint );
      mRadius = 0.5 * norm_2(HighPoint-LowPoint);   
    }

    SpatialBoundingBox(const TPointType& Center, const double& Radius) :  mCenter(Center), mRadius(Radius)
    {
      TPointType Side;
      Side[0] = mRadius;
      Side[1] = mRadius;

      if(Center.size()>2)
	Side[2] = mRadius;

      mHighPoint = Center + Side;
      mLowPoint  = Center - Side;
    }

    SpatialBoundingBox(const TPointType& Center, const double& Radius, const TPointType& Velocity) :  mCenter(Center), mRadius(Radius), mVelocity(Velocity)
    {
      TPointType Side(Center.size());
      Side[0] = 1.8 * mRadius;
      Side[1] = mRadius;

      if(Center.size()>2)
	Side[2] = 1.8 * mRadius;

      mHighPoint = Center + Side;
      mLowPoint  = Center - Side;

    }


   SpatialBoundingBox(ModelPart &rModelPart,const double& Radius)
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
    
      mCenter = 0.5*(Maximum+Minimum);
      mRadius = Radius + 0.5*(Maximum[0]-Minimum[0]);

      mVelocity = ZeroVector(3);

      TPointType Side(dimension);
      Side[0] = mRadius;
      Side[1] = mRadius;

      if(dimension>2)
	Side[2] = mRadius;

      mHighPoint = mCenter + Side;
      mLowPoint  = mCenter - Side;

    }

    /// Destructor.
    virtual ~SpatialBoundingBox() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    bool IsInside (const TPointType& Point, double & CurrentTime)
    {
      bool inside = true;

      TPointType Reference = mCenter + mVelocity * CurrentTime;
      
      if(norm_2((Reference-Point)) > 2 * mRadius)
	inside = false;

      Reference = mHighPoint + mVelocity * CurrentTime;

      for(unsigned int i=0; i<mCenter.size(); i++)
	{
	  if(Reference[i]<Point[i]){
	    inside = false;
	    break;
	  }
	}

      Reference = mLowPoint + mVelocity * CurrentTime;

      for(unsigned int i=0; i<mCenter.size(); i++)
	{
	  if(Reference[i]>Point[i]){
	    inside = false;
	    break;
	  }
	}

           
      return inside;
    }

    ///@}
    ///@name Access
    ///@{

    void Set(const TPointType& LowPoint, const TPointType& HighPoint)
    {
        mLowPoint  = LowPoint;
        mHighPoint = HighPoint;
    }

    TPointType const& HighPoint() const
    {
        return mHighPoint;
    }


    TPointType& HighPoint()
    {
        return mHighPoint;
    }

    TPointType const& LowPoint() const
    {
        return mLowPoint;
    }


    TPointType& LowPoint()
    {
        return mLowPoint;
    }

    TPointType& Center()
    {
        return mCenter;
    }

    double  Radius()
    {
        return mRadius;
    }

    TPointType& Velocity()
    {
        return mVelocity;
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
    std::vector<TPointType > GetVertices(double & CurrentTime)
    {
    
      std::vector<TPointType> vertices;

      TPointType Reference = mHighPoint + mVelocity * CurrentTime;
      
      double Side =  2.0 * (mHighPoint[0] - mCenter[0]);

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
      

      if( mCenter.size() > 2){

	Reference = mLowPoint + mVelocity * CurrentTime;
	
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


    /// Assignment operator.
    SpatialBoundingBox& operator=(SpatialBoundingBox const& rOther)
    {
        mHighPoint = rOther.mHighPoint;
        mLowPoint  = rOther.mLowPoint;
	mCenter    = rOther.mCenter;
	mRadius    = rOther.mRadius;
	mVelocity  = rOther.mVelocity;
        return *this;
    }


    /// Copy constructor.
    SpatialBoundingBox(SpatialBoundingBox const& rOther) :
        mHighPoint(rOther.mHighPoint),
        mLowPoint(rOther.mLowPoint),
	mCenter(rOther.mCenter),
	mRadius(rOther.mRadius),
	mVelocity(rOther.mVelocity)
    {
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
        rOStream << mHighPoint << " , " << mLowPoint;
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
    TPointType mHighPoint;
    TPointType mLowPoint;

    TPointType mCenter;
    double     mRadius;

    TPointType mVelocity;
  
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


