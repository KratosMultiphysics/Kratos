/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
/* *********************************************************   
*          
*   Last Modified by:    $Author: Nelson $
*   Date:                $Date: 2011-03-29 11:41:31 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_SEGMENT_2D_INCLUDED )
#define  KRATOS_SEGMENT_2D_INCLUDED
//System includes
//External includes
#include "boost/smart_ptr.hpp"
#include <cmath>

//Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
      class Segment2D
	{
	  
	  public:
	    
	  typedef array_1d<double,2 >  PointType;  
	  
	  Segment2D(){}
	  Segment2D(const PointType& P0, const PointType& P1)
	  {
	  noalias(mP0) = P0;
	  noalias(mP1) = P1;
	  ComputeCenterDirectionExtent();
	  }

	  ~Segment2D(){}

	  double mExtent;
	  PointType mP0;
	  PointType mP1;
	  PointType mCenter;
	  PointType mDirection;

	  
	  // Assignment operator.
	  Segment2D& operator=(const Segment2D& rOther)
	  {
	    mExtent    = rOther.mExtent;
	    mP0        = rOther.mP0;
	    mP1        = rOther.mP1;
	    mCenter    = rOther.mCenter;
	    mDirection = rOther.mDirection; 
	    return *this;
	  }

	  // Copy constructor.
	  Segment2D(const Segment2D& rOther)
	  {
	    *this =  rOther;
	  }
	  
	   
	  //--------------------------------------------------------------------------->
	  //--------------------------------------------------------------------------->
	  
	  void AssignPointsAndComputeParameters(const PointType& P0, const PointType& P1)
	  {
	    noalias(mP0) = P0;
	    noalias(mP1) = P1;
	    ComputeCenterDirectionExtent();
	  }
	  
	  //--------------------------------------------------------------------------->
	  //--------------------------------------------------------------------------->
	  
	  void ComputeCenterDirectionExtent ()
	  {
	    noalias(mCenter)    =  (0.50) * (mP0 + mP1);
	    noalias(mDirection) =  mP1 - mP0;
	    const double length =  (std::sqrt(inner_prod(mDirection, mDirection) ) );
	    mExtent             =  0.500 * length;
	    mDirection          =  (1.00/(length ) )* mDirection;
	  }

	  //--------------------------------------------------------------------------->
	  //--------------------------------------------------------------------------->
	  
	  inline double DotPerp (const PointType& vec)
	  {
	     return mDirection[0]*vec[1] - mDirection[1]*vec[0];
	  }
	  
	  //--------------------------------------------------------------------------->
	  //--------------------------------------------------------------------------->
	   
	  inline double DistPoint2Segment2D(const PointType& rPoint)
	  {
	      array_1d<double,2 > diff = rPoint - this->mCenter;
	      array_1d<double,2 > ClosestPoint1;
	      array_1d<double,2 > ClosestPoint0; 
	      double param = inner_prod(this->mDirection, diff);

	      if (-this->mExtent < param)
	      {
	         if (param < this->mExtent)
	          {
	            ClosestPoint1 = this->mCenter + param * this->mDirection;
	          }
	      else
	          {
	           ClosestPoint1 = this->mP1;
	          }
	      }
	      else
	      {
	          ClosestPoint1 = this->mP0;
	      }

	      ClosestPoint0 = rPoint;
	      diff          = ClosestPoint1 - ClosestPoint0;
	      return std::sqrt(inner_prod(diff, diff));
	      }
	      
         };
	 
///************************************************************************************************************
///************************************************************************************************************

      class IntersectionSegment2DToSegment2D
	 {
	   
	   public:
	     
	   typedef array_1d<double,2 >  PointType;  
	   typedef enum Intersect{IT_POINT   = 0, IT_SEGMENT, IT_EMPTY}  Intersect;
	    
	   
	   static Intersect IntersectSegment(
            PointType& Point, 
            vector<PointType >& Points0,
            vector<PointType >& Points1)
          {
	
	  KRATOS_TRY

	  double toler = 1E-14;
	  Point        = ZeroVector(2);
	  PointType parameter = ZeroVector(2);


	  Segment2D Segment0(Points0[0], Points0[1]);
	  Segment2D Segment1(Points1[0], Points1[1]);
	  Intersect IntersectionType = Classify(parameter, Segment0,  Segment1);
		
	  if (IntersectionType == IT_POINT)
	  {
	    // Test whether the line-line intersection is on the segments.
	    double a = std::fabs(parameter[0]) - Segment0.mExtent;
	    double b = std::fabs(parameter[1]) - Segment1.mExtent;
	    
	    if ( a<=toler && b<=toler)
	    {
	      Point    = Segment0.mCenter + parameter[0]*Segment0.mDirection;
	      //comprobando que el punto no sea el extremo del segmento
	      array_1d<double, 4 > aa; 
	      aa[0] = std::fabs(Points0(0)[0]- Point[0]);
	      aa[1] = std::fabs(Points0(0)[1]- Point[1]);
	      aa[2] = std::fabs(Points1(1)[0]- Point[0]);
	      aa[3] = std::fabs(Points1(1)[1]- Point[1]);

	      if( (aa[0]<toler && aa[1]<toler) || (aa[2]<toler && aa[3]<toler)){
		IntersectionType = IT_EMPTY;
	      }
	    }
	    else
	      IntersectionType = IT_EMPTY;
	  }


	  return IntersectionType;  // != IT_EMPTY;

	  KRATOS_CATCH("")

	  }
          

	  //--------------------------------------------------------------------------->
	  //--------------------------------------------------------------------------->
          
      static Intersect Classify(
      PointType& s,
      Segment2D& Segment0,
      Segment2D& Segment1
      )
      
      {
	
	KRATOS_TRY
     
	// The intersection of two lines is a solution to P0+s0*D0 = P1+s1*D1.
	// Rewrite this as s0*D0 - s1*D1 = P1 - P0 = Q.  If D0.Dot(Perp(D1)) = 0,
	// the lines are parallel.  Additionally, if Q.Dot(Perp(D1)) = 0, the
	// lines are the same.  If D0.Dot(Perp(D1)) is not zero, then
	// s0 = Q.Dot(Perp(D1))/D0.Dot(Perp(D1))
	// produces the point of intersection.  Also,
	// s1 = Q.Dot(Perp(D0))/D0.Dot(Perp(D1))

        double toler                   = 1E-12;   
	PointType originDiff = Segment1.mCenter - Segment0.mCenter;
	PointType diff       = originDiff;
	double D0DotPerpD1             = Segment0.DotPerp(Segment1.mDirection);

	if ( std::fabs(D0DotPerpD1) > toler)
	{
	  // Lines intersect in a single point.
	  double invD0DotPerpD1 = 1.00/D0DotPerpD1;
	  double diffDotPerpD0  = originDiff[0]*Segment0.mDirection[1] - originDiff[1]*Segment0.mDirection[0]; 
	  double diffDotPerpD1  = originDiff[0]*Segment1.mDirection[1] - originDiff[1]*Segment1.mDirection[0]; 
	  s[0] = diffDotPerpD1*invD0DotPerpD1;
	  s[1] = diffDotPerpD0*invD0DotPerpD1;

	  return IT_POINT;
	}


	// Lines are parallel.
	originDiff = originDiff * (1.00 / ( std::sqrt(inner_prod(originDiff, originDiff) ) ) );    
	PointType diffN = originDiff;
	

	double diffNDotPerpD1 = originDiff[0]*Segment1.mDirection[1] - originDiff[1]*Segment1.mDirection[0]; 
	if (std::fabs(diffNDotPerpD1) <= toler)
	{
	// Lines are colinear
	return IT_SEGMENT;
	}


	// Lines are parallel, but distinct.
	return IT_EMPTY;
	
	KRATOS_CATCH("")
	
        }
	     
    };
	 
	 
	 
}//namespace Kratos.

#endif /* KRATOS_SEGMENT_2D_INCLUDED  defined */

