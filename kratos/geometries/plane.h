/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2009-01-21 09:56:10 $
//   Revision:            $Revision: 1.14 $
//
//
#if !defined(KRATOS_PLANE_H_INCLUDED )
#define  KRATOS_PLANE_H_INCLUDED 



// System includes
#include <iostream>

// External includes
#include <boost/array.hpp> 
#include "utilities/math_utils.h"



namespace Kratos
{
 
    // The plane is represented as Dot(N,X) = c where N is a unit-length
    // normal vector, c is the plane constant, and X is any point on the
    // plane.  The user must ensure that the normal vector is unit length.
    class Plane
    {
      public:
	  
      Plane()
      {
	 noalias(mTriangleBary) = ZeroVector(3);
 	 noalias(mClosestPoint) = ZeroVector(3);  
         noalias(mNormal)       = ZeroVector(3);
	 mConstant = 0.00;
      }
      
     ~Plane(){}  
      
      //----------------------------------------------------------------------------
      Plane(const array_1d<double, 3>& normal, double& constant):
        mNormal(normal),
        mConstant(constant)
        {
        }
     //----------------------------------------------------------------------------
      Plane(const array_1d<double, 3>& normal, const array_1d<double, 3>& p):
       mNormal(normal)
      {
        mConstant = inner_prod(normal, p);
      }
     //----------------------------------------------------------------------------
     Plane(const array_1d<double, 3>& p0, 
	   const array_1d<double, 3>& p1,
           const array_1d<double, 3>& p2)
       {
        array_1d<double, 3> edge1 = p1 - p0;
        array_1d<double, 3> edge2 = p2 - p0;
        mNormal   =  MathUtils<double>::UnitCrossProduct(edge1,edge2);
        mConstant =  inner_prod(mNormal, p0);
      }

      
      void AssignPointsAndComputeParameters(array_1d<double, 3>& p0, array_1d<double, 3>& p1, array_1d<double, 3>& p2)
      {
	array_1d<double, 3> edge1 = p1 - p0;
	array_1d<double, 3> edge2 = p2 - p0;
	mNormal   =  MathUtils<double>::UnitCrossProduct(edge1,edge2);
	mConstant =  inner_prod(mNormal, p0);
      }

      //----------------------------------------------------------------------------
      /// Calcula la distancia del punto al plano  
      double DistanceTo(const array_1d<double, 3>& p) 
      {
         return  inner_prod(mNormal,p) - mConstant;
      }
      
      /// computa la distancia de un punto a un triangulo 3D
     double DistPoint3Triangle3(
                  array_1d<double, 3>& rPoint,
                  array_1d<double, 3>& p0, 
                  array_1d<double, 3>& p1, 
	          array_1d<double, 3>& p2 
     )
    {
	    array_1d<double, 3> diff  = p0 - rPoint;
	    array_1d<double, 3> edge0 = p1 - p0;
	    array_1d<double, 3> edge1 = p2 - p0;
	    
	    double a00 = inner_prod(edge0, edge0); 
	    double a01 = inner_prod(edge0, edge1);
	    double a11 = inner_prod(edge1, edge1); 
	    double b0  = inner_prod(diff,  edge0);
	    double b1  = inner_prod(diff,  edge1);
	    double c   = inner_prod(diff,  diff); 
	    double det = std::fabs(a00*a11 - a01*a01);
	    double s = a01*b1 - a11*b0;
	    double t = a01*b0 - a00*b1;
	    double sqrDistance = 0.00;

	    if (s + t <= det)
	    {
	    if (s < 0.00)
	    {
		if (t < 0.00 )  // region 4
		{
		    if (b0 < 0.00)
		    {
			t = 0.00;
			if (-b0 >= a00)
			{
			    s = 1.00;
			    sqrDistance = a00 + (2.00)*b0 + c;
			}
			else
			{
			    s = -b0/a00;
			    sqrDistance = b0*s + c;
			}
		    }
		    else
		    {
			s = 0.00;
			if (b1 >= 0.00)
			{
			    t = 0.00;
			    sqrDistance = c;
			}
			else if (-b1 >= a11)
			{
			    t = 1.00;
			    sqrDistance = a11 + (2.00)*b1 + c;
			}
			else
			{
			    t = -b1/a11;
			    sqrDistance = b1*t + c;
			}
		    }
		}
		else  // region 3
		{
		    s = 0.00;
		    if (b1 >= 0.00)
		    {
			t = 0.00;
			sqrDistance = c;
		    }
		    else if (-b1 >= a11)
		    {
			t = 1.00;
			sqrDistance = a11 + (2.00)*b1 + c;
		    }
		    else
		    {
			t = -b1/a11;
			sqrDistance = b1*t + c;
		    }
		}
	    }
	    else if (t < 0.00)  // region 5
	    {
		t = 0.00;
		if (b0 >= 0.00)
		{
		    s = 0.00;
		    sqrDistance = c;
		}
		else if (-b0 >= a00)
		{
		    s = 1.00;
		    sqrDistance = a00 + (1.00)*b0 + c;
		}
		else
		{
		    s = -b0/a00;
		    sqrDistance = b0*s + c;
		}
	    }
	    else  // region 0
	    {
		// minimum at interior point
		double invDet = (1.00)/det;
		s *= invDet;
		t *= invDet;
		sqrDistance = s*(a00*s + a01*t + (2.00)*b0) +
		    t*(a01*s + a11*t + (2.00)*b1) + c;
	    }
	    }
	    else
	    {
	    double tmp0, tmp1, numer, denom;

	    if (s < 0.00)  // region 2
	    {
		tmp0 = a01 + b0;
		tmp1 = a11 + b1;
		if (tmp1 > tmp0)
		{
		    numer = tmp1 - tmp0;
		    denom = a00 - (2.00)*a01 + a11;
		    if (numer >= denom)
		    {
			s = 1.00;
			t = 0.00;
			sqrDistance = a00 + (2.00)*b0 + c;
		    }
		    else
		    {
			s = numer/denom;
			t = 1.00 - s;
			sqrDistance = s*(a00*s + a01*t + (2.00)*b0) +
			    t*(a01*s + a11*t + (2.00)*b1) + c;
		    }
		}
		else
		{
		    s = 0.00;
		    if (tmp1 <= 0.00)
		    {
			t = 1.00;
			sqrDistance = a11 + (2.00)*b1 + c;
		    }
		    else if (b1 >= 0.00)
		    {
			t = 0.00;
			sqrDistance = c;
		    }
		    else
		    {
			t = -b1/a11;
			sqrDistance = b1*t + c;
		    }
		}
	    }
	    else if (t < 0.00)  // region 6
	    {
		tmp0 = a01 + b1;
		tmp1 = a00 + b0;
		if (tmp1 > tmp0)
		{
		    numer = tmp1 - tmp0;
		    denom = a00 - (2.00)*a01 + a11;
		    if (numer >= denom)
		    {
			t = 1.00;
			s = 0.00;
			sqrDistance = a11 + (2.00)*b1 + c;
		    }
		    else
		    {
			t = numer/denom;
			s = 1.00 - t;
			sqrDistance = s*(a00*s + a01*t + (2.00)*b0) +
			    t*(a01*s + a11*t + (2.00)*b1) + c;
		    }
		}
		else
		{
		    t = 0.00;
		    if (tmp1 <= 0.00)
		    {
			s = 1.00;
			sqrDistance = a00 + (2.00)*b0 + c;
		    }
		    else if (b0 >= 0.00)
		    {
			s = 0.00;
			sqrDistance = c;
		    }
		    else
		    {
			s = -b0/a00;
			sqrDistance = b0*s + c;
		    }
		}
	    }
	    else  // region 1
	    {
		numer = a11 + b1 - a01 - b0;
		if (numer <= 0.00)
		{
		    s = 0.00;
		    t = 1.00;
		    sqrDistance = a11 + (2.00)*b1 + c;
		}
		else
		{
		    denom = a00 - (2.00)*a01 + a11;
		    if (numer >= denom)
		    {
			s = 1.00;
			t = 0.00;
			sqrDistance = a00 + (2.00)*b0 + c;
		    }
		    else
		    {
			s = numer/denom;
			t = 1.00 - s;
			sqrDistance = s*(a00*s + a01*t + (2.00)*b0) +
			    t*(a01*s + a11*t + (2.00)*b1) + c;
		    }
		}
	    }
	    }

	    // Account for numerical round-off error.
	    if (sqrDistance < 0.00)
	    {
	      sqrDistance = 0.00;
	    }

	    //mClosestPoint0 = *mPoint;
	    mClosestPoint = p0 + s*edge0 + t*edge1;
	    mTriangleBary[0] = (double)1 - s - t;
	    mTriangleBary[1] = s;
	    mTriangleBary[2] = t;
 	    
	    return sqrDistance;
}

      
      
      //----------------------------------------------------------------------------
      int WhichSide(const array_1d<double, 3>& p) 
      {
	  double distance = DistanceTo(p);
	  if (distance < 0.00)
	      return -1;
	  else if (distance > 0.00)
	      return +1;
	  else
	      return 0;
      }
 
 	 array_1d<double, 3> mTriangleBary;
 	 array_1d<double, 3> mClosestPoint; 
         array_1d<double, 3> mNormal;
	 double mConstant;
    };
}
#endif 

