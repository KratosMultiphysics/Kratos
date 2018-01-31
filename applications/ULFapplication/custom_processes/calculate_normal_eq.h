//
//   Project Name:        Kratos
//   Last Modified by:    $Author: ajarauta $
//   Date:                $Date: 2007-11-06 12:34:26 $
//   Revision:            $Revision: 1.4 $
//
//  this process save structural elements in a separate list

#if !defined(CALCULATE_NORMAL_EQ_INCLUDED )
#define  CALCULATE_NORMAL_EQ_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "processes/process.h"
#include "includes/condition.h"
#include "includes/element.h"
#include "ULF_application.h"



namespace Kratos
{
  class CalculateNormalEq
  : public Process
  {
  public:
    void CalculateNormalEq3D(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	
	double pi = 3.14159265359;
	double theta_eq = ThisModelPart.GetProcessInfo()[CONTACT_ANGLE_STATIC];
	double theta_rad = theta_eq*pi/180.0;
	double theta = 0.0;
	
	array_1d<double,3> normal_eq = ZeroVector(3);
	array_1d<double,3> normal_xy = ZeroVector(3);
	normal_xy[2] = 1.0;
	array_1d<double,2> temp = ZeroVector(2);
	
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; im++)
	{
	    //Find the neighbours of TRIPLE_POINT at the boundary
	    if (im->FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	    {
		temp[0] = im->FastGetSolutionStepValue(NORMAL_TP_X);
		temp[1] = im->FastGetSolutionStepValue(NORMAL_TP_Y);
		temp = NormalizeVec2D(temp);
		normal_eq[0] = sin(theta_rad)*temp[0];
		normal_eq[1] = sin(theta_rad)*temp[1];
		normal_eq[2] = cos(theta_rad);
		normal_eq = NormalizeVec3D(normal_eq);	      
		im->FastGetSolutionStepValue(NORMAL_EQ) = normal_eq;
		
		im->FastGetSolutionStepValue(NORMAL_CL_EQ_X) = -cos(theta_rad)*temp[0];
		im->FastGetSolutionStepValue(NORMAL_CL_EQ_Y) = -cos(theta_rad)*temp[1];		
		im->FastGetSolutionStepValue(NORMAL_CL_EQ_Z) = sin(theta_rad);
		im->FastGetSolutionStepValue(NORMAL_CL_EQ) = NormalizeVec3D(im->FastGetSolutionStepValue(NORMAL_CL_EQ));

	    }
	    else
		im->FastGetSolutionStepValue(NORMAL_EQ) = ZeroVector(3);
	}	
	KRATOS_CATCH("")
    }    
    
   
    
    array_1d<double,2> Vector2D(const double x0, const double y0, const double x1, const double y1)
    {
      array_1d<double,2> r01;
      r01[0] = x1 - x0;
      r01[1] = y1 - y0;
      return r01;
    }    
    
    array_1d<double,3> Vector3D(const double x0, const double y0, const double z0,
		  const double x1, const double y1, const double z1)
    {
      array_1d<double,3> r01;
      r01[0] = x1 - x0;
      r01[1] = y1 - y0;
      r01[2] = z1 - z0;
      return r01;
    }
    
    double Norm2D(const array_1d<double,2>& a)
    {
      return sqrt(a[0]*a[0] + a[1]*a[1]);
    }
    
    double Norm3D(const array_1d<double,3>& a)
    {
      return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    }
    
    array_1d<double,2> NormalizeVec2D(array_1d<double,2>& input)
    {
      double norm = Norm2D(input);
      if (norm != 0.0)
      {
	input[0] /= norm;
	input[1] /= norm;
      }
      return input;
    }       
    
    array_1d<double,3> NormalizeVec3D(array_1d<double,3>& input)
    {
      double norm = Norm3D(input);
      if (norm != 0.0)
      {
	input[0] /= norm;
	input[1] /= norm;
	input[2] /= norm;
      }
      return input;
    }    
    
   array_1d<double,3> SumVecs3D(const array_1d<double,3>& a, const array_1d<double,3>& b)
    {
      array_1d<double,3> c;
      for (unsigned int i = 0; i < 3; i++)
      {
	c[i] = a[i] + b[i];
      }
      return c;
    }  
    
    double Angle2vecs3D(const array_1d<double,3>& a, const array_1d<double,3>& b)
    {
      double pi = 3.14159265359;
      double norm_a = Norm3D(a);
      double norm_b = Norm3D(b);
      double temp = 0.0;
      if (norm_a*norm_b > -0.00000001 && norm_a*norm_b < 0.00000001)
	temp = 0.0;
      else
	temp = DotProduct3D(a,b)/(norm_a*norm_b);
      return (acos(temp)*180.0/pi);
    }    
    
    double DotProduct3D(const array_1d<double,3>& a, const array_1d<double,3>& b)
    {
      return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    }    
    
    private:
    
  }; // Class CalculateNormalEq


}  // namespace Kratos.

#endif // KRATOS_CALCULATE_NORMAL_EQ_INCLUDED  defined 


