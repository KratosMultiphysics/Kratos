//
//   Project Name:        Kratos
//   Last Modified by:    $Author: ajarauta $
//   Date:                $Date: 2007-11-06 12:34:26 $
//   Revision:            $Revision: 1.4 $
//
//  this process save structural elements in a separate list

#if !defined(CALCULATE_ADHESION_FORCE_INCLUDED )
#define  CALCULATE_ADHESION_FORCE_INCLUDED



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
  class CalculateAdhesionForce
  : public Process
  {
  public:
    /*
    void CalculateAdhesionForce2D(ModelPart& ThisModelPart)
    {
	//KRATOS_TRY
	
	//KRATOS_CATCH("")
    }
    */
    void CalculateAdhesionForce3D(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	
	double gamma = ThisModelPart.GetProcessInfo()[SURFTENS_COEFF]; //Surface tension coefficient [N m-1]
	double pi = 3.14159265359;
	double theta_eq = ThisModelPart.GetProcessInfo()[CONTACT_ANGLE_STATIC];
	double theta, sign_force, theta_rad, cos_t;
	array_1d<double,2> nu0 = ZeroVector(2);
	array_1d<double,2> vel = ZeroVector(2);
	array_1d<double,2> n_tp = ZeroVector(2);
	
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ;
	    im != ThisModelPart.NodesEnd() ; im++)
	    {
	      //Find the neighbours of TRIPLE_POINT at the boundary
	      if ((im->FastGetSolutionStepValue(TRIPLE_POINT))*1000 != 0.0)
	      {
		nu0[0] = im->FastGetSolutionStepValue(NORMAL_TP_X);
		nu0[1] = im->FastGetSolutionStepValue(NORMAL_TP_Y);
		NormalizeVec2D(nu0);
		theta = im->FastGetSolutionStepValue(CONTACT_ANGLE);
		theta_rad = (theta)*pi/180.0;
		cos_t = cos(theta_rad);
		cos_t = sqrt(cos_t*cos_t);
		vel[0] = im->FastGetSolutionStepValue(VELOCITY_X);
		vel[1] = im->FastGetSolutionStepValue(VELOCITY_Y);	
		if (vel[0] != 0.0 && vel[1] != 0.0)
		{
// 		    sign_force = fsign(theta,theta_eq);
		    sign_force = DotProduct2D(nu0,vel);
		}
		else
		{
		    n_tp[0] = im->FastGetSolutionStepValue(NORMAL_TP_X);
		    n_tp[1] = im->FastGetSolutionStepValue(NORMAL_TP_Y);
		    sign_force = DotProduct2D(nu0,n_tp);
		}
		im->FastGetSolutionStepValue(ADHESION_FORCE_X) = -sign_force*gamma*cos_t*(im->FastGetSolutionStepValue(NODAL_LENGTH))*nu0[0];
		im->FastGetSolutionStepValue(ADHESION_FORCE_Y) = -sign_force*gamma*cos_t*(im->FastGetSolutionStepValue(NODAL_LENGTH))*nu0[1];
		im->FastGetSolutionStepValue(ADHESION_FORCE_Z) = 0.0;
	      }
	      else
	      {
		im->FastGetSolutionStepValue(ADHESION_FORCE_X) = 0.0;
		im->FastGetSolutionStepValue(ADHESION_FORCE_Y) = 0.0;
		im->FastGetSolutionStepValue(ADHESION_FORCE_Z) = 0.0;		
	      }
	      
	    }
	KRATOS_CATCH("")
    }    
    
    void Vector2D(const double x0, const double y0, const double x1, const double y1, array_1d<double,2>& r01)
    {
      r01[0] = x1 - x0;
      r01[1] = y1 - y0;
    }    
    
    void Vector3D(const double x0, const double y0, const double z0,
		  const double x1, const double y1, const double z1, array_1d<double,3>& r01)
    {
      r01[0] = x1 - x0;
      r01[1] = y1 - y0;
      r01[2] = z1 - z0;
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
    
   double SumVecs3D(const array_1d<double,3>& a, const array_1d<double,3>& b, array_1d<double,3>& c)
    {
      for (unsigned int i = 0; i < 3; i++)
      {
	c[i] = a[i] + b[i];
      }
    }    
    
    double DotProduct2D(const array_1d<double,2>& a, const array_1d<double,2>& b)
    {
      return (a[0]*b[0] + a[1]*b[1]);
    }    
    
   double fsign(const double& a, const double& b)
    {
      double c = a - b;
      if(c*c < 0.00000000000000001)
	return 0;
      else
	return c/sqrt(c*c);
    }        
    
    private:
    
  }; // Class CalculateAdhesionForce


}  // namespace Kratos.

#endif // KRATOS_CALCULATE_ADHESION_FORCE_INCLUDED defined 


