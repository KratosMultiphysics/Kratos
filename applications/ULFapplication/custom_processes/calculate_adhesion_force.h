//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Alex Jarauta
//  Co-author  :     Elaf Mahrous


#if !defined(CALCULATE_ADHESION_FORCE_INCLUDED )
#define  CALCULATE_ADHESION_FORCE_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include <pybind11/pybind11.h>
#include "includes/define.h"
#include "includes/define_python.h"

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
	calculating adhesion force for droplets, based on the given domain configuration


*/


  class CalculateAdhesionForce
  : public Process
  {
  public:
      
///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateAdhesionForce);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CalculateAdhesionForce()
    {
    }

    /// Destructor.
    ~CalculateAdhesionForce() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    //	void operator()()
    //	{
    //		MergeParts();
    //	}


    ///@}
    ///@name Operations
    ///@{
    
    
    void CalculateAdhesionForce3D(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	
	double gamma = ThisModelPart.GetProcessInfo()[SURFACE_TENSION_COEF]; //Surface tension coefficient [N m-1]
	double pi = 3.14159265359;
	///////double theta_eq = ThisModelPart.GetProcessInfo()[CONTACT_ANGLE_STATIC];
	double theta, sign_force, theta_rad, cos_t;
	array_1d<double,2> nu0 = ZeroVector(2);
	array_1d<double,2> vel = ZeroVector(2);
	array_1d<double,2> n_tp = ZeroVector(2);
	
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ;
	    im != ThisModelPart.NodesEnd() ; ++im)
	    {
	      //Find the neighbours of TRIPLE_POINT at the boundary
	      if ((im->FastGetSolutionStepValue(TRIPLE_POINT))*1000 != 0.0)
	      {
		nu0[0] = im->FastGetSolutionStepValue( NORMAL_TRIPLE_POINT_X);
		nu0[1] = im->FastGetSolutionStepValue( NORMAL_TRIPLE_POINT_Y);
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
		    n_tp[0] = im->FastGetSolutionStepValue( NORMAL_TRIPLE_POINT_X);
		    n_tp[1] = im->FastGetSolutionStepValue( NORMAL_TRIPLE_POINT_Y);
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
    std::string Info() const override
    {
        return "CalculateAdhesionForce";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CalculateAdhesionForce";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
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

    /// Assignment operator.
//		CalculateAdhesionForce& operator=(CalculateAdhesionForce const& rOther);

    /// Copy constructor.
//		CalculateAdhesionForce(CalculateAdhesionForce& rOther);


    ///@}

}; // Class CalculateAdhesionForce

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateAdhesionForce& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateAdhesionForce& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



}  // namespace Kratos.



#endif // KRATOS_CALCULATE_ADHESION_FORCE_INCLUDED defined 


