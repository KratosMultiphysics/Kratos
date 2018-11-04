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



#if !defined(CALCULATE_NORMAL_EQ_INCLUDED )
#define  CALCULATE_NORMAL_EQ_INCLUDED



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
	calculate normal equilibrium based on given geometries, contact angle, tripple point. It is used in droplet dynamis.


*/

  class CalculateNormalEq
  : public Process
  {
  public:

///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateNormalEq);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CalculateNormalEq()
    {	
    }

    /// Destructor.
    ~CalculateNormalEq() override
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
    

    void CalculateNormalEq3D(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	
	double pi = 3.14159265359;
	double theta_eq = ThisModelPart.GetProcessInfo()[CONTACT_ANGLE_STATIC];
	double theta_rad = theta_eq*pi/180.0;
	///////double theta = 0.0;
	
	array_1d<double,3> normal_eq = ZeroVector(3);
	array_1d<double,3> normal_xy = ZeroVector(3);
	normal_xy[2] = 1.0;
	array_1d<double,2> temp = ZeroVector(2);
	
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; ++im)
	{
	    //Find the neighbours of TRIPLE_POINT at the boundary
	    if (im->FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	    {
		temp[0] = im->FastGetSolutionStepValue(NORMAL_TRIPLE_POINT_X);
		temp[1] = im->FastGetSolutionStepValue(NORMAL_TRIPLE_POINT_Y);
		temp = NormalizeVec2D(temp);
		normal_eq[0] = sin(theta_rad)*temp[0];
		normal_eq[1] = sin(theta_rad)*temp[1];
		normal_eq[2] = cos(theta_rad);
		normal_eq = NormalizeVec3D(normal_eq);	      
		im->FastGetSolutionStepValue(NORMAL_EQUILIBRIUM) = normal_eq;
		
		im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE_EQUILIBRIUM_X) = -cos(theta_rad)*temp[0];
		im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE_EQUILIBRIUM_Y) = -cos(theta_rad)*temp[1];		
		im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE_EQUILIBRIUM_Z) = sin(theta_rad);
		im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE_EQUILIBRIUM) = NormalizeVec3D(im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE_EQUILIBRIUM));

	    }
	    else
		im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE_EQUILIBRIUM) = ZeroVector(3);
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
        return "AssignSurfaceTensionConditions";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignSurfaceTensionConditions";
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
//		CalculateNormalEq& operator=(CalculateNormalEq const& rOther);

    /// Copy constructor.
//		CalculateNormalEq(CalculateNormalEq const& rOther);


    ///@}

}; // Class CalculateNormalEq

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateNormalEq& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateNormalEq& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



}  // namespace Kratos.


#endif // KRATOS_CALCULATE_NORMAL_EQ_INCLUDED  defined 


