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

#if !defined(FIND_TRIPLE_POINT_CONDITION_INCLUDED )
#define  FIND_TRIPLE_POINT_CONDITION_INCLUDED



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
//#include "custom_conditions/Surface_Tension2D.h"



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
	find triple point for sessile droplets, based on the given geometries.
*/

  class FindTriplePoint
  : public Process
  {
  public:

 ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindTriplePoint);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FindTriplePoint()
    {	
    }

    /// Destructor.
    ~FindTriplePoint() override
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


    void FindTriplePoint2D(ModelPart& ThisModelPart)
    {
      KRATOS_TRY

      int is_free = 0;
      int is_struct = 0;
      int is_lagin = 0;
      
      for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; im++)
      {
	  if (im->FastGetSolutionStepValue(IS_STRUCTURE) != 0.0)
	  {
	      is_free = 0;
	      is_struct = 0;
	      is_lagin = 0;
	      WeakPointerVector< Node<3> >& neighb = im->GetValue(NEIGHBOUR_NODES);
	      for (unsigned int i = 0; i < neighb.size(); i++)
	      {
		  if (neighb[i].FastGetSolutionStepValue(IS_BOUNDARY) != 0.0)
		  {
		    if (neighb[i].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
			is_free++;
		    if (neighb[i].FastGetSolutionStepValue(IS_STRUCTURE) != 0.0)
			is_struct++;
		    if (neighb[i].FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) != 0.0)
			is_lagin++;
		  }
	      }
	      if (is_free == 1 && is_struct == 1)
		  im->FastGetSolutionStepValue(TRIPLE_POINT) = 1.0;
	  }
      }
      
      
      KRATOS_CATCH("")
    }
    
    void FindTriplePoint3D(ModelPart& ThisModelPart)
    {
      KRATOS_TRY
      
      unsigned int num_tp = 0;
      unsigned int num_fs = 0;
      array_1d<double,3> An = ZeroVector(3);
      
      for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; im++)
      {
	  num_fs = num_tp = 0;
//  	  if ((im->FastGetSolutionStepValue(TRIPLE_POINT) != 0.0) && (im->FastGetSolutionStepValue(CONTACT_ANGLE) == 0.0) && (im->FastGetSolutionStepValue(VELOCITY_X) != 0.0))
	  if ((im->FastGetSolutionStepValue(TRIPLE_POINT) != 0.0) && (im->FastGetSolutionStepValue(CONTACT_ANGLE) < 1e-15) && (im->FastGetSolutionStepValue(VELOCITY_X) != 0.0))
	      im->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
	  
	  An = im->FastGetSolutionStepValue(NORMAL);
	  NormalizeVec3D(An);
	  if (im->FastGetSolutionStepValue(TRIPLE_POINT) != 0.0 && An[2] < -0.99)
	  {
	      im->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
	      im->FastGetSolutionStepValue(CONTACT_ANGLE) = 0.0;
	      im->FastGetSolutionStepValue(FORCE) = ZeroVector(3);
	  }
	  
	  if (im->FastGetSolutionStepValue(IS_STRUCTURE) != 0.0)
	  {
	    if( An[0] > 0.01 || An[0] < -0.01 || An[1] > 0.01 || An[1] < -0.01)
	    {

	      im->FastGetSolutionStepValue(TRIPLE_POINT) = 1.0;
	    }
	    else
		im->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
	  }
      }

      KRATOS_CATCH("")
    }    

    double Norm3D(const array_1d<double,3>& a)
    {
      return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    }

    double DotProduct3D(const array_1d<double,3>& a, const array_1d<double,3>& b)
    {
      return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    }   
    
    void NormalizeVec3D(array_1d<double,3>& input)
    {
      double norm = Norm3D(input);
      if (norm != 0.0)
      {
	input[0] /= norm;
	input[1] /= norm;
	input[2] /= norm;
      }
    } 
    
    private:
    
  }; // Class FindTriplePoint


}  // namespace Kratos.

#endif // KRATOS_CREATE_INTERFACE_CONDITIONS_PROCESS_INCLUDED  defined 
