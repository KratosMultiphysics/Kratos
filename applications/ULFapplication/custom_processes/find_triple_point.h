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

#include <fstream>
#include <vector>
#include <ctime>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include <limits>

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
                        im->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
                        im->FastGetSolutionStepValue(CONTACT_ANGLE) = 0.0;
                        
                        if (im->FastGetSolutionStepValue(IS_WALL) != 0.0)
                        {
                            
                            
                            WeakPointerVector< Condition >& neighb_e = im->GetValue(NEIGHBOUR_CONDITIONS);
                            //Loop over faces -> find faces that two IS_FREE_SURFACE nodes
                            for (unsigned int i = 0; i < neighb_e.size(); i++)
                            {
                                int neighnum_st = 0;
                                int neighnum_fs = 0;
                                int visited = 0;	
                                int idx_j = 6;
                                int idx_k = 7;
                                for (unsigned int k = 0; k < neighb_e[i].GetGeometry().size(); k++)
                                {
                                    neighnum_st += neighb_e[i].GetGeometry()[k].FastGetSolutionStepValue(IS_WALL);
                                    neighnum_fs += neighb_e[i].GetGeometry()[k].FastGetSolutionStepValue(IS_FREE_SURFACE);
                                }
                                
                                if (neighnum_fs > 0)
                                {
                                    for (unsigned int j = 0; j < neighb_e[i].GetGeometry().size() ; j++)
                                    {
                                        if (neighb_e[i].GetGeometry()[j].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
                                            {
                                            if (visited < 1e-15)
                                            {
                                                idx_j = j;
                                            }
                                            else
                                            {	
                                                idx_k = j;
                                            }
                                            visited++;
                                            }
                                        }
                                    
                                    if (neighnum_fs > 0)
                                    {
                                        im->FastGetSolutionStepValue(TRIPLE_POINT) = 1.0;
                                    }
                                    else
                                    {
                                        im->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
                                        im->FastGetSolutionStepValue(CONTACT_ANGLE) = 0.0;
                                    }
                                }
                            }
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
                    
                    int is_free = 0;
                    int is_struct = 0;
                    int is_lagin = 0;
                    int is_boundary = 0;
      

                    for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; im++)
                    {
                        im->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
                        im->FastGetSolutionStepValue(CONTACT_ANGLE) = 0.0;
                        
                        if (im->FastGetSolutionStepValue(IS_WALL) != 0.0)
                        {
                            
                            
                            WeakPointerVector< Condition >& neighb_e = im->GetValue(NEIGHBOUR_CONDITIONS);
                            //Loop over faces -> find faces that two IS_FREE_SURFACE nodes
                            for (unsigned int i = 0; i < neighb_e.size(); i++)
                            {
                                int neighnum_st = 0;
                                int neighnum_fs = 0;
                                int visited = 0;	
                                int idx_j = 6;
                                int idx_k = 7;
                                for (unsigned int k = 0; k < neighb_e[i].GetGeometry().size(); k++)
                                {
                                    neighnum_st += neighb_e[i].GetGeometry()[k].FastGetSolutionStepValue(IS_WALL);
                                    neighnum_fs += neighb_e[i].GetGeometry()[k].FastGetSolutionStepValue(IS_FREE_SURFACE);
                                }
                                
                                if (neighnum_fs > 0)
                                {
                                    for (unsigned int j = 0; j < neighb_e[i].GetGeometry().size() ; j++)
                                    {
                                        if (neighb_e[i].GetGeometry()[j].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
                                            {
                                            if (visited < 1e-15)
                                            {
                                                idx_j = j;
                                            }
                                            else
                                            {	
                                                idx_k = j;
                                            }
                                            visited++;
                                            }
                                        }
                                    
                                    if (neighnum_fs > 0)
                                    {
                                        im->FastGetSolutionStepValue(TRIPLE_POINT) = 1.0;
                                    }
                                    else
                                    {
                                        im->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
                                        im->FastGetSolutionStepValue(CONTACT_ANGLE) = 0.0;
                                    }
                                }
                            }
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
