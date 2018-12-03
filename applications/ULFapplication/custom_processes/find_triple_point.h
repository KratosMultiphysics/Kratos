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
                    
                    int is_free = 0;
                    int is_struct = 0;
                    int is_lagin = 0;
                    int is_boundary = 0;
      

                    for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; im++)
                    {
                        im->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
                        im->FastGetSolutionStepValue(CONTACT_ANGLE) = 0.0;
                        
                        
                        
// // //                         double xi = im->X();
// // //                         double yi = im->Y();
// // //                         double zi = im->Z();
// // //                         
// // //                         
// // //                         double max_x, min_x, max_y, min_y, min_z;
// // //                         
// // //                         max_x = - 10000.0;
// // //                         min_x = 100000.0;
// // //                         max_y= - 100000.0;
// // //                         min_y = 100000.0;
// // //                         min_z = 100000.0;
// // //                         
// // //                         
// // //                         if(xi < min_x)
// // //                             min_x = xi;
// // // 
// // //                         if(xi > max_x)
// // //                             max_x = xi;
// // //                         
// // //                         
// // //                         
// // //                         if(yi < min_y)
// // //                             min_x = yi;
// // // 
// // //                         if(yi > max_y)
// // //                             max_y = yi;
// // //                         
// // //                         
// // //                         
// // //                         if(zi < min_z)
// // //                             min_z = zi;
// // //                         
// // //                         
// // //                         if (zi == min_z)
// // //                             {
// // //                                 im->FastGetSolutionStepValue(IS_STRUCTURE) = 1.0;
// // //                             }
                        
//                         int sum = 0;
                        
//                         x0=im->FastGetSolutionStepValue(DISPLACEMENT_X);
//                         
// 
// 
//                         if(ThisModelPart.Nodes().begin() == true)
//                         {
//                             min_x = x0;
//                             max_x = x0;
//                         }
//                         else
//                         {
//                             if(x0 < min_x)
//                                 min_x = x0;
// 
//                             if(x0 > max_x)
//                                 max_x = x0;
//                         }
//                         
//                     
//                         
//                         if (x0 <= max_x && x0 >= min_x)
//                         {
//                             for(ModelPart::NodesContainerType::iterator iy = ThisModelPart.NodesBegin() ; iy != ThisModelPart.NodesEnd() ; iy++)
//                             {
//                                 y0=iy->FastGetSolutionStepValue(DISPLACEMENT_Y);
//                                 
//                                 if(ThisModelPart.Nodes().begin() == true )
//                                 {
//                                     min_y = y0;
//                                     max_y = y0;
//                                 }
//                                 else
//                                 {
//                                     if(y0 < min_y)
//                                         min_y = y0;
// 
//                                     if(y0 > max_y)
//                                         max_y = y0;
//                                 }
//                         
//                                 if (y0 <= max_y && y0 >= min_y)
//                                 {
//                                     for(ModelPart::NodesContainerType::iterator iz = ThisModelPart.NodesBegin() ; iz != ThisModelPart.NodesEnd() ; iz++)
//                                     {
//                                         z0=iz->FastGetSolutionStepValue(DISPLACEMENT_Z);
//                                         
//                                         if(ThisModelPart.Nodes().begin() == true)
//                                         {
//                                             min_z = z0;
//                                         }
//                                         else
//                                         {
//                                             if(z0 < min_z)
//                                                 min_z = z0;
//                                         }
//                                         
//                                         if (z0 <= min_z)
//                                         {
//                                             iz->FastGetSolutionStepValue(IS_STRUCTURE) = 1.0;
//                                         }
// //                                         sum+=z0;
//                                     }
//                                 }
// //                                 sum+=y0;
//                             }
//                             
//                             
//                         }

//                         sum+=x0;
            
                        
                        if (im->FastGetSolutionStepValue(IS_STRUCTURE) != 0.0)
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
                                    neighnum_st += neighb_e[i].GetGeometry()[k].FastGetSolutionStepValue(IS_STRUCTURE);
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
