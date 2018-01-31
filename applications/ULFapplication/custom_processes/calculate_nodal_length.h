//
//   Project Name:        Kratos
//   Last Modified by:    $Author: ajarauta $
//   Date:                $Date: 2007-11-06 12:34:26 $
//   Revision:            $Revision: 1.4 $
//
//  this process save structural elements in a separate list

#if !defined(CALCULATE_NODAL_LENGTH_INCLUDED )
#define  CALCULATE_NODAL_LENGTH_INCLUDED



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
  class CalculateNodalLength
  : public Process
  {
  public:
    void CalculateNodalLength2D(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	
	double x0,y0,x1,y1,x2,y2,norm10,norm20;
	int neighnum = 0;
	array_1d<double,2> r10 = ZeroVector(2);
	array_1d<double,2> r20 = ZeroVector(2);
	
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; im++)
	    {
	      //Find the neighbours of TRIPLE_POINT at the boundary
	      if ((im->FastGetSolutionStepValue(TRIPLE_POINT))*1000 != 0.0)
	      {
		WeakPointerVector< Node<3> >& neighb = im->GetValue(NEIGHBOUR_NODES);
		x0 = im->X();
		y0 = im->Y();
		x1 = 0.0;
		y1 = 0.0;
		x2 = 0.0;
		y2 = 0.0;
		neighnum = 0;
		for (unsigned int i = 0; i < neighb.size(); i++)
		{
		  if (neighb[i].FastGetSolutionStepValue(IS_BOUNDARY) != 0.0)
		  {
		    if (neighnum == 0)
		    {
		      x1 = neighb[i].X();
		      y1 = neighb[i].Y();
		    }
		    if (neighnum == 1)
		    {
		      x2 = neighb[i].X();
		      y2 = neighb[i].Y();
		    }
		    neighnum++;
		  }
		}
		
		//Obtain the vectors pointing from node 0 to node 1 (r10) and from 0 to 2 (r20)
		Vector2D(x0,y0,x1,y1,r10);
		Vector2D(x0,y0,x2,y2,r20);
		norm10 = Norm2D(r10);
		norm20 = Norm2D(r20);
		im->FastGetSolutionStepValue(NODAL_LENGTH) = 0.5*(norm10+norm20);
	      }
	      
	      //Repeat for IS_FREE_SURFACE nodes
	      //Find the neighbours of IS_FREE_SURFACE at the boundary
	      if (im->FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
	      {
		WeakPointerVector< Node<3> >& neighb = im->GetValue(NEIGHBOUR_NODES);
		x0 = im->X();
		y0 = im->Y();
		x1 = 0.0;
		y1 = 0.0;
		x2 = 0.0;
		y2 = 0.0;
		neighnum = 0;
		for (unsigned int i = 0; i < neighb.size(); i++)
		{
		  if ((neighb[i].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0) || ((im->FastGetSolutionStepValue(TRIPLE_POINT))*1000 != 0.0 ))
		  {
		    if (neighnum == 0)
		    {
		      x1 = neighb[i].X();
		      y1 = neighb[i].Y();
		    }
		    if (neighnum == 1)
		    {
		      x2 = neighb[i].X();
		      y2 = neighb[i].Y();
		    }
		    neighnum++;
		  }
		}
		
		//Obtain the vectors pointing from node 0 to node 1 (r10) and from 0 to 2 (r20)
		Vector2D(x0,y0,x1,y1,r10);
		Vector2D(x0,y0,x2,y2,r20);
		norm10 = Norm2D(r10);
		norm20 = Norm2D(r20);
		im->FastGetSolutionStepValue(NODAL_LENGTH) = 0.5*(norm10+norm20);
	      }	   

		
	    }
	    
	    
	KRATOS_CATCH("")
    }
    
    void CalculateNodalLength3D(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	

	
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ;
	    im != ThisModelPart.NodesEnd() ; im++)
	    {
	      //Find the neighbours of TRIPLE_POINT at the boundary
	      if ((im->FastGetSolutionStepValue(TRIPLE_POINT))*1000 != 0.0)
	      {
		WeakPointerVector< Node<3> >& neighb = im->GetValue(NEIGHBOUR_NODES);
		double x0 = im->X();
		double y0 = im->Y();
		double z0 = im->Z();
		double x1 = 0.0;
		double y1 = 0.0;
		double z1 = 0.0;
		double x2 = 0.0;
		double y2 = 0.0;
		double z2 = 0.0;
		int neighnum = 0;
		for (unsigned int i = 0; i < neighb.size(); i++)
		{
		  if ((neighb[i].FastGetSolutionStepValue(TRIPLE_POINT))*1000.0 != 0.0)
		  {
		    if (neighb[i].X() != x0 || neighb[i].Y() != y0 || neighb[i].Z() != z0)
		    {
		      if (neighnum == 0)
		      {
			x1 = neighb[i].X();
			y1 = neighb[i].Y();
			z1 = neighb[i].Z();
		      }
		      else
		      {
			x2 = neighb[i].X();
			y2 = neighb[i].Y();
			z2 = neighb[i].Z();		      
		      }
		      neighnum++;
		    }
		  }
		}
		
		//Obtain the vectors pointing from node 0 to node 1 (r10) and from 0 to 2 (r20)
		array_1d<double,3> r10 = ZeroVector(3);
		array_1d<double,3> r20 = ZeroVector(3);
		Vector3D(x0,y0,z0,x1,y1,z1,r10);
		Vector3D(x0,y0,z0,x2,y2,z2,r20);
		double norm10 = Norm3D(r10);
		double norm20 = Norm3D(r20);
		
		im->FastGetSolutionStepValue(NODAL_LENGTH) = 0.5*(norm10+norm20);
// 		KRATOS_WATCH(im->FastGetSolutionStepValue(NODAL_LENGTH))
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
    
    private:
    
  }; // Class CalculateNodalLength


}  // namespace Kratos.

#endif // KRATOS_CREATE_INTERFACE_CONDITIONS_PROCESS_INCLUDED  defined 


