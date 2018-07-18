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



#if !defined(CALCULATE_NODAL_LENGTH_INCLUDED )
#define  CALCULATE_NODAL_LENGTH_INCLUDED



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
	calculate nodal length_h, used for droplet dynamics 


*/


  class CalculateNodalLength
  : public Process
  {
  public:

///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateNodalLength);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CalculateNodalLength()
    {	
    }

    /// Destructor.
    ~CalculateNodalLength() override
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

    void CalculateNodalLength2D(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	
	double x0,y0,x1,y1,x2,y2,norm10,norm20;
	int neighnum = 0;
	array_1d<double,2> r10 = ZeroVector(2);
	array_1d<double,2> r20 = ZeroVector(2);
	
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; ++im)
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
	    im != ThisModelPart.NodesEnd() ; ++im)
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
        return "CalculateNodalLength";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CalculateNodalLength";
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
//		CalculateNodalLength& operator=(CalculateNodalLength const& rOther);

    /// Copy constructor.
//		CalculateNodalLength(CalculateNodalLength const& rOther);


    ///@}

}; // Class CalculateNodalLength

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateNodalLength& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateNodalLength& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



}  // namespace Kratos.


#endif // KRATOS_CREATE_INTERFACE_CONDITIONS_PROCESS_INCLUDED  defined 


