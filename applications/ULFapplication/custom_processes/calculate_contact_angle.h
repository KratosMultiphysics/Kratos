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


#if !defined(CALCULATE_CONTACT_ANGLE_INCLUDED )
#define  CALCULATE_CONTACT_ANGLE_INCLUDED



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
	calculate contact angle for every tiple point, based on the given domain configuration


*/

  class CalculateContactAngle
  : public Process
  {
  public:

///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateContactAngle);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CalculateContactAngle()
    {	
    }

    /// Destructor.
    ~CalculateContactAngle() override
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

   void CalculateContactAngle2D(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	
	double theta = 0.0;
	double pi = 3.14159265359;
	
// 	ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();
// 	double theta_s = CurrentProcessInfo[CONTACT_ANGLE_STATIC];
// 	double delta_t = CurrentProcessInfo[DELTA_TIME];
// 	double theta_adv = theta_s + 10.0;
// 	double theta_rec = theta_s - 10.0;
	
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ;
	    im != ThisModelPart.NodesEnd() ; im++)
	    {
	      //Find the neighbours of TRIPLE_POINT at the boundary
	      if ((im->FastGetSolutionStepValue(TRIPLE_POINT))*1000 != 0.0)
	      {
		WeakPointerVector< Node<3> >& neighb = im->GetValue(NEIGHBOUR_NODES);
		double x0 = im->X();
		double y0 = im->Y();
		double x1 = 0.0;
		double y1 = 0.0;
		double x2 = 0.0;
		double y2 = 0.0;
// 		double vx1 = 0.0;
		int neighnum = 0;
		for (unsigned int i = 0; i < neighb.size(); i++)
		{
		  if (neighb[i].FastGetSolutionStepValue(IS_BOUNDARY) != 0.0)
		  {
//  		    if (neighnum == 0)
		    if (neighnum < 1e-15)
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
// 		    if (neighb[i].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
// 		    {
// 		      vx1 = neighb[i].FastGetSolutionStepValue(VELOCITY_X);
// 		    }
		  }
		}
		
		//Obtain the vectors pointing from node 0 to node 1 (r10) and from 0 to 2 (r20)
		array_1d<double,2> r10 = ZeroVector(2);
		array_1d<double,2> r20 = ZeroVector(2);
		r10[0] = x1 - x0;
		r10[1] = y1 - y0;
		r20[0] = x2 - x0;
		r20[1] = y2 - y0;
		double norm10 = sqrt(r10[0]*r10[0] + r10[1]*r10[1]);
		double norm20 = sqrt(r20[0]*r20[0] + r20[1]*r20[1]);
		
		double costheta = (r10[0]*r20[0] + r10[1]*r20[1])/(norm10*norm20);
		theta = acos(costheta)*180.0/pi;
		im->FastGetSolutionStepValue(CONTACT_ANGLE) = theta;
		
// 		if(theta > theta_adv || theta < theta_rec)
// 		{
// 		  im->FastGetSolutionStepValue(DISPLACEMENT_X) = vx1*delta_t;
// 		}
	      }
	      else
		im->FastGetSolutionStepValue(CONTACT_ANGLE) = 0.0;
	    }
	    /*
	    //Now we check that there are only TWO nodes with CONTACT_ANGLE != 0.0
	    int num_trip = 0;
	    double min_x = 0.0;
	    double max_x = 0.0;
	    for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; im++)
	    {
		if (im->FastGetSolutionStepValue(CONTACT_ANGLE) != 0.0)
		{
		    ++num_trip;
		    if (num_trip < 2)
		    {
			min_x = im->X();
			max_x = min_x;
		    }
		}
	    }
	    if (num_trip > 2)
	    {
		for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; im++)
		{
		    if (im->FastGetSolutionStepValue(CONTACT_ANGLE) != 0.0)
		    {
			if (im->X() < min_x) //left triple point
			{
			    min_x = im->X();
			}
			if (im->X() > max_x) //right triple point
			{
			    max_x = im->X();
			}
		    }
		}
		//Now we have the coordinates of the true left and right triple points -> remove others
		for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; im++)
		{
		    if (im->FastGetSolutionStepValue(CONTACT_ANGLE) != 0.0)
		    {
			if (im->X() >  min_x && im->X() < max_x)
			      im->FastGetSolutionStepValue(CONTACT_ANGLE) = 0.0;
		    }
		}
	    }
	    */
	
	KRATOS_CATCH("")
    }
    
    void CalculateContactAngle3D(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	
	double theta = 0.0;
	double pi = 3.14159265359;
	double theta_eq = ThisModelPart.GetProcessInfo()[CONTACT_ANGLE_STATIC];
	double theta_rad = theta_eq*pi/180.0;

	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; im++)
	{
	      if (im->FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		///////////////////////////////////////////////////////////////////////////////////////////////
		// For node i, you find all its neighbor conditions. This is a set of triangles (faces)
		// around node i. Only faces with two IS_FREE_SURFACE nodes will be used to obtain contact angle
		// For every face, cross product with the two adjacent nodes j and k gives the vector normal
		// to the face -> ALWAYS do it such that this vector points outwards the domain:
		// if (rij X rik) · NORMAL_GEOMETRIC < 0 -> vector is pointing inwards -> do rik X rij
		//////////////////////////////////////////////////////////////////////////////////////////
		
		double xi = im->X();
		double yi = im->Y();
		double zi = im->Z();
		double xj, yj, zj, xk, yk, zk;
		xj = yj = zj = xk = yk = zk = 0.0;		
		int idx_j = 6;
		int idx_k = 7;
		
		array_1d<double,3> rij = ZeroVector(3);
		array_1d<double,3> rik = ZeroVector(3);
		array_1d<double,3> cross_prod_ijk = ZeroVector(3);
		array_1d<double,3> normal_geom = ZeroVector(3);
		array_1d<double,3> normal_tp = ZeroVector(3);
		normal_geom = im->FastGetSolutionStepValue(NORMAL_GEOMETRIC);
		double dot_prod = 0.0;
		
		im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE) = ZeroVector(3);
		array_1d<double,2> aux = ZeroVector(2);
		array_1d<double,3> temp = ZeroVector(3);
		
		int neighnum_tp = 0;
		int neighnum_caf = 0; //caf == contact angle face
		int visited = 0;
		double num_faces = 0.0;
		WeakPointerVector< Condition >& neighb_faces = im->GetValue(NEIGHBOUR_CONDITIONS);
		//Loop over faces -> find faces that two IS_FREE_SURFACE nodes
		for (unsigned int i = 0; i < neighb_faces.size(); i++)
		{
		  neighnum_tp = 0;
		  neighnum_caf = 0;
		  visited = 0;
		  for (unsigned int k = 0; k < neighb_faces[i].GetGeometry().size(); k++)
		  {
		    neighnum_tp += neighb_faces[i].GetGeometry()[k].FastGetSolutionStepValue(TRIPLE_POINT);
		    neighnum_caf += neighb_faces[i].GetGeometry()[k].FastGetSolutionStepValue(IS_FREE_SURFACE);
		  }
		  
		  num_faces++;
		  
		  if (neighnum_caf == 2)
		  {
		    
		    for (unsigned int j = 0; j < neighb_faces[i].GetGeometry().size() ; j++)
		    {
			if (neighb_faces[i].GetGeometry()[j].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
			{
//  			  if (visited == 0)
			  if (visited < 1e-15)
			  {
			    xj = neighb_faces[i].GetGeometry()[j].X();
			    yj = neighb_faces[i].GetGeometry()[j].Y();
			    zj = neighb_faces[i].GetGeometry()[j].Z();
			    idx_j = j;
			  }
			  else
			  {
			    xk = neighb_faces[i].GetGeometry()[j].X();
			    yk = neighb_faces[i].GetGeometry()[j].Y();
			    zk = neighb_faces[i].GetGeometry()[j].Z();	
			    idx_k = j;
			  }
			  visited++;
			}
		    }
		  }
		  
		  if (neighnum_caf == 1)
		  {
		    for (unsigned int j = 0; j < neighb_faces[i].GetGeometry().size() ; j++)
		    {
		      if (neighb_faces[i].GetGeometry()[j].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
		      {
			  xj = neighb_faces[i].GetGeometry()[j].X();
			  yj = neighb_faces[i].GetGeometry()[j].Y();
			  zj = neighb_faces[i].GetGeometry()[j].Z();
		      }
		      if (neighb_faces[i].GetGeometry()[j].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0 && neighb_faces[i].GetGeometry()[j].Id() != im->Id())
		      {
			  xk = neighb_faces[i].GetGeometry()[j].X();
			  yk = neighb_faces[i].GetGeometry()[j].Y();
			  zk = neighb_faces[i].GetGeometry()[j].Z();	
		      }
		    }
		  }
		  
		  //Now we have neighbor nodes of considered face -> obtain normal vector
		  Vector3D(xi,yi,zi,xj,yj,zj,rij);
		  Vector3D(xi,yi,zi,xk,yk,zk,rik);
		  CrossProduct3D(rij, rik, cross_prod_ijk);
		  dot_prod = DotProduct3D(cross_prod_ijk,normal_geom);
		  if (dot_prod < 0.0)
		      CrossProduct3D(rik, rij, cross_prod_ijk);
		  //We add this vector to the normal_tp
		  normal_tp[0] += cross_prod_ijk[0];
		  normal_tp[1] += cross_prod_ijk[1];
		  normal_tp[2] += cross_prod_ijk[2];
		  //Finally, contact angle is obtained -> it is simply acos(dotproduct(cross_prod_ijk,(0,0,1))/norm(cross_prod_ijk))
		  NormalizeVec3D(cross_prod_ijk);
		  theta = acos(cross_prod_ijk[2])*180.0/pi;
		  im->FastGetSolutionStepValue(CONTACT_ANGLE) += theta;
		  
		  //Normal to the contact line
		  if (neighnum_caf == 1)
		      im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE) += rij;
		  if (neighnum_caf == 2)
		      im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE) += 0.5*(rij+rik);		  
		  
		}
		NormalizeVec3D(normal_tp);
		im->FastGetSolutionStepValue(NORMAL_TRIPLE_POINT_X) = normal_tp[0];
		im->FastGetSolutionStepValue(NORMAL_TRIPLE_POINT_Y) = normal_tp[1];
		im->FastGetSolutionStepValue(NORMAL_TRIPLE_POINT_Z) = normal_tp[2];
		
		aux[0] = normal_tp[0];
		aux[1] = normal_tp[1];
		NormalizeVec2D(aux);
		NormalizeVec3D(im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE));
		temp = im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE);
		
		im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE_X) = cos(asin(temp[2]))*aux[0];
		im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE_Y) = cos(asin(temp[2]))*aux[1];
		NormalizeVec3D(im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE));
		
		if (num_faces != 0.0)
		    im->FastGetSolutionStepValue(CONTACT_ANGLE) /= num_faces;
		
		if((im->FastGetSolutionStepValue(CONTACT_ANGLE)) < 90.0)
		{
		    im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE_X) *= -1.0;
		    im->FastGetSolutionStepValue(NORMAL_CONTACT_LINE_Y) *= -1.0;
		}
	      }
	}
	
	KRATOS_CATCH("")
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// AUXILIAR FUNCTIONS ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
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
    
    double DotProduct2D(const array_1d<double,2>& a, const array_1d<double,2>& b)
    {
      return (a[0]*b[0] + a[1]*b[1]);
    }
    
    double DotProduct3D(const array_1d<double,3>& a, const array_1d<double,3>& b)
    {
      return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    }
    
    void CrossProduct3D(const array_1d<double,3>& a, const array_1d<double,3>& b, array_1d<double,3>& c)
    {
      c[0] = a[1]*b[2] - a[2]*b[1];
      c[1] = a[2]*b[0] - a[0]*b[2];
      c[2] = a[0]*b[1] - a[1]*b[0];
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
//		CalculateContactAngle& operator=(CalculateContactAngle const& rOther);

    /// Copy constructor.
//		CalculateContactAngle(CalculateContactAngle const& rOther);


    ///@}

}; // Class CalculateContactAngle

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateContactAngle& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateContactAngle& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



}  // namespace Kratos.

#endif // KRATOS_CREATE_INTERFACE_CONDITIONS_PROCESS_INCLUDED  defined 


