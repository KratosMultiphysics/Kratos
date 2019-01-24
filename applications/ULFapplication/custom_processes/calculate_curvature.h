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


#if !defined(CALCULATE_CURVATURE_INCLUDED )
#define  CALCULATE_CURVATURE_INCLUDED



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
#include "utilities/openmp_utils.h"
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
	calculate curvature for 2D and 3D.
	
	curvature for 2D is calcualted based on equation (14) in the follwoing reference:
	Jarauta A, Ryzhakov P, Secanell M, Waghmare PR, Pons-Prats J. Numerical study of droplet dynamics in a polymer electrolyte fuel cell gas channel using an embedded Eulerian-Lagrangian approach. Journal of Power Sources. 2016 Aug 15;323:201-12.
	
	curvature for 3D is based on Mayer approach in the follwoing reference:
	M. Meyer, M. Desbrun, P. Schröder, and A. H. Barr. Discrete differential geometry operators for triangulated 2-manifolds. Visualization and Math. III, pages 35–57, 2003.
*/

  class CalculateCurvature
  : public Process
  {
  public:

///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateCurvature);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CalculateCurvature()
    {	
    }

    /// Destructor.
    virtual ~CalculateCurvature()
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

    void CalculateCurvature2D(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	
	double theta_eq = ThisModelPart.GetProcessInfo()[CONTACT_ANGLE_STATIC];
	double pi = 3.14159265359;
	double theta_w = theta_eq*pi/180.0;
	double x0,y0,x1,y1,x2,y2;
	int neighnum = 0;
	
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; ++im)
	    {
	      if (im->FastGetSolutionStepValue(IS_INTERFACE) != 0.0 && im->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) < 1e-15)
	      {
		WeakPointerVector< Node<3> >& neighb = im->GetValue(NEIGHBOUR_NODES);
		x0 = im->X();
		y0 = im->Y();
		neighnum = 0;

		for (unsigned int i = 0; i < neighb.size(); i++)
		{
		    if (neighb[i].FastGetSolutionStepValue(IS_INTERFACE) != 0.0 || neighb[i].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
		    {
		      if (neighb[i].X() != x0 || neighb[i].Y() != y0)
		      {
 			if (neighnum == 0)
                        //if (neighnum < 1e-15)
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
		}
		double curv = CalculateCurv(x0,y0,x1,y1,x2,y2);
		  
		// Sign of the curvature
		double sign_curv = 1.0;
		double xc = 0.5*(x1 + x2);
		double yc = 0.5*(y1 + y2);
		int curv_pos = 0;
		int struct_nodes = 0;
		WeakPointerVector< Element >& neighb_els = im->GetValue(NEIGHBOUR_ELEMENTS);
		for (unsigned int i = 0; i < neighb_els.size(); i++)
		{
		  int intf_nodes = 0;
		  for (unsigned int j = 0; j < 3; j++)
		  {
		    if (neighb_els[i].GetGeometry()[j].FastGetSolutionStepValue(IS_INTERFACE) != 0.0 || neighb_els[i].GetGeometry()[j].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
		      intf_nodes++;
		    if (neighb_els[i].GetGeometry()[j].FastGetSolutionStepValue(IS_STRUCTURE) != 0.0)
		      struct_nodes++;
		  }
		  if (intf_nodes > 0)
		    neighb_els[i].GetValue(IS_WATER_ELEMENT) = 1.0;
		  //Calculate N at (xc,yc)
		  double x0_loc = neighb_els[i].GetGeometry()[0].X();
		  double y0_loc = neighb_els[i].GetGeometry()[0].Y();
		  double x1_loc = neighb_els[i].GetGeometry()[1].X();
		  double y1_loc = neighb_els[i].GetGeometry()[1].Y();
		  double x2_loc = neighb_els[i].GetGeometry()[2].X();
		  double y2_loc = neighb_els[i].GetGeometry()[2].Y();
		  
		  double x10 = x1_loc - x0_loc;
		  double y10 = y1_loc - y0_loc;
		  
		  double x20 = x2_loc - x0_loc;
		  double y20 = y2_loc - y0_loc;
		  
		  double detJ = x10 * y20-y10 * x20;
		  double totArea=0.5*detJ;
		  
		  array_1d<double,3> N = ZeroVector(3);
		  N[0] = CalculateVol(x1_loc,y1_loc,x2_loc,y2_loc,xc,yc)  / totArea;
		  N[1] = CalculateVol(x2_loc,y2_loc,x0_loc,y0_loc,xc,yc)  / totArea;
		  N[2] = CalculateVol(x0_loc,y0_loc,x1_loc,y1_loc,xc,yc)  / totArea;
		  
		  double eps_sign = -0.00000001;
		  if (N[0] >= eps_sign && N[1] >= eps_sign && N[2] >= eps_sign)
		  {
		    if (neighb_els[i].GetValue(IS_WATER_ELEMENT) == 1.0)
		    {
		      curv_pos++;
		    }
		  }
		}
 		if (curv_pos == 0)
                //if (curv_pos < 1e-15)
		{
		  sign_curv = -1.0;
		}
		
		im->FastGetSolutionStepValue(MEAN_CURVATURE_2D) = sign_curv*curv;
		
	      }
	      
	      //20140923 Computing curvature at the triple point!
	      // Curvature at triple point must be equal to the rate of change of normal vector,
	      // using neighbor IS_FREE_SURFACE and n = nw*cos_theta + tw*sin_theta
	      if (im->FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		x0 = im->X();
		y0 = im->Y();
		double sin_t = sin(theta_w+0.5*pi);
		double cos_t = cos(theta_w+0.5*pi);
		array_1d<double,2> n0 = ZeroVector(2);
		array_1d<double,2> n1 = ZeroVector(2);
		array_1d<double,2> struct_dist = ZeroVector(2);
		double dist_struct = 0.0;
		WeakPointerVector< Node<3> >& neighb = im->GetValue(NEIGHBOUR_NODES);
		for (unsigned int i = 0; i < neighb.size(); i++)
		{
		  if (neighb[i].FastGetSolutionStepValue(IS_INTERFACE) != 0.0)
		  {
		      x1 = neighb[i].X();
		      y1 = neighb[i].Y();
		      n1[0] = neighb[i].FastGetSolutionStepValue(NORMAL_X);
		      n1[1] = neighb[i].FastGetSolutionStepValue(NORMAL_Y);
// 		      n0[1] = cos_t;
		      n0[1] = sin_t;
		      if (n1[0] > 0.0)
			  n0[0] = -cos_t;
// 			  n0[0] = sin_t;			    
		      else
			  n0[0] = cos_t;
// 			  n0[0] = -sin_t;
		  }
		  if (neighb[i].FastGetSolutionStepValue(IS_STRUCTURE) != 0.0)
		  {
		      x2 = neighb[i].X();
		      y2 = neighb[i].Y();
		      Vector2D(x0,y0,x2,y2,struct_dist);
		      dist_struct = Norm2D(struct_dist);
		      if (dist_struct < 1.0e-5)
			neighb[i].Set(TO_ERASE,true);
		  }
		}
		array_1d<double,2> ds_vec = ZeroVector(2);
		Vector2D(x0,y0,x1,y1,ds_vec);
		double ds = Norm2D(ds_vec);
		NormalizeVec2D(n1);
		array_1d<double,2> dT = ZeroVector(2);
		dT[0] = n1[0] - n0[0];
		dT[1] = n1[1] - n0[1];
		double dT_norm = Norm2D(dT);
		double curv_tp = dT_norm/ds;
		// There is no need to compute curvature's sign: at triple point it's always positive!
		im->FastGetSolutionStepValue(MEAN_CURVATURE_2D) = curv_tp;			
	      }
	      
 	      if(im->FastGetSolutionStepValue(IS_STRUCTURE) != 0.0 && im->FastGetSolutionStepValue(TRIPLE_POINT)*1.0E8 == 0.0)
              //if(im->FastGetSolutionStepValue(IS_STRUCTURE) != 0.0 && im->FastGetSolutionStepValue(TRIPLE_POINT)*1.0E8 < 1e-15)
		  im->FastGetSolutionStepValue(MEAN_CURVATURE_2D) = 0.0;
	    }
	    
	KRATOS_CATCH("")
    }      
    
    void CalculateCurvature3D(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	
	double pi = 3.14159265;
	double hpi = 0.5*pi;
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ;
	    im != ThisModelPart.NodesEnd() ; ++im)
	    {
	      if (im->FastGetSolutionStepValue(FLAG_VARIABLE) > 0.999)// || im->FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		//////////////////////////////////////////////////////////////////////////////////////////
		// For node i, you find all its neighbor conditions. This is a set of triangles (faces)
		// around node i. All the nodes surrounding node i are its 1-ring neighborhood.
		// For every face, the two adjacent nodes are nodes j and k -> A_voronoi is obtained in 
		// every triangle and added to the total area.
		//////////////////////////////////////////////////////////////////////////////////////////
		
		double xi = im->X();
		double yi = im->Y();
		double zi = im->Z();
		double xj = 0.0;
		double yj = 0.0;
		double zj = 0.0;
		double xk = 0.0;
		double yk = 0.0;
		double zk = 0.0;	
		
		array_1d<double,3> rij = ZeroVector(3);
		array_1d<double,3> rji = ZeroVector(3);
		array_1d<double,3> rik = ZeroVector(3);
		array_1d<double,3> rki = ZeroVector(3);
		array_1d<double,3> rjk = ZeroVector(3);
		double norm_ij = 0.0;
		double norm_ik = 0.0;
		double norm_jk = 0.0;
		array_1d<double,3> cross_prod_ijk = ZeroVector(3);
		
		double area_ijk = 0.0;
		double alfa_ij = 0.0;
		double beta_ij = 0.0;	
		double theta_kij = 0.0;
		double theta_sum = 0.0;
		
		int neighnum = 0;
		double A_mixed = 0.0;
		array_1d<double,3> K_xi = ZeroVector(3);	
		double Delta_xi = 0.0;
		double Sp = 0.0;
				 
		WeakPointerVector< Condition >& neighb_faces = im->GetValue(NEIGHBOUR_CONDITIONS);
		//Loop over faces -> 1-ring neighborhood
		for (unsigned int i = 0; i < neighb_faces.size(); i++)
		{
		  int num_flag_var = 0;
		  num_flag_var += neighb_faces[i].GetGeometry()[0].FastGetSolutionStepValue(FLAG_VARIABLE);
		  num_flag_var += neighb_faces[i].GetGeometry()[1].FastGetSolutionStepValue(FLAG_VARIABLE);
		  num_flag_var += neighb_faces[i].GetGeometry()[2].FastGetSolutionStepValue(FLAG_VARIABLE);
		  if (num_flag_var > 2.999)
		  {
		    neighnum = 0;
		    for (unsigned int j = 0; j < neighb_faces[i].GetGeometry().size() ; j++)
		    {
		      if (neighb_faces[i].GetGeometry()[j].X() != xi || 
			neighb_faces[i].GetGeometry()[j].Y() != yi || neighb_faces[i].GetGeometry()[j].Z() != zi)
		      {
 			if (neighnum == 0)
                       // if (neighnum < 1e-15)
			{
			  xj = neighb_faces[i].GetGeometry()[j].X();
			  yj = neighb_faces[i].GetGeometry()[j].Y();
			  zj = neighb_faces[i].GetGeometry()[j].Z();
			}
			else
			{
			  xk = neighb_faces[i].GetGeometry()[j].X();
			  yk = neighb_faces[i].GetGeometry()[j].Y();
			  zk = neighb_faces[i].GetGeometry()[j].Z();		      
			}
			neighnum++;
		      }
		    }
		    
		      Vector3D(xi,yi,zi,xj,yj,zj,rij);
		      Vector3D(xj,yj,zj,xi,yi,zi,rji);
		      Vector3D(xi,yi,zi,xk,yk,zk,rik);
		      Vector3D(xk,yk,zk,xi,yi,zi,rki);
		      Vector3D(xj,yj,zj,xk,yk,zk,rjk);
		      norm_ij = Norm3D(rij);
		      norm_ik = Norm3D(rik);
		      norm_jk = Norm3D(rjk);
		      alfa_ij = Angle2vecs3D(rji,rjk);
		      theta_kij = Angle2vecs3D(rij,rik);
		      beta_ij = pi - alfa_ij - theta_kij;
		      CrossProduct3D(rij, rik, cross_prod_ijk);
		      area_ijk = 0.5*Norm3D(cross_prod_ijk);
		      if (alfa_ij > hpi || beta_ij > hpi || theta_kij > hpi)
		      {
		      //Obtuse triangle!
			if (theta_kij > hpi)
			  A_mixed += 0.5*area_ijk*1000.0;
			else
			  A_mixed += 0.25*area_ijk*1000.0;
		      }
		      else
		      {
			A_mixed += 0.125*cotan(alfa_ij)*norm_ik*norm_ik*1000.0;
			A_mixed += 0.125*cotan(beta_ij)*norm_ij*norm_ij*1000.0;		      
		      }
		    //Mean curvature normal operator:
		      K_xi[0] += (cotan(alfa_ij)*rki[0] + cotan(beta_ij)*rji[0])*1000.0;
		      K_xi[1] += (cotan(alfa_ij)*rki[1] + cotan(beta_ij)*rji[1])*1000.0;
		      K_xi[2] += (cotan(alfa_ij)*rki[2] + cotan(beta_ij)*rji[2])*1000.0;		    
		    //Sum of angles in node i for gaussian curvature
		      theta_sum += theta_kij;
		      Sp += (0.5*area_ijk - 0.125*tan(hpi - theta_kij)*norm_jk*norm_jk);
		  }
		}
		double A_mixed_norm = sqrt(A_mixed*A_mixed);
		if(A_mixed_norm < 1e-10)
		{
		  K_xi[0] = 0.0;
		  K_xi[1] = 0.0;
		  K_xi[2] = 0.0;
		}
		else
		{
		  //Now we have A_voronoi for node i.
		  K_xi[0] = (0.5/A_mixed)*K_xi[0];
		  K_xi[1] = (0.5/A_mixed)*K_xi[1];
		  K_xi[2] = (0.5/A_mixed)*K_xi[2];
		}
// 		KRATOS_WATCH(A_mixed)
// 		double normKxi = Norm3D(K_xi);
		array_1d<double,3> K_xi_norm = ZeroVector(3);
		K_xi_norm[0] = K_xi[0];
		K_xi_norm[1] = K_xi[1];
		K_xi_norm[2] = K_xi[2];
		NormalizeVec3D(K_xi_norm);
		im->FastGetSolutionStepValue(NORMAL_GEOMETRIC_X) = K_xi_norm[0];
		im->FastGetSolutionStepValue(NORMAL_GEOMETRIC_Y) = K_xi_norm[1];
		im->FastGetSolutionStepValue(NORMAL_GEOMETRIC_Z) = K_xi_norm[2];
		
		//Final step: kappa_H and kappa_G (mean curvature kappa_H is half the magnitude (norm) of K_xi vector)
		im->FastGetSolutionStepValue(MEAN_CURVATURE_3D) = Norm3D(K_xi);
		im->FastGetSolutionStepValue(GAUSSIAN_CURVATURE) = (2*pi - theta_sum)/Sp;
		
		//Principal curvatures, taking care of square root term:
		Delta_xi = ((im->FastGetSolutionStepValue(MEAN_CURVATURE_3D))*(im->FastGetSolutionStepValue(MEAN_CURVATURE_3D)) - im->FastGetSolutionStepValue(GAUSSIAN_CURVATURE));
		if(Delta_xi > 0.0)
		{
		  im->FastGetSolutionStepValue(PRINCIPAL_CURVATURE_1) = im->FastGetSolutionStepValue(MEAN_CURVATURE_3D) + sqrt(Delta_xi);
		  im->FastGetSolutionStepValue(PRINCIPAL_CURVATURE_2) = im->FastGetSolutionStepValue(MEAN_CURVATURE_3D) - sqrt(Delta_xi);
		}
		else
		{
		  im->FastGetSolutionStepValue(PRINCIPAL_CURVATURE_1) = im->FastGetSolutionStepValue(MEAN_CURVATURE_3D);
		  im->FastGetSolutionStepValue(PRINCIPAL_CURVATURE_2) = im->FastGetSolutionStepValue(MEAN_CURVATURE_3D);
		}
		
		A_mixed = 0.0;
		theta_sum = 0.0;
	      }
	      
	      //Now we correct curvature at triple points
	      if (im->FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		double neigh_fs = 0.0;
		double mean_curv = 0.0;
		WeakPointerVector< Node<3> >& neighb = im->GetValue(NEIGHBOUR_NODES);
		for (unsigned int i = 0; i < neighb.size(); i++)
		{
		  if (neighb[i].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
		  {
		    neigh_fs++;
		    mean_curv += neighb[i].FastGetSolutionStepValue(MEAN_CURVATURE_3D);
		  }
		}
		if(neigh_fs != 0.0)
		  im->FastGetSolutionStepValue(MEAN_CURVATURE_3D) = mean_curv/neigh_fs;
	      }
	      
 	      if ((im->FastGetSolutionStepValue(IS_STRUCTURE) != 0.0) && (im->FastGetSolutionStepValue(TRIPLE_POINT) == 0.0))
             // if ((im->FastGetSolutionStepValue(IS_STRUCTURE) != 0.0) && (im->FastGetSolutionStepValue(TRIPLE_POINT) < 1e-15))
		im->FastGetSolutionStepValue(MEAN_CURVATURE_3D) = 0.0;
	    }
	KRATOS_CATCH("")
    }
    
    void CalculateCurvatureContactLine(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ;
	    im != ThisModelPart.NodesEnd() ; ++im)
	    {
	      if (im->FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {		
		double xi = im->X();
		double yi = im->Y();
		double xj = 0.0;
		double yj = 0.0;
		double xk = 0.0;
		double yk = 0.0;
		
		int neighnum = 0;
		WeakPointerVector< Node<3> >& neighb = im->GetValue(NEIGHBOUR_NODES);
		for (unsigned int i = 0; i < neighb.size(); i++)
		{
		    if (neighb[i].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
		    {
 		      if (neighnum == 0)
                      //if (neighnum < 1e-15)
		      {
			xj = neighb[i].X();
			yj = neighb[i].Y();
		      }
		      else
		      {
			xk = neighb[i].X();
			yk = neighb[i].Y();
		      }
		      neighnum++;
		    }
		}
		
		double curv = CalculateCurv(xi,yi,xj,yj,xk,yk);
		if (neighnum > 2)
		    curv = 0.0;
		
		// Sign of the curvature
		double sign_curv = 1.0;
		double xc = 0.5*(xj + xk);
		double yc = 0.5*(yj + yk);
		array_1d<double,2> ric = ZeroVector(2);
		
		Vector2D(xi,yi,xc,yc,ric);
		
		array_1d<double,2> ni = ZeroVector(2);
		ni[0] = im->FastGetSolutionStepValue(NORMAL_GEOMETRIC_X);
		ni[1] = im->FastGetSolutionStepValue(NORMAL_GEOMETRIC_Y);
		double alfa = DotProduct2D(ric,ni);
		
		if (alfa > 0.0)
		    sign_curv = -1.0;
		
		im->FastGetSolutionStepValue(MEAN_CURVATURE_2D) = sign_curv*curv;
	      }		
	    }
	    
	KRATOS_CATCH("")
    }    
    
    
    void CalculatePrincipalDirections3D(ModelPart& ThisModelPart)
    {
	KRATOS_TRY
	
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ;
	    im != ThisModelPart.NodesEnd() ; ++im)
	    {
	      if (im->FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
	      {
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		// THERE IS A FUNCTION CALLED "FIND_NODAL_NEIGHBOURS_SURFACE_PROCESS" THAT THEORETICALLY
		// CAN OBTAIN NEIGHBOUR_NODES FROM NODES AT SURFACE. Now, for every node i in the surface,
		// we found its tangent plane and we take first neighbour j to compute rij -> first basis vector
		// Second basis vector is obtained by CrossProduct3D(n,rij)
		// Then we solve system 2x2 to find coordinates of rij expressed in base (u,v). Next step
		// is to solve 3x3 system of least squares and we obtain matrix B = (a b ; b c)
		// Find eigenvectors using analytical solution, normalize them and these are principal directions!
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		
		double xi = im->X();
		double yi = im->Y();
		double zi = im->Z();
		double xj = 0.0;
		double yj = 0.0;
		double zj = 0.0;
		double alfa = 0.0;
		double beta = 0.0;
		double kappaN_ij = 0.0;
		//double M = 0.0;
		//double N = 0.0;
		//double L = 0.0;
		array_1d<double,6> terms_func = ZeroVector(6);
		array_1d<double,10> terms_func_der = ZeroVector(10);
		WeakPointerVector< Node<3> >& neighb = im->GetValue(NEIGHBOUR_NODES);
		array_1d<double,3> dij = ZeroVector(3);
		
		double kappa_H = im->FastGetSolutionStepValue(MEAN_CURVATURE_3D);
		double kappa_G = im->FastGetSolutionStepValue(GAUSSIAN_CURVATURE);
		
		int option = 2; //choose 1 if you want to consider a + b = 2*kappa_H, choose 2 if a + c = 2*kappa_H
		
// 		boost::numeric::ublas::bounded_matrix<double, 3, 3> LHSmat = ZeroMatrix(3,3);
// 		array_1d<double,3> RHSvec = ZeroVector(3);
		boost::numeric::ublas::bounded_matrix<double, 2, 2> LHSmat = ZeroMatrix(2,2);
		array_1d<double,2> RHSvec = ZeroVector(2);		
		boost::numeric::ublas::bounded_matrix<double, 2, 2> B = ZeroMatrix(2,2);
		
		//STEP 1: find (u,v) basis of tangent plane at node i
		//Coordinates to create vector u (basis vector of tangent plane)
		double xu = neighb[0].X();
		double yu = neighb[0].Y();
		double zu = neighb[0].Z();
		array_1d<double,3> u_vec = ZeroVector(3);
		array_1d<double,3> v_vec = ZeroVector(3);
		array_1d<double,3> n_vec = ZeroVector(3);
		array_1d<double,2> eigen1 = ZeroVector(2);
		array_1d<double,2> eigen2 = ZeroVector(2);
		array_1d<double,3> eigen1_xyz = ZeroVector(3);
		array_1d<double,3> eigen2_xyz = ZeroVector(3);
		
		n_vec[0] = im->FastGetSolutionStepValue(NORMAL_GEOMETRIC_X);
		n_vec[1] = im->FastGetSolutionStepValue(NORMAL_GEOMETRIC_Y);
		n_vec[2] = im->FastGetSolutionStepValue(NORMAL_GEOMETRIC_Z);
		
		Vector3D(xi,yi,zi,xu,yu,zu,u_vec);
		FindUnitTangent3D(u_vec,n_vec);
		
		//The other basis vector v is just cross product between n and u
		CrossProduct3D(n_vec,u_vec,v_vec);
		NormalizeVec3D(v_vec);
		
		//STEP 2: For each neighbour node, find components of dij in base (u,v) and add term to system matrix
		int num_neighbs = 0;
		for (unsigned j = 0; j < neighb.size(); j++)
		{
		  num_neighbs++;
		}
		
		for (unsigned j = 0; j < neighb.size(); j++)
		{
                
                    double xj;
                    double yj;
                    double zj;
		  //first find unit tanjent of edge i-j
		  xj = neighb[j].X();
		  yj = neighb[j].Y();
		  zj = neighb[j].Z();
		  Vector3D(xi,yi,zi,xj,yj,zj,dij);
		  FindUnitTangent3D(dij, n_vec);
		  
		  //now find components in (u,v) base
		  FindComponentsUV(alfa,beta,dij,u_vec,v_vec);
		  
		  //terms are added to global system matrix and RHSvec
		  kappaN_ij = NormalCurvature3D(xi,yi,zi,xj,yj,zj,n_vec);
// 		  AddTermsToSys(kappaN_ij,kappa_H,alfa,beta,num_neighbs,LHSmat,RHSvec,option); //Includes versions depending on option below
// 		  AddTermsToEq(kappaN_ij,kappa_H,alfa,beta,num_neighbs,M,N,L);
		  AddTermsToEq2(kappaN_ij,kappa_H,alfa,beta,num_neighbs,terms_func,terms_func_der);
		}
		
		double a = 0.0;
		double b = 0.0;
		double c = 0.0;
		
		//STEP 3: global system matrix is full. Now we solve the system
		if(option == 1)
		{
		  //OPTION 1.1 - Meyer is right AND consider just first condition
		  SolveSys2x2(b,c,LHSmat,RHSvec);
		  a = 2.0*kappa_H - c;
		}
		else
		{
		  //OPTION 2.1 - consider just first condition
// 		  SolveSys3x3(a,b,c,LHSmat,RHSvec);
// 		  SolveSys2x2(a,b,LHSmat,RHSvec);
// 		  c = 2.0*kappa_H - a;
// 		  b = sqrt(a*(2.0*kappa_H - a) - kappa_G);
		  //OPTION 2.2 - consider both conditions
                NewtonMethod(terms_func,terms_func_der,a,kappa_H,kappa_G,kappaN_ij);
                c = 2.0*kappa_H - a;
                b = sqrt(a*(2.0*kappa_H - a) - kappa_G);
		}
		
		//Regardless option choice, fill the curvature tensor
		B(0,0) = a;
		B(1,0) = b;
		B(0,1) = b;
		B(1,1) = c;
		KRATOS_WATCH(option)
		KRATOS_WATCH(a)
		KRATOS_WATCH(b)
		KRATOS_WATCH(c)
		
		double cond_sum = 0.0;
		double cond_prod = a*c - b*b;
		double eps_cond = 1e-5;
		if (option == 1)
		  cond_sum = 0.5*(a + b);
		else
		  cond_sum = 0.5*(a + c);
		if ((cond_sum - kappa_H) > eps_cond || cond_prod != kappa_G)
		{
		  KRATOS_WATCH("Not fulfilling this condition!!!!!!!!!!!!!!!!")
		  KRATOS_WATCH(cond_prod)
		  KRATOS_WATCH(kappa_G)
		}
		
		//STEP 4: find eigenvectors of matrix B
		double T = a + c;
		double D = a*c - b*b;
		
		double lambda_1 = 0.5*T + sqrt(0.25*T*T - D);
		double lambda_2 = 0.5*T - sqrt(0.25*T*T - D);
		
		if (b != 0.0)
		{
// 		  eigen1[0] = lambda_1 - c;
// 		  eigen1[1] = b;
// 		  eigen2[0] = lambda_2 - c;
// 		  eigen2[1] = b;
		  //OR:
		  eigen1[0] = b;
		  eigen1[1] = lambda_1 - a;
		  eigen2[0] = b;
		  eigen2[1] = lambda_2 - a;		  
		}
		else
		{
		  KRATOS_WATCH("b is zero")
		  eigen1[0] = 1.0;
		  eigen1[1] = 0.0;
		  eigen2[0] = 0.0;
		  eigen2[1] = 1.0;
		}
		
		//Watch out! Eigenvectors are expressed in terms of base (u,v) -> transform to (x,y,z)
		eigen1_xyz[0] = eigen1[0]*u_vec[0] + eigen1[1]*v_vec[0];
		eigen1_xyz[1] = eigen1[0]*u_vec[1] + eigen1[1]*v_vec[1];
		eigen1_xyz[2] = eigen1[0]*u_vec[2] + eigen1[1]*v_vec[2];
		eigen2_xyz[0] = eigen2[0]*u_vec[0] + eigen2[1]*v_vec[0];
		eigen2_xyz[1] = eigen2[0]*u_vec[1] + eigen2[1]*v_vec[1];
		eigen2_xyz[2] = eigen2[0]*u_vec[2] + eigen2[1]*v_vec[2];
		NormalizeVec3D(eigen1_xyz);
		NormalizeVec3D(eigen2_xyz);
		im->FastGetSolutionStepValue(PRINCIPAL_DIRECTION_1_X) = eigen1_xyz[0];
		im->FastGetSolutionStepValue(PRINCIPAL_DIRECTION_1_Y) = eigen1_xyz[1];
		im->FastGetSolutionStepValue(PRINCIPAL_DIRECTION_1_Z) = eigen1_xyz[2];
		im->FastGetSolutionStepValue(PRINCIPAL_DIRECTION_2_X) = eigen2_xyz[0];
		im->FastGetSolutionStepValue(PRINCIPAL_DIRECTION_2_Y) = eigen2_xyz[1];
		im->FastGetSolutionStepValue(PRINCIPAL_DIRECTION_2_Z) = eigen2_xyz[2];
	      }
	    }
	KRATOS_CATCH("")
    }    
    
    double CalculateVol(const double x0, const double y0, const double x1, const double y1, const double x2, const double y2)
    {
	return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
    }
    
    double CalculateCurv(const double x0, const double y0, const double x1, const double y1, const double x2, const double y2)
    {
      array_1d<double,2> r01 = ZeroVector(2);
      array_1d<double,2> r20 = ZeroVector(2); 
      Vector2D(x1,y1,x0,y0,r01);
      Vector2D(x0,y0,x2,y2,r20);
      double norm01 = Norm2D(r01);
      double norm20 = Norm2D(r20);
      NormalizeVec2D(r01);
      NormalizeVec2D(r20);
      array_1d<double,2> r_diff = r01 - r20;
      double norm_diff = Norm2D(r_diff);
      return 2.0*norm_diff/(norm01 + norm20);
    }
    
    void IsObtuse(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1,
		  const double x2, const double y2, const double z2, bool& isobt, bool& isobt_i, double& alfa)
    {
      array_1d<double,3> r01 = ZeroVector(3);
      array_1d<double,3> r02 = ZeroVector(3);
      array_1d<double,3> r12 = ZeroVector(3);
      Vector3D(x0,y0,z0,x1,y1,z1,r01);
      Vector3D(x0,y0,z0,x2,y2,z2,r02);
      Vector3D(x1,y1,z1,x2,y2,z2,r12);      
      
      double hpi = 3.14159265*0.5;
      double alfa0 = 0.0; //angle at node 0
      double alfa1 = 0.0; //angle at node 1
      double alfa2 = 0.0; //angle at node 2
      
      alfa0 = Angle2vecs3D(r01,r02);
      r01[0] = -r01[0];
      r01[1] = -r01[1];
      r01[2] = -r01[2];
      alfa1 = Angle2vecs3D(r01,r12); //minus sign because alfa1 is formed by vectors r10 and r12
      alfa2 = Angle2vecs3D(r02,r12);
      
      alfa = alfa2;
      
      if (alfa0 > hpi || alfa1 > hpi || alfa2 > hpi)
	  isobt = true;
      if (alfa0 > hpi)
	  isobt_i = true;
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
    
    double Angle2vecs2D(const array_1d<double,2>& a, const array_1d<double,2>& b)
    {
      double norm_a = Norm2D(a);
      double norm_b = Norm2D(b);
      double temp = 0.0;
      if (norm_a*norm_b > -1.0e-20 && norm_a*norm_b < 1.0e-20)
	temp = 0.0;
      else
	temp = DotProduct2D(a,b)/(norm_a*norm_b);
      return acos(temp);
    }    
    
    double Angle2vecs3D(const array_1d<double,3>& a, const array_1d<double,3>& b)
    {
      double norm_a = Norm3D(a);
      double norm_b = Norm3D(b);
      double temp = 0.0;
//       if (norm_a*norm_b > -1.0e-20 && norm_a*norm_b < 1.0e-20)
      if (norm_a*norm_b == 0.0)
      //if (norm_a*norm_b < 1e-15)
	temp = 0.0;
      else
	temp = DotProduct3D(a,b)/(norm_a*norm_b);
      return acos(temp);
    }
    
    double cotan(const double& alfa)
    {
      return cos(alfa)/sin(alfa);
    }
    
    double NormalCurvature3D(const double xi, const double yi, const double zi,
		  const double xj, const double yj, const double zj, array_1d<double,3>& n)
    {
      double kappaN_ij = 0.0;
      array_1d<double,3> rji = ZeroVector(3);
      Vector3D(xj,yj,zj,xi,yi,zi,rji);
      double norm_rji = Norm3D(rji);
      
      kappaN_ij = DotProduct3D(rji,n);
      kappaN_ij = 2*kappaN_ij/(norm_rji*norm_rji);
      return kappaN_ij;
    }
    
    void FindUnitTangent3D(array_1d<double,3>& u_vec, array_1d<double,3>& n_vec)
    {
      double dotprod = DotProduct3D(u_vec,n_vec);
      u_vec[0] = (u_vec[0] - dotprod*n_vec[0]);
      u_vec[1] = (u_vec[1] - dotprod*n_vec[1]);
      u_vec[2] = (u_vec[2] - dotprod*n_vec[2]);
      double norm_uvec = Norm3D(u_vec);
      u_vec[0] /= norm_uvec;
      u_vec[1] /= norm_uvec;
      u_vec[2] /= norm_uvec;      
    }
    
    void FindComponentsUV(double& alfa, double& beta, array_1d<double,3>& dij, array_1d<double,3>& u_vec, array_1d<double,3>& v_vec)
    {
      boost::numeric::ublas::bounded_matrix<double, 2, 2> mat_syst;
      boost::numeric::ublas::bounded_matrix<double, 2, 2> mat_dx;
      boost::numeric::ublas::bounded_matrix<double, 2, 2> mat_dy;
      int num_case = 0;
      mat_syst(0,0) = u_vec[0];
      mat_syst(1,0) = u_vec[1];
      mat_syst(0,1) = v_vec[0];
      mat_syst(1,1) = v_vec[1];      
      double det_matsys = Det2x2(mat_syst);
      if (det_matsys*det_matsys < 0.00000001)
      {
	mat_syst(1,0) = u_vec[2];
	mat_syst(1,1) = v_vec[2];
	det_matsys = Det2x2(mat_syst);
	num_case++;
	if (det_matsys*det_matsys < 0.00000001)
	{
	  mat_syst(0,0) = u_vec[1];
	  mat_syst(0,1) = v_vec[1];
	  num_case++;
	}
      }
      
      if (num_case == 0)
      //if (num_case < 1e-15)
      {
	mat_dx(0,0) = dij[0];
	mat_dx(1,0) = dij[1];
	mat_dx(0,1) = v_vec[0];
	mat_dx(1,1) = v_vec[1];
	mat_dy(0,0) = u_vec[0];
	mat_dy(1,0) = u_vec[1];
	mat_dy(0,1) = dij[0];
	mat_dy(1,1) = dij[1];
      }
      if (num_case == 1)
      {
	mat_dx(0,0) = dij[0];
	mat_dx(1,0) = dij[2];
	mat_dx(0,1) = v_vec[0];
	mat_dx(1,1) = v_vec[2];
	mat_dy(0,0) = u_vec[0];
	mat_dy(1,0) = u_vec[2];
	mat_dy(0,1) = dij[0];
	mat_dy(1,1) = dij[2];
      }
      else
      {
	mat_dx(0,0) = dij[1];
	mat_dx(1,0) = dij[2];
	mat_dx(0,1) = v_vec[1];
	mat_dx(1,1) = v_vec[2];
	mat_dy(0,0) = u_vec[1];
	mat_dy(1,0) = u_vec[2];
	mat_dy(0,1) = dij[1];
	mat_dy(1,1) = dij[2];
      }
      
      double det_matdx = Det2x2(mat_dx);
      double det_matdy = Det2x2(mat_dy);
      
      alfa = det_matdx/det_matsys;
      beta = det_matdy/det_matsys;
    }
    
    void SolveSys2x2(double& a, double& b,boost::numeric::ublas::bounded_matrix<double, 2, 2>& LHSmat,
		     array_1d<double,2>& RHSvec)
    {
      boost::numeric::ublas::bounded_matrix<double, 2, 2> mat_syst = ZeroMatrix(2,2);
      boost::numeric::ublas::bounded_matrix<double, 2, 2> mat_dx = ZeroMatrix(2,2);
      boost::numeric::ublas::bounded_matrix<double, 2, 2> mat_dy = ZeroMatrix(2,2);
      
      CopyMatrix2x2(LHSmat,mat_syst);
      CopyMatrix2x2(LHSmat,mat_dx);
      CopyMatrix2x2(LHSmat,mat_dy);
      double det_matsys = Det2x2(mat_syst);
      if (det_matsys > -0.000000001 && det_matsys < 0.000000001)
	  det_matsys = 0.000000001;

      mat_dx(0,0) = RHSvec[0];
      mat_dx(1,0) = RHSvec[1];
      mat_dy(0,1) = RHSvec[0];
      mat_dy(1,1) = RHSvec[1];
      
      double det_matdx = Det2x2(mat_dx);
      double det_matdy = Det2x2(mat_dy);
      
      a = det_matdx/det_matsys;
      b = det_matdy/det_matsys;
    }
    
    void SolveSys3x3(double& a, double& b, double& c,boost::numeric::ublas::bounded_matrix<double, 3, 3>& LHSmat,
		     array_1d<double,3>& RHSvec)
    {
      boost::numeric::ublas::bounded_matrix<double, 3, 3> mat_syst = ZeroMatrix(3,3);
      boost::numeric::ublas::bounded_matrix<double, 3, 3> mat_dx = ZeroMatrix(3,3);
      boost::numeric::ublas::bounded_matrix<double, 3, 3> mat_dy = ZeroMatrix(3,3);
      boost::numeric::ublas::bounded_matrix<double, 3, 3> mat_dz = ZeroMatrix(3,3);
      
      CopyMatrix3x3(LHSmat,mat_syst);
      CopyMatrix3x3(LHSmat,mat_dx);
      CopyMatrix3x3(LHSmat,mat_dy);
      CopyMatrix3x3(LHSmat,mat_dz);
      double det_matsys = Det3x3(mat_syst);

      mat_dx(0,0) = RHSvec[0];
      mat_dx(1,0) = RHSvec[1];
      mat_dx(2,0) = RHSvec[2];
      mat_dy(0,1) = RHSvec[0];
      mat_dy(1,1) = RHSvec[1];
      mat_dy(2,1) = RHSvec[2];
      mat_dz(0,2) = RHSvec[0];
      mat_dz(1,2) = RHSvec[1];
      mat_dz(2,2) = RHSvec[2];
      
      double det_matdx = Det3x3(mat_dx);
      double det_matdy = Det3x3(mat_dy);
      double det_matdz = Det3x3(mat_dz);
      
      a = det_matdx/det_matsys;
      b = det_matdy/det_matsys;
      c = det_matdz/det_matsys;
    }    
    
//     void AddTermsToSys(double& kappaN_ij, double& alfa, double& beta, int& num_neighbs,
// 		       boost::numeric::ublas::bounded_matrix<double, 3, 3>& LHSmat, array_1d<double,3>& RHSvec)
    void AddTermsToSys(double& kappaN_ij, double& kappa_H, double& alfa, double& beta, int& num_neighbs,
		       boost::numeric::ublas::bounded_matrix<double, 2, 2>& LHSmat, array_1d<double,2>& RHSvec, const int& option)
    {
      double weight = 1.0/num_neighbs;
      /*
      LHSmat(0,0) += weight*alfa*alfa*alfa*alfa;
      LHSmat(1,0) += weight*alfa*alfa*alfa*beta;
      LHSmat(2,0) += weight*alfa*alfa*beta*beta;
      LHSmat(0,1) += 2.0*weight*alfa*alfa*alfa*beta;
      LHSmat(1,1) += 2.0*weight*alfa*alfa*beta*beta;
      LHSmat(2,1) += 2.0*weight*alfa*beta*beta*beta;
      LHSmat(0,2) += weight*alfa*alfa*beta*beta;
      LHSmat(1,2) += weight*alfa*beta*beta*beta;
      LHSmat(2,2) += weight*beta*beta*beta*beta;
      RHSvec[0] += weight*kappaN_ij*alfa*alfa;
      RHSvec[1] += weight*kappaN_ij*alfa*beta;
      RHSvec[2] += weight*kappaN_ij*beta*beta;
      */
      /*
      LHSmat(0,0) += weight*(alfa*alfa - 2.0*alfa*beta)*(alfa*alfa - 2.0*alfa*beta);
      LHSmat(0,1) += weight*(alfa*alfa - 2.0*alfa*beta)*beta*beta;
      LHSmat(1,0) += weight*(alfa*alfa - 2.0*alfa*beta)*beta*beta;
      LHSmat(1,1) += weight*beta*beta*beta*beta;
      RHSvec[0] += weight*(kappaN_ij - 4.0*kappa_H*alfa*beta)*(alfa*alfa - 2.0*alfa*beta);
      RHSvec[1] += weight*(kappaN_ij - 4.0*kappa_H*alfa*beta)*beta*beta;
      */
      if(option == 1)
      {
	//OPTION 1.1 - Meyer is right
	LHSmat(0,0) += weight*(-alfa*alfa - 2.0*alfa*beta)*(-alfa*alfa - 2.0*alfa*beta);
	LHSmat(0,1) += weight*beta*beta*(-alfa*alfa - 2.0*alfa*beta);
	LHSmat(1,0) += weight*(-alfa*alfa - 2.0*alfa*beta)*beta*beta;
	LHSmat(1,1) += weight*beta*beta*beta*beta;
	RHSvec[0] += weight*(kappaN_ij - 2.0*kappa_H*alfa*alfa)*(-alfa*alfa - 2.0*alfa*beta);
	RHSvec[1] += weight*(kappaN_ij - 2.0*kappa_H*alfa*alfa)*beta*beta;      
      }
      else
      {
	//OPTION 2.1 - Meyer is wrong
	LHSmat(0,0) += weight*(alfa*alfa - beta*beta)*(alfa*alfa - beta*beta);
	LHSmat(0,1) += 2.0*weight*alfa*beta*(alfa*alfa - beta*beta);
	LHSmat(1,0) += weight*alfa*beta*(alfa*alfa - beta*beta);
	LHSmat(1,1) += 2.0*weight*alfa*alfa*beta*beta;
	RHSvec[0] += weight*(kappaN_ij - 2.0*kappa_H*beta*beta)*(alfa*alfa - beta*beta);
	RHSvec[1] += weight*(kappaN_ij - 2.0*kappa_H*beta*beta)*alfa*beta;      
      }
    }
    /*
    void AddTermsToEq(double& kappaN_ij, double& kappa_H, double& alfa, double& beta, int& num_neighbs,
		      double& term_1, double& term_2, double& term_3, double& term_4, double& term_5, double& term_6)
    {
      double weight = 1.0/num_neighbs;
      M += weight*(alfa*alfa - beta*beta);
      N += weight*2.0*alfa*beta;
      L += weight*(kappaN_ij - 2.0*beta*beta*kappa_H);
    }
    */
    void AddTermsToEq2(double& kappaN_ij, double& kappa_H, double& alfa, double& beta, int& num_neighbs,
		      array_1d<double,6>& terms_vec, array_1d<double,10>& terms_vec_der)
    {
      double weight = 1.0/num_neighbs;
      terms_vec[0] += weight*(alfa*alfa - beta*beta)*(alfa*alfa - beta*beta);
      terms_vec[1] += weight*alfa*beta*(alfa*alfa - beta*beta);
      terms_vec[2] += 2.0*weight*alfa*beta*(alfa*alfa - beta*beta);
      terms_vec[3] += 2.0*weight*alfa*alfa*beta*beta;
      terms_vec[4] += weight*(2.0*kappa_H*beta*beta - kappaN_ij)*(alfa*alfa - beta*beta);
      terms_vec[5] += weight*weight*(2.0*kappa_H*beta*beta - kappaN_ij)*alfa*beta;
      
      terms_vec_der[0] += weight*(alfa*alfa - beta*beta)*(alfa*alfa - beta*beta);
      terms_vec_der[1] += weight*(alfa*alfa - beta*beta)*(alfa*alfa - beta*beta);
      terms_vec_der[2] += 2.0*weight*alfa*beta*(alfa*alfa - beta*beta);
      terms_vec_der[3] += weight*alfa*beta*(alfa*alfa - beta*beta);
      terms_vec_der[4] += 0.5*weight*alfa*beta*(alfa*alfa - beta*beta);
      terms_vec_der[5] += 3.0*weight*alfa*beta*(alfa*alfa - beta*beta);
      terms_vec_der[6] += 4.0*weight*alfa*alfa*beta*beta;
      terms_vec_der[7] += 2.0*weight*alfa*alfa*beta*beta;
      terms_vec_der[8] += 0.5*weight*(2.0*kappa_H*beta*beta - kappaN_ij)*(alfa*alfa - beta*beta)*(alfa*alfa - beta*beta);
      terms_vec_der[9] += 2.0*weight*(2.0*kappa_H*beta*beta - kappaN_ij)*alfa*beta;      
    }
    
    void SolveEq2ndDeg(double& A, double& B, double& C, double& a)
    {
      a = (-B - sqrt(B*B - 4.0*A*C))/(2.0*A);
    }
    
    void NewtonMethod(array_1d<double,6>& terms_vec, array_1d<double,10>& terms_vec_der, double& a,
		      double& kappa_H, double& kappa_G, double& kappaN_ij)
    {
      double x0 = 1.0;			//initial guess
      double x1 = 0.0;			//next guess or solution
      double fx = 0.0;			//function
      double dfx = 0.0;      		//function derivative
      double tol = 1.0e-7;		//tolerance for the solution
      double epsi = 1.0e-10;		//minimum value of function derivative
      unsigned int MaxIter = 20;	//maximum number of iterations
      bool foundsol = false;		//solution has been found
      for(unsigned int i = 0; i < MaxIter; i++)
      {
          
          double fx = 0.0;
          double dfx = 0.0;

	  fx = func_Newton(x0,terms_vec,kappa_H,kappa_G,kappaN_ij);
	  dfx = dxfunc_Newton(x0,terms_vec_der,kappa_H,kappa_G,kappaN_ij);
	  if (abs(dfx) < epsi)
	    x0 += 1.0;
	  else
	    x1 = x0 - fx/dfx;
	  
	  if(abs(x1-x0)/abs(x0) < tol)
	    foundsol = true;
	    
	  x1 = x0;
      }
      if (foundsol == false)
	KRATOS_WATCH("did not converge!!!!!!!!!!!!!!!!!!!")
    }
    
    double func_Newton(double& x, array_1d<double,6>& terms_vec, double& kappa_H, double& kappa_G, double& kappaN_ij)
    {
      double f_x;
      f_x = terms_vec[0]*x*f_a(x,kappa_H,kappa_G) + terms_vec[1]*2.0*(kappa_H - x)*x*sqrt(f_a(x,kappa_H,kappa_G)) + 
	    terms_vec[2]*f_a(x,kappa_H,kappa_G)*sqrt(f_a(x,kappa_H,kappa_G)) + terms_vec[3]*f_a(x,kappa_H,kappa_G) + 
	    terms_vec[4]*sqrt(f_a(x,kappa_H,kappa_G)) + terms_vec[5]*2.0*(kappa_H - a);
      return f_x;
    }
    
    double dxfunc_Newton(double& x, array_1d<double,10>& tvder, double& kappa_H, double& kappa_G, double& kappaN_ij)
    {
      double df_x;
      df_x = tvder[0]*f_a(x,kappa_H,kappa_G) + tvder[1]*x*df_a(x,kappa_H) - tvder[2]*x*sqrt(f_a(x,kappa_H,kappa_G)) + 
	      tvder[3]*2.0*(kappa_H - x)*sqrt(f_a(x,kappa_H,kappa_G)) + 
	      tvder[4]*2.0*(kappa_H - x)*x*0.5*df_a(x,kappa_H)/sqrt(f_a(x,kappa_H,kappa_G)) + 
	      tvder[5]*df_a(x,kappa_H)*sqrt(f_a(x,kappa_H,kappa_G)) + tvder[6]*f_a(x,kappa_H,kappa_G) + 
	      tvder[7]*2.0*(kappa_H - x)*df_a(x,kappa_H) + tvder[8]*sqrt(f_a(x,kappa_H,kappa_G))*df_a(x,kappa_H) - tvder[9];
      return df_x;
    }    
    
    double f_a(double& x, double& kappa_H, double& kappa_G)
    {
      return (x*(2.0*kappa_H - x) - kappa_G);
    }
    
    double df_a(double& x, double& kappa_H)
    {
      return (2.0*(kappa_H - x));
    }
    
    double Det2x2(boost::numeric::ublas::bounded_matrix<double, 2, 2>& input)
    {
      return (input(0,0)*input(1,1) - input(0,1)*input(1,0));
    }
    
    double Det3x3(boost::numeric::ublas::bounded_matrix<double, 3, 3>& input)
    {
      return (input(0,0)*input(1,1) - input(0,1)*input(1,0));
    }
    
    void NormalizeVec2D(array_1d<double,2>& input)
    {
      double norm = Norm2D(input);
      if (norm != 0.0)
      {
	input[0] /= norm;
	input[1] /= norm;
      }
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
    
    void CopyMatrix3x3(boost::numeric::ublas::bounded_matrix<double, 3, 3>& input,
		       boost::numeric::ublas::bounded_matrix<double, 3, 3>& output)
    {
      for (unsigned int i = 0; i < 3; i++)
      {
	for (unsigned int j = 0; j < 3; j++)
	{
	  output(i,j) = input(i,j);
	}
      }
    }
    
    void CopyMatrix2x2(boost::numeric::ublas::bounded_matrix<double, 2, 2>& input,
		       boost::numeric::ublas::bounded_matrix<double, 2, 2>& output)
    {
      for (unsigned int i = 0; i < 2; i++)
      {
	for (unsigned int j = 0; j < 2; j++)
	{
	  output(i,j) = input(i,j);
	}
      }
    }
    
    bool GetNeighbours(array_1d<double,2> n_node,WeakPointerVector< Node<3> >& neighb,double& x1,double& y1,double& x2,double& y2, int& i_neck, int& neighnum)
    {
      bool neck = false;
      double theta = 0.0;
      array_1d<double,2> n_vec;
      for (unsigned int i = 0; i < neighb.size(); i++)
      {
	n_vec = neighb[i].FastGetSolutionStepValue(NORMAL);
	if (neighb[i].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)// || neighb[i].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0 || neighb[i].FastGetSolutionStepValue(IS_STRUCTURE) != 0.0)
	{
	  theta = Angle2vecs2D(n_node,n_vec);
	  if (theta > 2.8 || theta < -2.8)
	  {
	    i_neck = i;
	    neck = true;
	  }
	  else
	  {
 	    if (neighnum == 0)
            //if (neighnum < 1e-15)
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
      }
      return neck;
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
    virtual std::string Info() const
    {
        return "CalculateCurvature";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "CalculateCurvature";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
//		CalculateCurvature& operator=(CalculateCurvature const& rOther);

    /// Copy constructor.
//		CalculateCurvature(CalculateCurvature const& rOther);


    ///@}

}; // Class CalculateCurvature

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateCurvature& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateCurvature& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



}  // namespace Kratos.


#endif // KRATOS_CREATE_INTERFACE_CONDITIONS_PROCESS_INCLUDED  defined 


