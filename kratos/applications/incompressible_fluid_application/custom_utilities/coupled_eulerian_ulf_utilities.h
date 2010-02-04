/*
==============================================================================
KratosULFApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2009-01-15 14:50:24 $
//   Revision:            $Revision: 1.12 $
//
//
//utilities for the updated lagrangian fluid

#if !defined(KRATOS_COUPLED_EUL_ULF_UTILITIES_INCLUDED )
#define  KRATOS_COUPLED_EUL_ULF_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "incompressible_fluid_application.h"

namespace Kratos
{
	class CoupledEulerianUlfUtils
	{
	public:
		typedef Node<3> NodeType;
		//**********************************************************************************************
		//**********************************************************************************************

		//function that stores the elements that represent the Lagrangian part within Eulerian mesh into a pasuedo_lag_part

		void SavePseudoLagPart(ModelPart& full_model_part, ModelPart& lagrangian_part, ModelPart& pseudo_lag_part)
		{
		KRATOS_TRY
		int n_int; 
//		int n_disabled;
		//clear reduced_model_part
		pseudo_lag_part.Conditions().clear();
		pseudo_lag_part.Elements().clear();
		pseudo_lag_part.Nodes().clear();

		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
			
			if (n_int==3)
				{
				pseudo_lag_part.AddElement(*(im.base()));
				}

		}
		
		for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; 
				in != full_model_part.NodesEnd() ; ++in)
		{	  
			
			double n_int=in->FastGetSolutionStepValue(IS_INTERFACE);
			
			if (n_int==1.0)
				{
				pseudo_lag_part.AddNode(*(in.base()));				
				}			

		}
		
		
		
		KRATOS_CATCH("")
		}
		//**********************************************************************************************
		//**********************************************************************************************
		// function that applies ADDITIONAL DIRICHLET CONDITIONS UPON THE nodes that lie at the INTERFACE between Lagrangain and Eulerian parts
		void ApplyProjDirichlet(ModelPart& full_model_part, ModelPart& aux_conditions_model_part)
		{
		KRATOS_TRY
		int n_int;

		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE,1);
			n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE,1);
			n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE,1);
			
			if (n_int==3)
				{
				im->GetGeometry()[0].FastGetSolutionStepValue(DISABLE)=false;
				im->GetGeometry()[1].FastGetSolutionStepValue(DISABLE)=false;
				im->GetGeometry()[2].FastGetSolutionStepValue(DISABLE)=false;
				}

		}
		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
			
			if (n_int==3)
				{
				im->GetGeometry()[0].FastGetSolutionStepValue(DISABLE)=true;
				im->GetGeometry()[1].FastGetSolutionStepValue(DISABLE)=true;
				im->GetGeometry()[2].FastGetSolutionStepValue(DISABLE)=true;
				}

		}

		
		
				
		
		
		//first we remove the Dirichlet conditions from the nodes that were defining the interface in the previous step:

		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE,1);
			n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE,1);
			n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE,1);
			
			if (n_int>0.0)
				{
								
				im->GetGeometry()[0].Free(VELOCITY_X);
				im->GetGeometry()[0].Free(VELOCITY_Y);

				im->GetGeometry()[1].Free(VELOCITY_X);
				im->GetGeometry()[1].Free(VELOCITY_Y);

				im->GetGeometry()[2].Free(VELOCITY_X);
				im->GetGeometry()[2].Free(VELOCITY_Y);

				im->GetGeometry()[0].Free(AUX_VEL_X);
				im->GetGeometry()[0].Free(AUX_VEL_Y);

				im->GetGeometry()[1].Free(AUX_VEL_X);
				im->GetGeometry()[1].Free(AUX_VEL_Y);

				im->GetGeometry()[2].Free(AUX_VEL_X);
				im->GetGeometry()[2].Free(AUX_VEL_Y);

				im->GetGeometry()[0].Free(PRESSURE);				
				im->GetGeometry()[1].Free(PRESSURE);
				im->GetGeometry()[2].Free(PRESSURE);
				}		

		}
		

		if (aux_conditions_model_part.Conditions().size()>0);
		{
			//and corerct the velocity of the elements that were "intersected" and stored and solved in aux_condition_model_part
			for(ModelPart::ConditionsContainerType::iterator ic = aux_conditions_model_part.ConditionsBegin() ; 
					ic != aux_conditions_model_part.ConditionsEnd() ; ++ic)
			{	  
			
				n_int=ic->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
				n_int+=ic->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
				n_int+=ic->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
				//elements with n_int=3 lie inside the fictitious domain, and therefore are of no interest
				if (n_int==1 || n_int==2)
			
					{
					//and now fix the velocity of the IS_INTERFACE nodes (those are lying inside of the fictitious domain)
					for (int i=0;i<2;i++)
						{
						if (ic->GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==1.0)
							{
							ic->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)=ic->GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL);
							//KRATOS_WATCH(im->GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL))
						
							//im->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)=0.0;

							ic->GetGeometry()[i].Fix(VELOCITY_X);
							ic->GetGeometry()[i].Fix(VELOCITY_Y);
							//im->GetGeometry()[i].Fix(PRESSURE);
						
							}
						}
					}
								
			
			

			}
		}
		
		
		KRATOS_CATCH("")
		}
		//**********************************************************************************************
		//**********************************************************************************************
		//origin - Lagrangian part, destination - Eulerian, AuxModelPart will store the "projected Dirichlet conditions"
		void FindInterface(ModelPart& origin_model_part, ModelPart& destination_model_part, ModelPart& AuxModelPart)
		{
			KRATOS_TRY  
			//at each step we need to make AuxModelPart empty
			AuxModelPart.Conditions().clear();
			AuxModelPart.Elements().clear();
			AuxModelPart.Nodes().clear();

			//FREE ALL THE AUX_VELS of the destination model part!
			for(ModelPart::NodesContainerType::iterator in = destination_model_part.NodesBegin() ; 
				in != destination_model_part.NodesEnd() ; ++in)
			{
				//reset the IS_INTERFACE to zero for all the nodes of interface conditions
				in->Free(AUX_VEL_X);
				in->Free(AUX_VEL_Y);
				in->FastGetSolutionStepValue(AUX_VEL_X)=0.0;				
				in->FastGetSolutionStepValue(AUX_VEL_Y)=0.0;

			}


			//first we need to find the biggest element size among all the elements of both the model parts
			// it will be set as the seacrh radius for the tree-search
			typedef Node<3> PointType;
			typedef Node<3>::Pointer PointPointerType;
			typedef std::vector<PointType::Pointer>           PointVector;
			typedef PointVector::iterator PointIterator;
			typedef std::vector<double>               DistanceVector;
			typedef std::vector<double>::iterator     DistanceIterator;


			

			// bucket types
			//typedef Bucket<3, PointType, ModelPart::NodesContainerType, PointPointerType, PointIterator, DistanceIterator > BucketType;
			//typedef Bins< 3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > StaticBins;
			// bucket types
			typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
		
				
			//*************
			// DynamicBins;	
			typedef Tree< KDTreePartition<BucketType> > kd_tree; //Kdtree;
			//typedef Tree< StaticBins > Bin; 			     //Binstree;
			unsigned int bucket_size = 50;
			
			//performing the interpolation - all of the nodes in this list will be preserved
			unsigned int max_results = 200;
			//PointerVector<PointType> res(max_results);
			//NodeIterator res(max_results);
			PointVector res(max_results);
			PointVector res1(max_results);
			DistanceVector res_distances(max_results);
			DistanceVector res_distances1(max_results);
			Node<3> work_point(0,0.0,0.0,0.0);
			array_1d<double, 3> temp;

			array_1d<double, 3> vel;
			array_1d<double, 3> vel1;
			array_1d<double, 3> vel2;


			double x0, y0, x1, y1, x2, y2, xc, yc, l0, l1, l2, l, radius;//, x0_orig, x1_orig, y0_orig, y1_orig;
			//PointVector IntersectionPoints;

			std::vector<array_1d<double,3> > IntersectionPoints;
			std::vector<array_1d<double,3> > IntersectionVel;
			//there should be just two intersection points per element
			IntersectionPoints.reserve(2);
			IntersectionVel.reserve(2);
			l=0.0;

		





			///////////////////////////////////////////////////////////////////////////////////////////////////////////
			//													//
			//			SEARCHING FOR THE LARGEST EDGE SIZE in either model part			//
			//													//
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			//search inside the first model part
			for(ModelPart::ElementsContainerType::iterator im = origin_model_part.ElementsBegin() ; 
				im != origin_model_part.ElementsEnd() ; ++im)
			{
			//we use the distance from the element center to the corner points
			x0=(im->GetGeometry()[0]).X();
			y0=(im->GetGeometry()[0]).Y();

			x1=(im->GetGeometry()[1]).X();
			y1=(im->GetGeometry()[1]).Y();

			x2=(im->GetGeometry()[2]).X();
			y2=(im->GetGeometry()[2]).Y();

			xc=0.333333333*(x0+x1+x2);
			yc=0.333333333*(y0+y1+y2);

			l0=sqrt((xc-x0)*(xc-x0)+(yc-y0)*(yc-y0));
			if (l0>l)
				l=l0;
			l1=sqrt((xc-x0)*(xc-x0)+(yc-y0)*(yc-y0));
			if (l1>l)
				l=l1;
			l2=sqrt((xc-x0)*(xc-x0)+(yc-y0)*(yc-y0));
			if (l2>l)
				l=l2;

			}
			//second model part
			for(ModelPart::ElementsContainerType::iterator im = destination_model_part.ElementsBegin() ; 
				im != destination_model_part.ElementsEnd() ; ++im)
			{
			//we use the distance from the element center to the corner points
			x0=im->GetGeometry()[0].X();
			y0=im->GetGeometry()[0].Y();

			x1=im->GetGeometry()[1].X();
			y1=im->GetGeometry()[1].Y();

			x2=im->GetGeometry()[2].X();
			y2=im->GetGeometry()[2].Y();

			xc=0.333333333*(x0+x1+x2);
			yc=0.333333333*(y0+y1+y2);

			l0=sqrt((xc-x0)*(xc-x0)+(yc-y0)*(yc-y0));
			if (l<l0)
				l=l0;
			l1=sqrt((xc-x0)*(xc-x0)+(yc-y0)*(yc-y0));
			if (l1>l)
				l=l1;
			l2=sqrt((xc-x0)*(xc-x0)+(yc-y0)*(yc-y0));
			if (l2>l)
				l=l2;

			}
			radius=l;
			KRATOS_WATCH("The search radius is");
			KRATOS_WATCH(radius);

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////IS_INTERFACE IDENTIFICATION -- can be optimized or changed later
		//save all the nodes of the conditions of interface	
				
		PointVector list_dest_cond_nodes;
		//list_dest_cond_nodes.reserve(AuxModelPart.Nodes().size());
		list_dest_cond_nodes.reserve(destination_model_part.Nodes().size());
		
		for(ModelPart::NodesContainerType::iterator in = destination_model_part.NodesBegin() ; 
				in != destination_model_part.NodesEnd() ; ++in)
			{
				//reset the IS_INTERFACE to zero for all the nodes of interface conditions
				in->FastGetSolutionStepValue(IS_INTERFACE)=0.0;				
				(list_dest_cond_nodes).push_back(*(in.base()));					
			}


		//create tree that contains all these nodes, i.e. the vertices of the interface elements
		kd_tree  nodes_tree2(list_dest_cond_nodes.begin(),list_dest_cond_nodes.end(), bucket_size);			
		
		//now we loop over the elements of origin and check if there are any nodes of destination, that lie inside
		//if the lie inside  - set some flag as TRUE (or some auxilliary var)
		//actually - only the elements that have at least one node of IS_BOUNDARY have to be checked
	
		for(ModelPart::ElementsContainerType::iterator im = origin_model_part.ElementsBegin() ; 
				im != origin_model_part.ElementsEnd() ; ++im)
		{
		//we use the distance from the element center to the corner points
		double xx0=im->GetGeometry()[0].X();
		double yy0=im->GetGeometry()[0].Y();

		double xx1=im->GetGeometry()[1].X();
		double yy1=im->GetGeometry()[1].Y();

		double xx2=im->GetGeometry()[2].X();
		double yy2=im->GetGeometry()[2].Y();

		double xxc=0.333333333*(xx0+xx1+xx2);
		double yyc=0.333333333*(yy0+yy1+yy2);

		work_point[0]=xxc;
		work_point[1]=yyc;
		work_point[2]=0.0;

		//number of IS_BOUNDARY
		int n_points_in_radius;

						
			n_points_in_radius = nodes_tree2.SearchInRadius(work_point, radius, res1.begin(),res_distances1.begin(), max_results);
			if (n_points_in_radius>=1.0)
				{
				
				//KRATOS_WATCH(n_points_in_radius)
				//KRATOS_WATCH(work_point)
				for(PointIterator i=res1.begin(); i!=res1.begin() + n_points_in_radius ; i++)
					{
					array_1d<double,3> N;
					bool is_inside=false;
					double dest_x = (*i)->X();
					double dest_y = (*i)->Y();
					//KRATOS_WATCH(dest_x)
					//KRATOS_WATCH(dest_y)
					is_inside=CalculatePosition( xx0, yy0, xx1, yy1, xx2, yy2, dest_x, dest_y, N);	
					//if destination node is inside the origin
					if (is_inside==true)
						{
						(*i)->FastGetSolutionStepValue(IS_INTERFACE)=1.0;						
						}
					}
						

				}
			

			
		}

		
		
		
			
		KRATOS_CATCH("")
		}

		//**********************************************************************************************
		//**********************************************************************************************
		void FindIntersectionOfEdges(ModelPart& origin_model_part, ModelPart& destination_model_part, ModelPart& AuxModelPart)
		{
			KRATOS_TRY  
			//at each step we need to make AuxModelPart empty
			AuxModelPart.Conditions().clear();
			AuxModelPart.Elements().clear();
			AuxModelPart.Nodes().clear();

			//FREE ALL THE AUX_VELS of the destination model part!
			for(ModelPart::NodesContainerType::iterator in = destination_model_part.NodesBegin() ; 
				in != destination_model_part.NodesEnd() ; ++in)
			{
				//reset the IS_INTERFACE to zero for all the nodes of interface conditions
				in->Free(AUX_VEL_X);
				in->Free(AUX_VEL_Y);
				in->FastGetSolutionStepValue(AUX_VEL_X)=0.0;				
				in->FastGetSolutionStepValue(AUX_VEL_Y)=0.0;

			}


			//first we need to find the biggest element size among all the elements of both the model parts
			// it will be set as the seacrh radius for the tree-search
			typedef Node<3> PointType;
			typedef Node<3>::Pointer PointPointerType;
			typedef std::vector<PointType::Pointer>           PointVector;
			typedef PointVector::iterator PointIterator;
			typedef std::vector<double>               DistanceVector;
			typedef std::vector<double>::iterator     DistanceIterator;


			

			// bucket types
			//typedef Bucket<3, PointType, ModelPart::NodesContainerType, PointPointerType, PointIterator, DistanceIterator > BucketType;
			//typedef Bins< 3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > StaticBins;
			// bucket types
			typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
		
				
			//*************
			// DynamicBins;	
			typedef Tree< KDTreePartition<BucketType> > kd_tree; //Kdtree;
			//typedef Tree< StaticBins > Bin; 			     //Binstree;
			unsigned int bucket_size = 50;
			
			//performing the interpolation - all of the nodes in this list will be preserved
			unsigned int max_results = 200;
			//PointerVector<PointType> res(max_results);
			//NodeIterator res(max_results);
			PointVector res(max_results);
			PointVector res1(max_results);
			DistanceVector res_distances(max_results);
			DistanceVector res_distances1(max_results);
			Node<3> work_point(0,0.0,0.0,0.0);
			array_1d<double, 3> temp;

			array_1d<double, 3> vel;
			array_1d<double, 3> vel1;
			array_1d<double, 3> vel2;


			double x0, y0, x1, y1, x2, y2, xc, yc, l0, l1, l2, l, radius, x0_orig, x1_orig, y0_orig, y1_orig;
			//PointVector IntersectionPoints;

			std::vector<array_1d<double,3> > IntersectionPoints;
			std::vector<array_1d<double,3> > IntersectionVel;
			//there should be just two intersection points per element
			IntersectionPoints.reserve(2);
			IntersectionVel.reserve(2);
			l=0.0;

		





			///////////////////////////////////////////////////////////////////////////////////////////////////////////
			//													//
			//			SEARCHING FOR THE LARGEST EDGE SIZE in either model part			//
			//													//
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			//search inside the first model part
			for(ModelPart::ElementsContainerType::iterator im = origin_model_part.ElementsBegin() ; 
				im != origin_model_part.ElementsEnd() ; ++im)
			{
			//we use the distance from the element center to the corner points
			x0=(im->GetGeometry()[0]).X();
			y0=(im->GetGeometry()[0]).Y();

			x1=(im->GetGeometry()[1]).X();
			y1=(im->GetGeometry()[1]).Y();

			x2=(im->GetGeometry()[2]).X();
			y2=(im->GetGeometry()[2]).Y();

			xc=0.333333333*(x0+x1+x2);
			yc=0.333333333*(y0+y1+y2);

			l0=sqrt((xc-x0)*(xc-x0)+(yc-y0)*(yc-y0));
			if (l0>l)
				l=l0;
			l1=sqrt((xc-x0)*(xc-x0)+(yc-y0)*(yc-y0));
			if (l1>l)
				l=l1;
			l2=sqrt((xc-x0)*(xc-x0)+(yc-y0)*(yc-y0));
			if (l2>l)
				l=l2;

			}
			//second model part
			for(ModelPart::ElementsContainerType::iterator im = destination_model_part.ElementsBegin() ; 
				im != destination_model_part.ElementsEnd() ; ++im)
			{
			//we use the distance from the element center to the corner points
			x0=im->GetGeometry()[0].X();
			y0=im->GetGeometry()[0].Y();

			x1=im->GetGeometry()[1].X();
			y1=im->GetGeometry()[1].Y();

			x2=im->GetGeometry()[2].X();
			y2=im->GetGeometry()[2].Y();

			xc=0.333333333*(x0+x1+x2);
			yc=0.333333333*(y0+y1+y2);

			l0=sqrt((xc-x0)*(xc-x0)+(yc-y0)*(yc-y0));
			if (l<l0)
				l=l0;
			l1=sqrt((xc-x0)*(xc-x0)+(yc-y0)*(yc-y0));
			if (l1>l)
				l=l1;
			l2=sqrt((xc-x0)*(xc-x0)+(yc-y0)*(yc-y0));
			if (l2>l)
				l=l2;

			}
			radius=l;
			KRATOS_WATCH("The search radius is");
			KRATOS_WATCH(radius);

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////IS_INTERFACE IDENTIFICATION -- can be optimized or changed later
		//save all the nodes of the conditions of interface	
			
		PointVector list_dest_cond_nodes;
		//list_dest_cond_nodes.reserve(AuxModelPart.Nodes().size());
		list_dest_cond_nodes.reserve(destination_model_part.Nodes().size());
		
		for(ModelPart::NodesContainerType::iterator in = destination_model_part.NodesBegin() ; 
				in != destination_model_part.NodesEnd() ; ++in)
			{
				//reset the IS_INTERFACE to zero for all the nodes of interface conditions
				in->FastGetSolutionStepValue(IS_INTERFACE)=0.0;				
				(list_dest_cond_nodes).push_back(*(in.base()));					
			}


		//create tree that contains all these nodes, i.e. the vertices of the interface elements
		kd_tree  nodes_tree2(list_dest_cond_nodes.begin(),list_dest_cond_nodes.end(), bucket_size);			
		
		//now we loop over the elements of origin and check if there are any nodes of destination, that lie inside
		//if the lie inside  - set some flag as TRUE (or some auxilliary var)
		//actually - only the elements that have at least one node of IS_BOUNDARY have to be checked
	
		for(ModelPart::ElementsContainerType::iterator im = origin_model_part.ElementsBegin() ; 
				im != origin_model_part.ElementsEnd() ; ++im)
		{
		//we use the distance from the element center to the corner points
		double xx0=im->GetGeometry()[0].X();
		double yy0=im->GetGeometry()[0].Y();

		double xx1=im->GetGeometry()[1].X();
		double yy1=im->GetGeometry()[1].Y();

		double xx2=im->GetGeometry()[2].X();
		double yy2=im->GetGeometry()[2].Y();

		double xxc=0.333333333*(xx0+xx1+xx2);
		double yyc=0.333333333*(yy0+yy1+yy2);

		work_point[0]=xxc;
		work_point[1]=yyc;
		work_point[2]=0.0;

		//number of IS_BOUNDARY
		int n_points_in_radius;

						
			n_points_in_radius = nodes_tree2.SearchInRadius(work_point, radius, res1.begin(),res_distances1.begin(), max_results);
			if (n_points_in_radius>=1.0)
				{
				
				//KRATOS_WATCH(n_points_in_radius)
				//KRATOS_WATCH(work_point)
				for(PointIterator i=res1.begin(); i!=res1.begin() + n_points_in_radius ; i++)
					{
					array_1d<double,3> N;
					bool is_inside=false;
					double dest_x = (*i)->X();
					double dest_y = (*i)->Y();
					//KRATOS_WATCH(dest_x)
					//KRATOS_WATCH(dest_y)
					is_inside=CalculatePosition( xx0, yy0, xx1, yy1, xx2, yy2, dest_x, dest_y, N);	
					//if destination node is inside the origin
					if (is_inside==true)
						{
						(*i)->FastGetSolutionStepValue(IS_INTERFACE)=1.0;						
						}
					}
						

				}
			

			
		}

		
		///////////////////////////////////////////////////////////////////////////////////////////////

			///////////////////////////////////////////////////////////////////////////////////////////////////////////
			//													//
			//			FINISHED SEARCHING FOR THE LARGEST EDGE SIZE					//
			//													//
			//////////////////////////////////////////////////////////////////////////////////////////////////////////


			//here we shall store the centers of the elements of Eulerian mesh, that belong to the interface with Lagrangian part
			PointVector list_of_nodes;
			list_of_nodes.reserve(origin_model_part.Nodes().size());
		
			//number of IS_INTERFACE nodes
			//unsigned int n_int=0;
			//number of IS_BOUNDARY nodes
			//unsigned int n_b=0;

			ModelPart:: ConditionsContainerType list_of_conditions;

			for(ModelPart::NodesContainerType::iterator in = origin_model_part.NodesBegin() ; 
				in != origin_model_part.NodesEnd() ; ++in)
			{
			//we are interested only in the nodes that belong to the boundary edges of the Origin Model Part
			// there mesher creates conditions. So we shall store these nodes
			if ((in->GetValue(NEIGHBOUR_CONDITIONS)).size()>0.0)
				{
				(list_of_nodes).push_back(*(in.base()));	
				}
			//(list_of_nodes).push_back(ic->GetGeometry()(1));
			
			}


			//create tree
			kd_tree  nodes_tree1(list_of_nodes.begin(),list_of_nodes.end(), bucket_size);
			// and now we start the search, only inside those elements of destiantion_model_part, that have two nodes of IS_INTERFACE
			for(ModelPart::ElementsContainerType::iterator im = destination_model_part.ElementsBegin() ; 
				im != destination_model_part.ElementsEnd() ; ++im)
			{
			
			IntersectionPoints.clear();
			IntersectionVel.clear();


			
				
				x0=im->GetGeometry()[0].X();
				y0=im->GetGeometry()[0].Y();

				x1=im->GetGeometry()[1].X();
				y1=im->GetGeometry()[1].Y();

				x2=im->GetGeometry()[2].X();
				y2=im->GetGeometry()[2].Y();

				xc=0.333333333*(x0+x1+x2);
				yc=0.333333333*(y0+y1+y2);
			

				work_point[0]=xc;
				work_point[1]=yc;
				work_point[2]=0.0;

				double n_points_in_radius = nodes_tree1.SearchInRadius(work_point, radius, res.begin(),res_distances.begin(), max_results);
				array_1d<double,3> N;
				bool is_inside=false;
				bool is_inside1=false;
				bool is_inside2=false;
				Matrix system(2,2);
				Matrix inv_sys(2,2);
				double det;
				Vector rhs(2);
				Vector solution(2);
				if (n_points_in_radius>=1)
				{										
					//find to which edge (condition) of the origin mesh this point belongs
				  for(PointIterator i=res.begin(); i!=res.begin() + static_cast<int>(n_points_in_radius) ; i++)
						{
						WeakPointerVector< Condition >& neighb_conds = (*i)->GetValue(NEIGHBOUR_CONDITIONS);
						
						for (unsigned int ii=0;ii<neighb_conds.size();ii++)
							{
							//here we shall find whether 
							x0_orig=neighb_conds[ii].GetGeometry()[0].X();
							y0_orig=neighb_conds[ii].GetGeometry()[0].Y();
							
							x1_orig=neighb_conds[ii].GetGeometry()[1].X();
							y1_orig=neighb_conds[ii].GetGeometry()[1].Y();							
	
							//we also store the nodal velocity of the conditions verices, in odrer to store finally
							//the velocity of the intersection points
							vel1=neighb_conds[ii].GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
							vel2=neighb_conds[ii].GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
							

							//check whether both points are inside the element, or outside or one is in, on out (last one is the case of interest)
							is_inside1=CalculatePosition( x0, y0, x1, y1, x2, y2, x0_orig, y0_orig, N);							is_inside2=CalculatePosition( x0, y0, x1, y1, x2, y2, x1_orig, y1_orig,
N);				
							//if one is inside, and another - outside, the origin condition is intersecting destination edge
							//we have to find the intersection point.. must check for all three edges of the destination
							if ((is_inside1+is_inside2)==1.0)
								{
								is_inside=false;
								//first edge check
								system(0,0)=y0_orig-y1_orig;
								system(1,0)=y0-y1;
								system(0,1)=x1_orig-x0_orig;
								system(1,1)=x1-x0;

								rhs[0]=y0_orig*(x1_orig-x0_orig)-x0_orig*(y1_orig-y0_orig);
								rhs[1]=y0*(x1-x0)-x0*(y1-y0);
								MathUtils<double>::InvertMatrix(system, inv_sys, det);
								
								if (fabs(det)>0.00000000001)
								{							
								solution=prod(inv_sys, rhs);
								
								//and now make sure that the intersection is lying inside the destination element
								is_inside=IsOnEdge( x0, y0, x1, y1, x2, y2, solution[0], solution[1], x0_orig,  x1_orig, y0_orig, y1_orig,N);
								if (is_inside==true)
									{
									temp[0]=solution[0];
									temp[1]=solution[1];
									temp[2]=0.0;
									//push back the solution, only if it is not there yet
									if (!IsAlreadyInList(temp,IntersectionPoints))
										{
										IntersectionPoints.push_back(temp);
										//velocity of the intersection point, knowing vel1 and vel2 of the vertices of the condition: vel=N1(intersection)*v1+N2(intersection)*v2

										vel=vel2*sqrt((solution[0]-x0_orig)*(solution[0]-x0_orig)+(solution[1]-y0_orig)*(solution[1]-y0_orig))/sqrt((x1_orig-x0_orig)*(x1_orig-x0_orig)+(y1_orig-y0_orig)*(y1_orig-y0_orig));
										
										vel+=vel1*sqrt((solution[0]-x1_orig)*(solution[0]-x1_orig)+(solution[1]-y1_orig)*(solution[1]-y1_orig))/sqrt((x1_orig-x0_orig)*(x1_orig-x0_orig)+(y1_orig-y0_orig)*(y1_orig-y0_orig));

										IntersectionVel.push_back(vel);
										//KRATOS_WATCH("Found intersection 1, coords are")
										//KRATOS_WATCH(solution)
										//KRATOS_WATCH("Vel at intersection is")
										//KRATOS_WATCH(vel)
										
										}
																	
									}
								}
								//HERE WE SHOULD PROBABLY IMPOSE THE DIRICHLET!
								is_inside=false;
								//second edge check
								system(0,0)=y0_orig-y1_orig;
								system(1,0)=y0-y2;
								system(0,1)=x1_orig-x0_orig;
								system(1,1)=x2-x0;

								rhs[0]=y0_orig*(x1_orig-x0_orig)-x0_orig*(y1_orig-y0_orig);
								rhs[1]=y0*(x2-x0)-x0*(y2-y0);
								MathUtils<double>::InvertMatrix(system, inv_sys, det);
								
								if (fabs(det)>0.00000000001)	
								{							
								solution=prod(inv_sys, rhs);
								
								//and now make sure that the intersection is lying inside the destination element
								is_inside=IsOnEdge( x0, y0, x1, y1, x2, y2, solution[0], solution[1], x0_orig, x1_orig, y0_orig, y1_orig,N);

								if (is_inside==true)
									{									
									temp[0]=solution[0];
									temp[1]=solution[1];
									temp[2]=0.0;
									//push back the solution, only if it is not there yet
									if (!IsAlreadyInList(temp,IntersectionPoints))
										{
										IntersectionPoints.push_back(temp);
//velocity of the intersection point, knowing vel1 and vel2 of the vertices of the condition: vel=N1(intersection)*v1+N2(intersection)*v2
										vel=vel2*sqrt((solution[0]-x0_orig)*(solution[0]-x0_orig)+(solution[1]-y0_orig)*(solution[1]-y0_orig))/sqrt((x1_orig-x0_orig)*(x1_orig-x0_orig)+(y1_orig-y0_orig)*(y1_orig-y0_orig));
										
										vel+=vel1*sqrt((solution[0]-x1_orig)*(solution[0]-x1_orig)+(solution[1]-y1_orig)*(solution[1]-y1_orig))/sqrt((x1_orig-x0_orig)*(x1_orig-x0_orig)+(y1_orig-y0_orig)*(y1_orig-y0_orig));

										IntersectionVel.push_back(vel);
										//KRATOS_WATCH("Found intersection 2, coords are")
										//KRATOS_WATCH(solution)
										//KRATOS_WATCH("Vel at intersection is")
										//KRATOS_WATCH(vel)
										
										}
									}
								}
								//HERE WE SHOULD PROBABLY IMPOSE THE DIRICHLET!
								is_inside=false;
								//third edge check
								system(0,0)=y0_orig-y1_orig;
								system(1,0)=y1-y2;
								system(0,1)=x1_orig-x0_orig;
								system(1,1)=x2-x1;

								rhs[0]=y0_orig*(x1_orig-x0_orig)-x0_orig*(y1_orig-y0_orig);
								rhs[1]=y1*(x2-x1)-x1*(y2-y1);
								MathUtils<double>::InvertMatrix(system, inv_sys, det);
								if (fabs(det)>0.00000000001)
								{							
								solution=prod(inv_sys, rhs);
								
								//and now make sure that the intersection is lying inside the destination element
								is_inside=IsOnEdge( x0, y0, x1, y1, x2, y2, solution[0], solution[1], x0_orig, x1_orig, y0_orig, y1_orig, N);
  
								if (is_inside==true)
									{									
									temp[0]=solution[0];
									temp[1]=solution[1];
									temp[2]=0.0;
									//push back the solution, only if it is not there yet
									if (!IsAlreadyInList(temp,IntersectionPoints))
										{
										IntersectionPoints.push_back(temp);
//velocity of the intersection point, knowing vel1 and vel2 of the vertices of the condition: vel=N1(intersection)*v1+N2(intersection)*v2
										vel=vel2*sqrt((solution[0]-x0_orig)*(solution[0]-x0_orig)+(solution[1]-y0_orig)*(solution[1]-y0_orig))/sqrt((x1_orig-x0_orig)*(x1_orig-x0_orig)+(y1_orig-y0_orig)*(y1_orig-y0_orig));
										
										vel+=vel1*sqrt((solution[0]-x1_orig)*(solution[0]-x1_orig)+(solution[1]-y1_orig)*(solution[1]-y1_orig))/sqrt((x1_orig-x0_orig)*(x1_orig-x0_orig)+(y1_orig-y0_orig)*(y1_orig-y0_orig));

										IntersectionVel.push_back(vel);
										//KRATOS_WATCH("Found intersection 3, coords are")

										//KRATOS_WATCH(solution)
										//KRATOS_WATCH("Vel at intersection is")
										//KRATOS_WATCH(vel)
										
										}
									}
								}

								//HERE WE SHOULD PROBABLY IMPOSE THE DIRICHLET!
								}

								
							}

						}

					}
					
			
				
			//KRATOS_WATCH(im->GetId())
			//KRATOS_WATCH(IntersectionPoints[0])
			//if (IntersectionPoints.size()!=0)
			//		{
					//KRATOS_WATCH(IntersectionPoints.size())
					//KRATOS_WATCH(IntersectionVel.size())
			//		}
			if (IntersectionPoints.size()==2)
			//if there are two intersection points 
			{
				int id=AuxModelPart.Conditions().size()+1;
				Properties::Pointer properties = destination_model_part.GetMesh().pGetProperties(1);	
				Geometry< Node<3> >::Pointer geom = im->pGetGeometry();
				//checks if the intersection is too close to any of the vertices
				bool bad_intersection=false;
				//create new condition
				
				array_1d<double,3> Point1;
				array_1d<double,3> Point2;
				Point1=IntersectionPoints[0];
				Point2=IntersectionPoints[1];
				//push it back to the proj_dirichlet_model_part (AuxModelPart we call it))
				//check that both intersections do not belong to the same edge
				int which_edge1=100;
				int which_edge2=100;
				CalculateN_at_Point(im->GetGeometry(), Point1[0], Point1[1], N);
				for (int i=0;i<3;i++)
				{
				if (N[i]<0.00000000000001)
					which_edge1=i;
				}
				if ( (N[0]<0.2 && N[1]<0.2) || (N[0]<0.2 && N[2]<0.2) || (N[2]<0.2 && N[1]<0.2))
					bad_intersection=true;

				CalculateN_at_Point(im->GetGeometry(), Point2[0], Point2[1], N);
				for (int i=0;i<3;i++)
				{
				if (N[i]<0.00000000000001)
					which_edge2=i;
				}
				if ( (N[0]<0.2 && N[1]<0.2) || (N[0]<0.2 && N[2]<0.2) || (N[2]<0.2 && N[1]<0.2))
					bad_intersection=true;

				//add only if two intersections intersect two different edges, and no intersection is close to the vertex
				
				if (which_edge1!=which_edge2 && bad_intersection==false) 
					{
					Condition::Pointer p_condition(new ProjDirichletCond(id, geom,properties, IntersectionPoints[0], IntersectionPoints[1], IntersectionVel[0], IntersectionVel[1]));
					AuxModelPart.Conditions().push_back(p_condition);
					}
			
				/*
				aaaa=IntersectionPoints[0];
				if (aaaa[0]==0.5 && aaaa[1]==0.4)
					AuxModelPart.Conditions().push_back(p_condition);

				bbbb=IntersectionPoints[1];
				if (bbbb[0]==0.5 && bbbb[1]==0.4)
					AuxModelPart.Conditions().push_back(p_condition);

				if (bbbb[0]==0.5 && bbbb[1]==0.5)
					AuxModelPart.Conditions().push_back(p_condition);
			
				*/
			}
			
						
			}
		
		
		
			
		KRATOS_CATCH("")
		}
		///////////////////////////////////////////////////////////////////////////
		////////	SUBDOMAIN DISABLING			//////////////////
		///////////////////////////////////////////////////////////////////////////
		void DisableSubdomain(ModelPart& full_model_part, ModelPart& reduced_model_part)
		{
		KRATOS_TRY
		int n_int;
		int n_disabled;
		//clear reduced_model_part
		reduced_model_part.Conditions().clear();
		reduced_model_part.Elements().clear();
		reduced_model_part.Nodes().clear();

		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
			
			if (n_int==3)
				{
				im->GetGeometry()[0].FastGetSolutionStepValue(DISABLE)=true;
				im->GetGeometry()[1].FastGetSolutionStepValue(DISABLE)=true;
				im->GetGeometry()[2].FastGetSolutionStepValue(DISABLE)=true;
				}

		}
		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			n_disabled=im->GetGeometry()[0].FastGetSolutionStepValue(DISABLE);
			n_disabled+=im->GetGeometry()[1].FastGetSolutionStepValue(DISABLE);
			n_disabled+=im->GetGeometry()[2].FastGetSolutionStepValue(DISABLE);
			
			if (n_disabled<3)
				{
				reduced_model_part.AddElement(*(im.base()));
				//reduced_model_part.AddNode(im->GetGeometry()[0]);				
				}
			if (n_disabled>3)
				KRATOS_ERROR(std::logic_error,  "Number of DISABLE flags cant exceed number of the element nodes.... " , "");

		}
	
		for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; 
				in != full_model_part.NodesEnd() ; ++in)
		{	  
			
			n_disabled=in->FastGetSolutionStepValue(DISABLE);
			
			if (n_disabled==0.0)
				{
				reduced_model_part.AddNode(*(in.base()));
				//reduced_model_part.AddNode(im->GetGeometry()[0]);				
				}
			if (n_disabled>3)
				KRATOS_ERROR(std::logic_error,  "Number of DISABLE flags cant exceed number of the element nodes.... " , "");

		}

		/*
		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			n_disabled=im->GetGeometry()[0].FastGetSolutionStepValue(DISABLE);
			n_disabled+=im->GetGeometry()[1].FastGetSolutionStepValue(DISABLE);
			n_disabled+=im->GetGeometry()[2].FastGetSolutionStepValue(DISABLE);
			
			if (n_disabled==1 || n_disabled==2)
				{
				for (int i=0;i<3;i++)
					{

					if (im->GetGeometry()[i].FastGetSolutionStepValue(DISABLE)==1)
						{
						reduced_model_part.AddNode(im->GetGeometry()(i));				

						}
	
					}
				}

		}
		*/
		/*
		for(ModelPart::NodesContainerType::iterator in = reduced_model_part.NodesBegin() ; 
				in != reduced_model_part.NodesEnd() ; ++in)
		{
		if (in->Y()>0.00001  && in->Y()<0.09999)
			{
			if (in->FastGetSolutionStepValue(DISABLE)!=1.0)
				{
				in->Free(VELOCITY_X);
				in->Free(VELOCITY_Y);
				}
			}
		}
		*/
		/*
		for(ModelPart::NodesContainerType::iterator in = reduced_model_part.NodesBegin() ; 
				in != reduced_model_part.NodesEnd() ; ++in)
		{	
			if (in->FastGetSolutionStepValue(DISABLE,1)==1 && in->FastGetSolutionStepValue(DISABLE)!=1)
				{
				in->Free(VELOCITY_X);
				in->Free(VELOCITY_Y);
				}

		}
		*/
				
	
		
		KRATOS_CATCH("")
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//													   //
		//			AUXILIARY FUNCTIONS								   //
		//													   //
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		inline void CalculateN_at_Point(Element::GeometryType& geom, const double xc, const double yc, array_1d<double,3>& N_at_c)
	{
			//first we calculate the area of the whole triangle			
			double x10 = geom[1].X() - geom[0].X();
			double y10 = geom[1].Y() - geom[0].Y();
			
			double x20 = geom[2].X() - geom[0].X();
			double y20 = geom[2].Y() - geom[0].Y();
			
			double detJ = x10 * y20-y10 * x20;
			double totArea=0.5*detJ;			
			//and now we calculate the areas of three respective triangle, that (xc,yc) divide the original one into
			// xc, 0, 1
			double x0c = geom[0].X() - xc ;
			double y0c = geom[0].Y() - yc ;
			
			double x1c = geom[1].X() - xc;
			double y1c = geom[1].Y() - yc;

			double x2c = geom[2].X() - xc;
			double y2c = geom[2].Y() - yc;
			//xc, 0, 1
			detJ= x0c * y1c - y0c * x1c;
			double Area2 = 0.5*detJ;
			//xc, 0, 2
			detJ= x0c * y2c - y0c * x2c;
			double Area1 = 0.5*detJ;
			//xc, 1, 2
			detJ= x1c * y2c - y1c * x2c;
			double Area0 = 0.5*detJ;

			if (totArea<0.00000000000000001)
				KRATOS_ERROR(std::logic_error,  "Your element Proj DIrichlet Cond has a zero area!!!! " , "");
			//and now we fill in the array of shape functions values:
			// 1 0 2
			N_at_c[0]=fabs(Area0/totArea);
			N_at_c[1]=fabs(Area1/totArea);
			N_at_c[2]=fabs(Area2/totArea);
			if (  (N_at_c[0]<0.05 && N_at_c[1]<0.05) || (N_at_c[0]<0.05 && N_at_c[2]<0.05) || (N_at_c[2]<0.05 && N_at_c[1]<0.05))			
			KRATOS_WATCH("Dangerous VERTICES!!!")		
			//KRATOS_ERROR(std::logic_error,  "Too close to the node is the INTERSECTION!!!! " , "")

	}
		//////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
inline double CalculateVol(	const double x0, const double y0,
						const double x1, const double y1,
    						const double x2, const double y2
					  )
		{
			return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
		}
		//////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		inline bool IsAlreadyInList(array_1d<double,3>& current_point, std::vector<array_1d<double,3> >& IntersectionPointsList)
		{
		for (unsigned int i=0;i<IntersectionPointsList.size();i++)
			{
			//temp=IntersectionPointsList[i];
			if (std::equal(current_point.begin(), current_point.end(), IntersectionPointsList[i].begin()))
				{				
				return true;
				}
			}
		return false;
	
		}
		//////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		inline bool CalculatePosition(	const double x0, const double y0,
						const double x1, const double y1,
   						const double x2, const double y2,
						const double xc, const double yc,
						array_1d<double,3>& N		
					  )
		{
			double area = CalculateVol(x0,y0,x1,y1,x2,y2);
			double inv_area = 0.0;
			if(area < 0.000000000001)
			  {
				KRATOS_ERROR(std::logic_error,"element with zero area found","");
			  }
			else
			  {
				inv_area = 1.0 / area;
			  }
			
			  
			  N[0] = CalculateVol(x1,y1,x2,y2,xc,yc) * inv_area;
			  N[1] = CalculateVol(x2,y2,x0,y0,xc,yc) * inv_area;
			  N[2] = CalculateVol(x0,y0,x1,y1,xc,yc) * inv_area;
			  
/*			  N[0] = CalculateVol(x0,y0,x1,y1,xc,yc) * inv_area;
			N[1] = CalculateVol(x1,y1,x2,y2,xc,yc) * inv_area;
			N[2] = CalculateVol(x2,y2,x0,y0,xc,yc) * inv_area;*/
			
			if(N[0] > 0.0 && N[1] > 0.0 && N[2] > 0.0 && N[0] < 1.0 && N[1] < 1.0 && N[2] < 1.0) //if the xc yc is inside the triangle
			//if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
				return true;
			
			return false;
		}	
		//////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		inline bool IsOnEdge(	const double x0, const double y0,
						const double x1, const double y1,
   						const double x2, const double y2,
						const double sol_x, const double sol_y,
						const double x0_orig, const double x1_orig,
						const double y0_orig, const double y1_orig,
						array_1d<double,3>& N		
					  )
		{
			double area = CalculateVol(x0,y0,x1,y1,x2,y2);
			double inv_area = 0.0;
			if(area < 0.000000000001)
			  {
				KRATOS_ERROR(std::logic_error,"element with zero area found","");
			  }
			else
			  {
				inv_area = 1.0 / area;
			  }
			
			  
			  N[0] = CalculateVol(x1,y1,x2,y2,sol_x,sol_y) * inv_area;
			  N[1] = CalculateVol(x2,y2,x0,y0,sol_x,sol_y) * inv_area;
			  N[2] = CalculateVol(x0,y0,x1,y1,sol_x,sol_y) * inv_area;
			  
			//the intersection of the lines should be within the element and also inside the line segment (origin condition) defined by x0_orig, x1_orig	
			//if ((N[0]==0.0 || N[1]==0.0 || N[2]==0.0) && ( (x1_orig-sol_x)*(sol_x-x0_orig)>0.0 || (y1_orig-sol_y)*(sol_y-y0_orig)>0.0) )
			//double t1=(x1_orig-sol_x)*(sol_x-x0_orig);		
			//double t2=(y1_orig-sol_y)*(sol_y-y0_orig);
			if( (N[0] >= -0.0000000000001 && N[1] >= -0.00000000000001 && N[2] >= -0.0000000000001 && N[0] <= 1.00000000000001 && N[1] <= 1.00000000000001 && N[2] <= 1.00000000000001) && ((x1_orig-sol_x)*(sol_x-x0_orig)>0.0 || (y1_orig-sol_y)*(sol_y-y0_orig)>0.0) )
				return true;
						
			return false;
		}	

	private:

		//aux vars
		static boost::numeric::ublas::bounded_matrix<double,3,3> msJ; //local jacobian
		static boost::numeric::ublas::bounded_matrix<double,3,3> msJinv; //inverse jacobian
		static array_1d<double,3> msc; //center pos
		static array_1d<double,3> ms_rhs; //center pos
		

	};

}  // namespace Kratos.

#endif // KRATOS_COUPLED_EUL_ULF_UTILITIES_INCLUDED  defined 


