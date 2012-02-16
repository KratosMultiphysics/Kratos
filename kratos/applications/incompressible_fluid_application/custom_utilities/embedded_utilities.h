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
//   Last Modified by:    $Author: pavel $
//   Date:                $Date: 2009-01-15 14:50:24 $
//   Revision:            $Revision: 1.12 $
//
//

#if !defined(KRATOS_EMBEDDED_UTILITIES_INCLUDED )
#define  KRATOS_EMBEDDED_UTILITIES_INCLUDED



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
	class EmbeddedUtils
	{
	public:
		typedef Node<3> NodeType;
		//**********************************************************************************************
		//**********************************************************************************************
		// function that applies ADDITIONAL DIRICHLET CONDITIONS UPON THE nodes that lie at the INTERFACE between Lagrangain and Eulerian parts
		/*void ApplyProjDirichlet(ModelPart& full_model_part, ModelPart& aux_conditions_model_part)
		{
		KRATOS_TRY
		
		KRATOS_ERROR(std::logic_error,  "Use the ApplyProjDirichlet process instead .... " , "");
		
		KRATOS_CATCH("")
		}
		*/
		//this function stores the volume elements intersected by the embedded "skin" into a
		//separate model part
		//these elements are distinguished by containing nodes with positive and negative "distances" from the intersection
		
		/* THIS FUNCTION IS NOT WORKING AND IS NOT NECESSARY.. left it commented just in case-....
		void SaveInterfaceElemsModelPart(ModelPart& reduced_model_part, ModelPart& model_part)
		{
		model_part.Conditions().clear();
		model_part.Elements().clear();
		model_part.Nodes().clear();
		double distance=0.0;
		double distance_aux=0.0;
		bool intersected_elem=false;
		  for(ModelPart::ElementsContainerType::iterator im = reduced_model_part.ElementsBegin() ; 
				  im != reduced_model_part.ElementsEnd() ; ++im)
		    {
		     for (std::size_t i=0;i<im->GetGeometry().size();i++)
			  {
			  distance=distance_aux;
			  distance_aux=im->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
			  if (distance*distance_aux<0.0)
			      {
			      intersected_elem=true;
			      }
		     	  }
		     if (intersected_elem==true)
		       model_part.AddElement(*(im.base()));		  
		    }
		    for(ModelPart::NodesContainerType::iterator in = reduced_model_part.NodesBegin() ; 
				in != reduced_model_part.NodesEnd() ; ++in)
		      {
			distance=in->FastGetSolutionStepValue(DISTANCE);
			bool node_of_intersected_elem=false;
			WeakPointerVector< Node<3> >& neighb_nodes = in->GetValue(NEIGHBOUR_NODES); 
			for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++) 
				{ 					
				 distance_aux=i->FastGetSolutionStepValue(DISTANCE);
				 if (distance*distance_aux<0.0)
				   node_of_intersected_elem=true;
				}
				if (node_of_intersected_elem==true)
				  model_part.AddNode(*(in.base()));
				  
				
		      }
		KRATOS_WATCH(reduced_model_part)
		KRATOS_WATCH(model_part)
		  
		}
		*/
		////////////////////////////////////////////////////////////////////////////////////////
		void CreateIntersConditions(ModelPart& model_part, ModelPart& interface_conditions_model_part )
		{
		KRATOS_TRY

		interface_conditions_model_part.Conditions().clear();
		interface_conditions_model_part.Elements().clear();
		interface_conditions_model_part.Nodes().clear();
		//reset the IS_INTERFACE flag - if the distance is negative (i.e. the fluid node lies inside of the solid object and is fictitious - set 1 otherwise 0
		for(ModelPart::NodesContainerType::iterator in = model_part.NodesBegin() ;  in != model_part.NodesEnd() ; ++in)
		{
		if (in->FastGetSolutionStepValue(DISTANCE)<0.0)
			in->FastGetSolutionStepValue(IS_INTERFACE)=1.0;
		else
			in->FastGetSolutionStepValue(IS_INTERFACE)=0.0;
		}

		for(ModelPart::ElementsContainerType::iterator im = model_part.ElementsBegin() ;  im != model_part.ElementsEnd() ; ++im)
		{
		//intersection Points - at most we shall consider 4 points
		std::vector<array_1d<double,3> > IntersectionPoints;
		std::vector<array_1d<double,3> > IntersectionVel;

		IntersectionPoints.reserve(4);
		IntersectionVel.reserve(4);

		array_1d<double,3> Point;
		//identify 		
		int intersection_count=0;
		//iterating over edges (01 02 03 12 13 23)
		for (int i=0;i<3;i++)
			{
			for (int j=i+1;j<4;j++)
				{
				//std::cout<<"edge ij "<<i<<j<<std::endl;	
				double d0=im->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
				double d1=im->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE);	

				//if the product of distances of two nodes is negative - the edge is crossed
				if (d0*d1<0.0)
					{
					std::cout<<"Intersected edge "<<i<<j<<std::endl;	

					double x0=im->GetGeometry()[i].X();
					double x1=im->GetGeometry()[j].X();
					double k=0.0;
					double b=0.0;
					double x_inters=0.0;
					double y_inters=0.0;
					double z_inters=0.0;
					if (x1!=x0)
					  {
					  k=(d1-d0)/(x1-x0);
					  b=d0-k*x0;
					  x_inters=-b/k;
					  }
					else 
					  x_inters=x0;

					double y0=im->GetGeometry()[i].Y();
					double y1=im->GetGeometry()[j].Y();
					if (y1!=y0)
					  {
					  k=(d1-d0)/(y1-y0);
					  b=d0-k*y0;
					  y_inters=-b/k;
					  }
					else
					  y_inters=y0;

					double z0=im->GetGeometry()[i].Z();
					double z1=im->GetGeometry()[j].Z();
					if (z1!=z0)
					  {
					  k=(d1-d0)/(z1-z0);
					  b=d0-k*z0;
 					  z_inters=-b/k;
					  }
					else
					  z_inters=z0;

					Point[0]=x_inters;
					Point[1]=y_inters;
					Point[2]=z_inters;
		
					//KRATOS_WATCH(Point)

					IntersectionPoints.push_back(Point);
					intersection_count++;
					}

				}
			}
		
		//if the element is intersected by the embedded skin, create condition
		if (intersection_count!=0)
			{
			//the size of array should represent actual intersection number
			IntersectionPoints.resize(intersection_count);
			//KRATOS_WATCH(IntersectionPoints.size())
					
		
			//for now we assume zero velocity of the structure TO BE completed later...
			array_1d<double,3> ZeroVel=ZeroVector(3);
			for (unsigned int i=0;i<IntersectionVel.size();i++)
				IntersectionVel.push_back(ZeroVel);

			//////////////////////////////////////////////////////////////////////
			Geometry< Node<3> >::Pointer geom = im->pGetGeometry();
			Properties::Pointer properties = model_part.GetMesh().pGetProperties(1);	
			
			int id=interface_conditions_model_part.Conditions().size()+1;

			if (IntersectionPoints.size()==3)
				{
				Condition::Pointer p_condition(new ProjDirichletCond3D(id, geom,properties, IntersectionPoints[0], IntersectionPoints[1], IntersectionPoints[2], IntersectionVel[0], IntersectionVel[1], IntersectionVel[2]));
				interface_conditions_model_part.Conditions().push_back(p_condition);
				}
			else if (IntersectionPoints.size()==4)
				{
				Condition::Pointer p_condition(new ProjDirichletCond3D(id, geom,properties, IntersectionPoints[0], IntersectionPoints[1], IntersectionPoints[2], IntersectionPoints[3], IntersectionVel[0], IntersectionVel[1], 			IntersectionVel[2], IntersectionVel[3]));
			
				interface_conditions_model_part.Conditions().push_back(p_condition);
				}
			else 		
				KRATOS_ERROR(std::logic_error,  "Strange number of intersections - neither 3 nor 4 - check  CreateIntersectionConditions function" , "");
			
			}
		//end loop over elements
		}
		KRATOS_WATCH(interface_conditions_model_part)
		KRATOS_CATCH("")


		}
	
		
		///////////////////////////////////////////////////////////////////////////
		////////	SUBDOMAIN DISABLING			//////////////////
		///////////////////////////////////////////////////////////////////////////
		void DisableSubdomain(ModelPart& full_model_part, ModelPart& reduced_model_part)
		{
		KRATOS_TRY
		/*
		std::size_t n_int=0;
		std::size_t n_disabled=0;
		//clear reduced_model_part
		reduced_model_part.Conditions().clear();
		reduced_model_part.Elements().clear();
		reduced_model_part.Nodes().clear();
		// reset the DISABLED flag
		for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; 
				in != full_model_part.NodesEnd() ; ++in)
		{	  
			in->FastGetSolutionStepValue(DISABLE)=false;
		}
		//internal elements will be always disabled
		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			n_int=0;
			for (std::size_t i=0;i<im->GetGeometry().size();i++)
				n_int+=im->GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE);				
				

			if (n_int==im->GetGeometry().size())
				{
					
				  for (std::size_t i=0;i<im->GetGeometry().size();i++)
				      im->GetGeometry()[i].FastGetSolutionStepValue(DISABLE)=true;				
				}

		}
		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
	            n_disabled=0;

		    for (std::size_t i=0;i<im->GetGeometry().size();i++)
			n_disabled+=im->GetGeometry()[i].FastGetSolutionStepValue(DISABLE);
			

		if (n_disabled>0 && n_disabled<im->GetGeometry().size())
			{
			reduced_model_part.AddElement(*(im.base()));
			//reduced_model_part.AddNode(im->GetGeometry()[0]);				
			}
		else if (n_disabled>im->GetGeometry().size())
			{
			KRATOS_ERROR(std::logic_error,  "Number of DISABLE flags cant exceed number of the element nodes... " , "");
			}
			

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
			//and also add the fictitious nodes of the interface elements, i.e. the ones that have neighbor nodes, NOT DISABLED
			else 
				{
				WeakPointerVector< Node<3> >& neighb_nodes = in->GetValue(NEIGHBOUR_NODES); 
					unsigned int count=0;
					for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++) 
					{ 					
							count+=	i->FastGetSolutionStepValue(DISABLE);						
					}
				if (count<neighb_nodes.size()) //i.e. if not all the neighbor nodes are disabled, add the node 
					{
					reduced_model_part.AddNode(*(in.base()));
					//KRATOS_WATCH("ADDING THE NODE FOR WEAK IMPOSITION OF INTERFACE BOUNDARY CONDITION");
					}
				//set the velocity at the "interior" node (the one that is completely fictitious) to zero
				
				}
			*/

		int n_int;
		int n_disabled;
		//clear reduced_model_part
		reduced_model_part.Conditions().clear();
		reduced_model_part.Elements().clear();
		reduced_model_part.Nodes().clear();

		reduced_model_part.Conditions().reserve(full_model_part.Conditions().size());
		reduced_model_part.Elements().reserve(full_model_part.Elements().size());
		reduced_model_part.Nodes().reserve(full_model_part.Nodes().size());

		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			
			n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[3].FastGetSolutionStepValue(IS_INTERFACE);
			
			if (n_int==4)
				{
				im->GetGeometry()[0].FastGetSolutionStepValue(DISABLE)=true;
				im->GetGeometry()[1].FastGetSolutionStepValue(DISABLE)=true;
				im->GetGeometry()[2].FastGetSolutionStepValue(DISABLE)=true;
				im->GetGeometry()[3].FastGetSolutionStepValue(DISABLE)=true;
				}

		}
		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			n_disabled=im->GetGeometry()[0].FastGetSolutionStepValue(DISABLE);
			n_disabled+=im->GetGeometry()[1].FastGetSolutionStepValue(DISABLE);
			n_disabled+=im->GetGeometry()[2].FastGetSolutionStepValue(DISABLE);
			n_disabled+=im->GetGeometry()[3].FastGetSolutionStepValue(DISABLE);
			
			if (n_disabled<4)
				{
				reduced_model_part.Elements().push_back(*(im.base()));
				//reduced_model_part.AddNode(im->GetGeometry()[0]);				
				}
			
			//these are the nodes of intersected elements that lie on the fictitious part and are used for applying the Dirichlet boundary conditions
			if (n_disabled>0 && n_disabled<4)
				{
				for (int i=0;i<4;i++)
					{

					if (im->GetGeometry()[i].FastGetSolutionStepValue(DISABLE)==1)
						{
						reduced_model_part.Nodes().push_back(im->GetGeometry()(i));	
						}	
					}
				}

		}
		//and now we add all the nodes of "real" elements
		for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; 
				in != full_model_part.NodesEnd() ; ++in)
		{	  
			
			n_disabled=in->FastGetSolutionStepValue(DISABLE);
			
			if (n_disabled==0.0)
				{
				reduced_model_part.Nodes().push_back(*(in.base()));
				}
			

		}

		reduced_model_part.Nodes().Unique();

		for(ModelPart::PropertiesContainerType::iterator i_properties = full_model_part.PropertiesBegin() ; 
				i_properties != full_model_part.PropertiesEnd() ; ++i_properties)
		{	  
			reduced_model_part.AddProperties(*(i_properties.base()));			
		}

		for(ModelPart::ConditionsContainerType::iterator i_condition = full_model_part.ConditionsBegin() ; 
				i_condition != full_model_part.ConditionsEnd() ; ++i_condition)
		{	  
			reduced_model_part.AddCondition(*(i_condition.base()));			
		}		

		
				
		KRATOS_CATCH("")
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////
		void ApplyProjDirichlet(ModelPart& full_model_part)
		{
		KRATOS_TRY
		unsigned int n_old_int;
					
		//first we remove the Dirichlet conditions from the nodes that were defining the interface in the previous step:
		for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; 
				in != full_model_part.NodesEnd() ; ++in)
		{	
		n_old_int=in->FastGetSolutionStepValue(IS_INTERFACE,1);

		//make the velocity free at the nodes that are not the "Dangerosu ones"
		if (n_old_int>0.0 && in->FastGetSolutionStepValue(IS_INTERFACE)!=100.0)
			{
			//KRATOS_WATCH("OLD INTERFACEEEEEEEEEEEEEEEE!!!!!!!!!!!!!!!!!!!!!!1")
			in->FastGetSolutionStepValue(DISABLE)=false;
			in->Free(VELOCITY_X);
			in->Free(VELOCITY_Y);
			in->Free(VELOCITY_Z);

			in->Free(AUX_VEL_X);
			in->Free(AUX_VEL_Y);
			in->Free(AUX_VEL_Z);

			in->Free(PRESSURE);
			}
		//fix the IS_INt to 1 additionally at the nodes that are too close to the interface
		if (in->FastGetSolutionStepValue(IS_INTERFACE)==100.0)
				{
				in->FastGetSolutionStepValue(IS_INTERFACE)=1.0;
				KRATOS_WATCH("BAD NODE IS")
				KRATOS_WATCH(in->GetId())
				}		
				
		}

		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
					im != full_model_part.ElementsEnd() ; ++im)
		{
				double n_int=0.0;
				//iterate over the of nodes in the element (interface element)
				for (unsigned int i=0;i<im->GetGeometry().size();i++)
					{
					n_int=im->GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE);
					// if the node is lying on fictitious side - apply the project Dirichlet condition
					if (n_int==1.0)
					//if (((ic->GetGeometry()[i]).GetDof(AUX_VEL_X)).IsFixed())
						{						
						im->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)=im->GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL);
						im->GetGeometry()[i].Fix(VELOCITY_X);
						im->GetGeometry()[i].Fix(VELOCITY_Y);
						im->GetGeometry()[i].Fix(VELOCITY_Z);
						}
					}
	
		}
			
		
		
		
		KRATOS_CATCH("")
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//													   //
		//			AUXILIARY FUNCTIONS								   //
		//													   //
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////
			
		
		
		

	private:

		
		

	};

}  // namespace Kratos.

#endif // KRATOS_EMBEDDED_UTILITIES_INCLUDED  defined 


