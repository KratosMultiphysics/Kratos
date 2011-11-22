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
		
		void SaveInterfaceElemsModelPart(ModelPart& reduced_model_part, ModelPart& interface_model_part)
		{
		interface_model_part.Conditions().clear();
		interface_model_part.Elements().clear();
		interface_model_part.Nodes().clear();
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
		       interface_model_part.AddElement(*(im.base()));		  
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
				  reduced_model_part.AddNode(*(in.base()));
				  
				
		      }
		  
		}
		
		///////////////////////////////////////////////////////////////////////////
		////////	SUBDOMAIN DISABLING			//////////////////
		///////////////////////////////////////////////////////////////////////////
		void DisableSubdomain(ModelPart& full_model_part, ModelPart& reduced_model_part)
		{
		KRATOS_TRY
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
				/*
				else 
					{
					in->FastGetSolutionStepValue(VELOCITY_X)=0.0;
					in->FastGetSolutionStepValue(VELOCITY_Y)=0.0;
					in->FastGetSolutionStepValue(VELOCITY_Z)=1.0;
					}
				*/
				}

		}

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

#endif // KRATOS_EMBEDDED_UTILITIES_INCLUDED  defined 


