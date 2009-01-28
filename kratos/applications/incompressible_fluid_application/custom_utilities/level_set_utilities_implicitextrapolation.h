/*
==============================================================================
KratosIncompressibleFluidApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
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
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2009-01-14 16:24:49 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_LEVEL_SET_UTILITIES_IMPLICIT_EXTRAPOLATION_INCLUDED )
#define  KRATOS_LEVEL_SET_UTILITIES_IMPLICIT_EXTRAPOLATION_INCLUDED



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
//#include "geometries/tetrahedra_3d_4.h"
#include "incompressible_fluid_application.h"



#include "spatial_containers/spatial_containers.h"


namespace Kratos
{
	
	class LevelSetUtilitiesImplicitExtrapolation
	{
	public:

	//***************************************************************************************
	//***************************************************************************************
	void RegenerateFluidModelPart( ModelPart& base_model_part,
				       ModelPart& fluid_model_part,
	   				Variable<double>& rDistanceVar,
                                        const double acceptable_distance, // = 0.0,
					unsigned int dim
				     )
	{
		KRATOS_TRY
				
		//clearing the model part
		fluid_model_part.Nodes().clear();
		fluid_model_part.Elements().clear();
		fluid_model_part.Conditions().clear();
		
		//clear flags
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			inode->FastGetSolutionStepValue(IS_FLUID) = 0.0;
			inode->FastGetSolutionStepValue(IS_BOUNDARY) = 0.0;
		}	
		
		//looping over the elements to identify
		std::vector<unsigned int> node_is_fluid(dim+1);
		std::vector<unsigned int> node_is_acceptable(dim+1);
		for( ModelPart::ElementsContainerType::iterator iel = base_model_part.ElementsBegin();
				   iel != base_model_part.ElementsEnd();
				   iel++)
		{
			Geometry< Node<3> >& geom = iel->GetGeometry();
			
			//identify nodes inside the fluid
			unsigned int n_fluid_nodes = 0;
			unsigned int n_acceptable = 0;
			for(unsigned int i =0; i<dim+1; i++)
			{
				double dist = geom[i].FastGetSolutionStepValue(rDistanceVar);
				
				if(dist < acceptable_distance) //if we are inside the fluid
				{
					
					n_acceptable++;
					
					node_is_acceptable[i] = 1;
					if(dist < 0.0)
					{
						n_fluid_nodes++;
						node_is_fluid[i] = 1;
					}
					else
						node_is_fluid[i] =0;
				}
				else 
				{
					node_is_acceptable[i] = 0;
					node_is_fluid[i] = 0;
				}
			}
			
			//mark:
			//IS_FLUID = 1, IS_BOUNDARY = 1 inside the proper fluid domain
			//IS_FLUID = -1, IS_BOUNDARY = 0 on the extension
			if( n_acceptable > 0) // the element is a inside the "fluid" or inside its extension
			{
				//add the element to the fluid domain
				fluid_model_part.Elements().push_back( *(iel.base() ) );
				
				if(n_fluid_nodes > 0) //element of fluid
				{
					for(unsigned int i =0; i<dim+1; i++)
					{
						//mark all of its nodes as fluid
						geom[i].FastGetSolutionStepValue(IS_FLUID) = 1.0;
		
						//mark as boundary the nodes for which distancei is > 0
						if( node_is_fluid[i] == 0)
							geom[i].FastGetSolutionStepValue(IS_BOUNDARY) = 1.0;
					}	
				}
				else //element in the extension part of the domain
				{
					for(unsigned int i =0; i<dim+1; i++)
					{
						//mark all of its nodes as fluid
						if(geom[i].FastGetSolutionStepValue(IS_BOUNDARY) != 1.0)
							geom[i].FastGetSolutionStepValue(IS_FLUID) = -1.0;
		
						//mark as boundary the nodes for which distancei is > 0
						if( node_is_acceptable[i] == 0)
							geom[i].FastGetSolutionStepValue(IS_BOUNDARY) = -1.0;
					}	
				}					
			}
				
		}
		
	
		//add the nodes which were marked as fluid 
		//mark extrapolation part with IS_FLUID=0
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			if( inode->FastGetSolutionStepValue(IS_FLUID) == 1.0)
			{
				fluid_model_part.Nodes().push_back( *(inode.base() ) );
			}
			else if( inode->FastGetSolutionStepValue(IS_FLUID) == -1.0 )
			{
/*				inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
				inode->FastGetSolutionStepValue(PRESS_PROJ) = ZeroVector(3);*/
				inode->FastGetSolutionStepValue(IS_FLUID) = 0.0;
				fluid_model_part.Nodes().push_back( *(inode.base() ) );
			}		
		}	
				
		KRATOS_CATCH("");
	}
	
	//***************************************************************************************
	//***************************************************************************************
	void FluidDistanceComputation_FromBoundary( 	ModelPart::NodesContainerType& rNodes,
					Variable<double>& rIsVisitedVar,
     					Variable<double>& rDistanceVar
			       )
	{
		KRATOS_TRY
				
	
				
		for( ModelPart::NodesContainerType::iterator inode = rNodes.begin(); inode != rNodes.end(); inode++)	
		{
			
			if(inode->FastGetSolutionStepValue(IS_BOUNDARY) == 1)
			{
				double dist =  inode->FastGetSolutionStepValue(rDistanceVar) ;
				inode->FastGetSolutionStepValue(rDistanceVar) = -dist;
				inode->GetValue(rIsVisitedVar) = 1;
			}	
			else
			{	
				double dist =  inode->FastGetSolutionStepValue(rDistanceVar) ;
				if( dist > 0.0 )
					inode->GetValue(rIsVisitedVar) = 1;
				else
					inode->GetValue(rIsVisitedVar) = 0;
			}
		}
				
		KRATOS_CATCH("");
	}
	
	//***************************************************************************************
	//***************************************************************************************
	void PrepareForInternalFluidDistanceComputation( 	ModelPart& fluid_model_part,
			Variable<double>& rIsVisitedVar,
  			 Variable<double>& rDistanceVar
						       )
	{
		KRATOS_TRY
				
		//reset is_visited on all of the fluid nodes
		for( ModelPart::NodesContainerType::iterator inode = fluid_model_part.NodesBegin();
		inode != fluid_model_part.NodesEnd();
		inode++)	
		{
			inode->GetValue(rIsVisitedVar) = 0;
		}
				
		//select all the elements with at least one node on the boundary
		for( ModelPart::ElementsContainerType::iterator iel = fluid_model_part.ElementsBegin();
			iel != fluid_model_part.ElementsEnd();
			iel++)
		{
			Geometry< Node<3> >& geom = iel->GetGeometry();
			
			//identify elements with at least one node on boundary
			double  n_boundary_nodes = 0;
			for(unsigned int i =0; i<geom.size(); i++)
			{
				n_boundary_nodes += geom[i].FastGetSolutionStepValue(IS_BOUNDARY);
			}
			 
			//fix them and change sign to the distance
			if(n_boundary_nodes > 0)
			{
				for(unsigned int i =0; i<geom.size(); i++)
				{
					geom[i].GetValue(rIsVisitedVar) = 1;
					double dist =  geom[i].FastGetSolutionStepValue(rDistanceVar) ;
					geom[i].FastGetSolutionStepValue(rDistanceVar) = -dist;
				}
					

			}
		}
			
				
		KRATOS_CATCH("");
	}
	
	//***************************************************************************************
	//***************************************************************************************
	void MarkNodesAsVisited( 	ModelPart::NodesContainerType& rNodes,
					Variable<double>& rIsVisitedVar
				  )
	{
		KRATOS_TRY
				
	
				
		//extrapolate velocities from the closesest point
		for( ModelPart::NodesContainerType::iterator inode = rNodes.begin();
				inode != rNodes.end();
				inode++)	
		{
			if( inode->FastGetSolutionStepValue(IS_FLUID) == 1.0)
			{
				inode->GetValue(rIsVisitedVar) = 1;
			}
			else
				inode->GetValue(rIsVisitedVar) = 0;
		}
				
		KRATOS_CATCH("");
	}
	
	//***************************************************************************************
	//***************************************************************************************
	void SetDistanceToNegative( 	ModelPart::NodesContainerType& rNodes,
					Variable<double>& rDistanceVar
			       )
	{
		KRATOS_TRY
				
	
		for( ModelPart::NodesContainerType::iterator inode = rNodes.begin();
		inode != rNodes.end();
		inode++)	
		{
			double dist =  inode->FastGetSolutionStepValue(rDistanceVar) ;
			
			inode->FastGetSolutionStepValue(rDistanceVar) = -dist;
		}
				
		KRATOS_CATCH("");
	}
	
	

	
	//***************************************************************************************
	//***************************************************************************************
	//generate a model part containing the extrapolation domain
	//fix the velocity on the boundary
	void ExtrapolateVelocities(ModelPart& base_model_part,
					      ModelPart& fluid_model_part,
	   				     const double acceptable_distance,
					      Variable<double>& rDistanceVar,
	   				      Variable< array_1d<double,3> >& rVelVar,
	      				      Variable< array_1d<double,3> >& rConvVelVar
					     )
	{
		KRATOS_TRY
				
		array_1d<double,3> aux;	
						
		//reset velocities
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			if(inode->FastGetSolutionStepValue(IS_BOUNDARY) != -1.0 && inode->FastGetSolutionStepValue(DISTANCE) > acceptable_distance)
			{
				noalias( inode->FastGetSolutionStepValue(rConvVelVar) ) = ZeroVector(3);
				inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
				noalias( inode->FastGetSolutionStepValue(PRESS_PROJ) ) = ZeroVector(3);
					
					//i should NOT touch the bounary!!
				if(inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
					noalias(inode->FastGetSolutionStepValue(rVelVar)) = ZeroVector(3);
					
			}
		}

		
		
		
		
		
		
		
		
		
		
		
		
		for( ModelPart::NodesContainerType::iterator inode = fluid_model_part.NodesBegin();
				   inode != fluid_model_part.NodesEnd();
				   inode++)	
		{
//			KRATOS_WATCH("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb");
			if( inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1.0 ) //not on wall
			{
				inode->FastGetSolutionStepValue(rConvVelVar) = inode->FastGetSolutionStepValue(rVelVar);
			}
			else //extrapolate to the walls
			{
				//make the average of the velocity of the neighbouring fluid nodes
				double n_fluid_neigh = 0.0;
				noalias(aux) = ZeroVector(3);
				WeakPointerVector< Node<3> >& neighb_nodes = inode->GetValue(NEIGHBOUR_NODES); 
				for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++) 
				{ 
					if(i->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
					{
						noalias(aux) += i->FastGetSolutionStepValue(rVelVar);
						n_fluid_neigh += 1.0;
					}
				} 
			
				if(n_fluid_neigh != 0.0)
					aux /= n_fluid_neigh;
				
				//remove the component normal to the wall
/*				const array_1d<double,3>& normal = inode->FastGetSolutionStepValue(NORMAL);
				double An2 = inner_prod(normal,normal);
										
				double temp = inner_prod(aux,normal);
				noalias(aux) -= (temp / An2) * normal;*/
				
				//save the extrapolated convection velocity
				inode->FastGetSolutionStepValue(rConvVelVar) = aux;
			}
			
			
		}		
			
		//set to zero the velocity on the extension part of the domain
		for( ModelPart::NodesContainerType::iterator inode = fluid_model_part.NodesBegin();
				   inode != fluid_model_part.NodesEnd();
				   inode++)	
		{
//			KRATOS_WATCH("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb");
			if( inode->FastGetSolutionStepValue(IS_FLUID) == 0.0 && inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1.0) //not on wall
			{
				inode->FastGetSolutionStepValue(rVelVar) = ZeroVector(3);
			}
		}		
	
		KRATOS_CATCH("");
	}
	
	void ApplyFluidProperties(ModelPart::NodesContainerType& rNodes, 
				  bool fix_pressures,
				  const double fluid1_density,
      				  const double fluid1_viscosity,
				  const double fluid2_density,
      				  const double fluid2_viscosity,
				  Vector& body_force
				 )
	{
		//apply the material density as liked
		for(ModelPart::NodesContainerType::iterator inode = rNodes.begin(); inode!=rNodes.end();  inode++)
		{
			if(inode->FastGetSolutionStepValue(IS_FLUID) == 1)
			{
				inode->FastGetSolutionStepValue(DENSITY) = fluid1_density;
				inode->FastGetSolutionStepValue(VISCOSITY) = fluid1_viscosity;	
//				noalias(inode->FastGetSolutionStepValue(BODY_FORCE)) = body_force;	
			}
			else
			{
				inode->FastGetSolutionStepValue(DENSITY) = fluid2_density;
				inode->FastGetSolutionStepValue(VISCOSITY) = fluid2_viscosity;
//				noalias(inode->FastGetSolutionStepValue(BODY_FORCE)) = ZeroVector(3);	
			}
		}
		
		//apply the pressure fixity if required
		if(fix_pressures == true)
		{
KRATOS_WATCH("fixing pressures");
			for(ModelPart::NodesContainerType::iterator inode = rNodes.begin(); inode!=rNodes.end();  inode++)
			{
				if(inode->FastGetSolutionStepValue(IS_BOUNDARY) == -1)
				{
// 					std::cout << "fixed pressure -- ID = " << inode->Id() << std::endl;
					inode->Fix(PRESSURE);
					inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
				}
				else
				{
					inode->Free(PRESSURE);
				}
			}		
		}
		
		
		//extend fluid velocities outside of the fluid1 part - this will be corrected during the solution
		
// 		//******************************************************************
// 		//extrapolate velocities "explicit" by vicinity
// ***
			//defintions for spatial search
			typedef Node<3> PointType;
			typedef Node<3>::Pointer PointTypePointer;
			typedef std::vector<PointType::Pointer>           PointVector;
			typedef std::vector<PointType::Pointer>::iterator PointIterator;
			typedef std::vector<double>               DistanceVector;
			typedef std::vector<double>::iterator     DistanceIterator;


			//creating an auxiliary list for the new nodes 
			PointVector list_of_new_nodes;

// 			KRATOS_WATCH("STARTING KDTREE CONSTRUCTION");
			//starting calculating time of construction of the kdtree
// 			boost::timer kdtree_construction;

			//*************
			// Bucket types
// 			   typedef Bucket< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
// 			   typedef Bins< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
// 			   typedef BinsDynamic< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > DynamicBins;
			//*************
			// DynamicBins;
			
// 			   typedef Tree< KDTreePartition<BucketType> > tree; 		//Kdtree;
// 			   typedef Tree< OCTreePartition<BucketType> > tree; 		//Octree;
// 			   typedef Tree< StaticBins > tree; 		     		//Binstree;
// 			   typedef Tree< KDTreePartition<StaticBins> > tree; 		//KdtreeBins;
// 			   typedef typename KdtreeBins::Partitions SubPartitions;
// 			   typedef Tree< OCTreePartition<StaticBins> > tree; 		//OctreeBins;
/*			
			   typedef Bins< TDim, PointType, stdPointVector> stdBins;
			   typedef Tree< Bins<TDim,PointType,stdPointVector> > tree; 	//stdStaticBins;*/
// ***

// 		//create a list of the nodes on the boundary
// 		PointVector surface_nodes;
// 		array_1d<double,3> aux;			
// 		for( ModelPart::NodesContainerType::iterator inode = rNodes.begin();
// 				   inode != rNodes.end();
// 				   inode++)	
// 		{
// 		
// 			//add boundary nodes to spatial database
// 			//NOTE: the node on the "advancing front" (boundary+structure) is excluded!!
// 			if(inode->FastGetSolutionStepValue(IS_BOUNDARY) == 1 && inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
// 				surface_nodes.push_back( *(inode.base() ) );
// 		}
// 				
// 		//create a spatial database
// 		unsigned int bucket_size = 2;
// 		kd_tree  nodes_tree2(surface_nodes.begin(),surface_nodes.end(),bucket_size);
// 		array_1d<double,3> aux_vect;	
// 		//extrapolate velocities from the closesest point
// 		for( ModelPart::NodesContainerType::iterator inode = rNodes.begin();
// 				   inode != rNodes.end();
// 				   inode++)	
// 		{
// 			
// 			if(inode->FastGetSolutionStepValue(IS_FLUID) == 0 )
// 			{
// 				if(inode->FastGetSolutionStepValue(IS_BOUNDARY) == 0 )
// 				{
// 					//find closest point from octtree
// 					double closest_distance = 0.0;
// 					PointTypePointer p_closest = nodes_tree2.SearchNearestPoint(*inode , closest_distance);
// 					
// 					//extrapolate velocity EXCEPT on STRUCTURE nodes
// 					if( inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
// 					{
// 						array_1d<double,3>& vel = inode->FastGetSolutionStepValue(VELOCITY);
// 						vel = (*it_closest)->FastGetSolutionStepValue(VELOCITY);
// 					}
// 										
// 				}
// 			}
// 		}
		//******************************************************************
			
	}
	
	private:
		



	};

}  // namespace Kratos.

#endif // KRATOS_LEVEL_SET_UTILITIES_IMPLICIT_EXTRAPOLATION_INCLUDED  defined 


