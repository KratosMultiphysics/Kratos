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
//   Date:                $Date: 2009-01-13 16:40:58 $
//   Revision:            $Revision: 1.24 $
//
//


#if !defined(KRATOS_LEVEL_SET_UTILITIES_INCLUDED )
#define  KRATOS_LEVEL_SET_UTILITIES_INCLUDED



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
	template< class T, std::size_t dim >
			class PointDistance{
				public:
					double operator()( T const& p1, T const& p2 ){
						double dist = 0.0;
						for( std::size_t i = 0 ; i < dim ; i++){
							double tmp = p1[i] - p2[i];
							dist += tmp*tmp;
						}
						return sqrt(dist);
					}
			};
// 	template<std::size_t TDim>
	class LevelSetUtilities
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
//			inode->FastGetSolutionStepValue(IS_FREE_SURFACE) = 0.0;
			inode->FastGetSolutionStepValue(IS_BOUNDARY) = 0.0;
		}
		
		//looping over the elements to identify the fluid nodes 
		std::vector<unsigned int> node_is_fluid(dim+1);
		for( ModelPart::ElementsContainerType::iterator iel = base_model_part.ElementsBegin();
				   iel != base_model_part.ElementsEnd();
				   iel++)
		{
			Geometry< Node<3> >& geom = iel->GetGeometry();
			
			//identify nodes inside the fluid
			unsigned int n_fluid_nodes = 0;
			for(unsigned int i =0; i<dim+1; i++)
			{
				double dist = geom[i].FastGetSolutionStepValue(rDistanceVar);
				
				if(dist < acceptable_distance) //if we are inside the fluid
				{
					n_fluid_nodes++;	//number of fluid nodes in an element
					node_is_fluid[i] = 1;
				}
				else 
				{
					node_is_fluid[i] = 0;
				}
			}
			
			if( n_fluid_nodes > 0) // the element is a "fluid" element (that is the element as at least one fluid node;)
			{
				//add the element to the fluid domain
				fluid_model_part.Elements().push_back( *(iel.base() ) );
				
				for(unsigned int i =0; i<dim+1; i++)
				{
					//mark all of its nodes as fluid
					geom[i].FastGetSolutionStepValue(IS_FLUID) = 1.0;
	
					//mark as boundary the nodes for which distancei is > 0
					if( node_is_fluid[i] == 0)
						geom[i].FastGetSolutionStepValue(IS_BOUNDARY) = 1.0;
				}				
			}
				
		}
		
		//structure nodes can not be of boundary
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			if( inode->FastGetSolutionStepValue(IS_STRUCTURE) == 1.0)
			{
				inode->FastGetSolutionStepValue(IS_BOUNDARY) = 0.0;
			}
		}
			
		//add the nodes which were marked as fluid
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			if( inode->FastGetSolutionStepValue(IS_FLUID) == 1.0)
			{
				fluid_model_part.Nodes().push_back( *(inode.base() ) );
			}
		}
			
/*		//extrapolate the pressure inside the domain to the boundary
		array_1d<double,3> aux;
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			//KRATOS_WATCH("	ENTRAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA		********");
			if( inode->FastGetSolutionStepValue(IS_BOUNDARY) == 1.0)
			{
				WeakPointerVector< Node<3> >& neighb_nodes = inode->GetValue(NEIGHBOUR_NODES); 

				WeakPointerVector< Node<3> >::iterator extrapolation_node = neighb_nodes.begin();
				double min_dist = extrapolation_node->FastGetSolutionStepValue(DISTANCE);
				
				//find the node farther inside the domain
				for( WeakPointerVector< Node<3> >::iterator ineigh =	neighb_nodes.begin(); ineigh != neighb_nodes.end(); ineigh++) 
				{
					double dist = extrapolation_node->FastGetSolutionStepValue(DISTANCE);
					
					if(dist < min_dist)
					{
						min_dist = dist;
						extrapolation_node = ineigh;
					}					
				}	
				
				//vector oriented as the distance
				aux[0] = inode->Y() - extrapolation_node->Y();
				aux[2] = inode->X() - extrapolation_node->X();
				aux[1] = inode->Z() - extrapolation_node->Z();
				
				const array_1d<double,3>& proj = extrapolation_node->FastGetSolutionStepValue(PRESS_PROJ);
				double factor = inner_prod(aux,proj);
				double aux_norm = norm_2(proj); // * inode->FastGetSolutionStepValue(DENSITY) ;
				
				inode->FastGetSolutionStepValue(PRESSURE) = extrapolation_node->FastGetSolutionStepValue(PRESSURE) + factor * aux_norm ;
				
				
				
				
			}
		}	*/	
				
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
				inode->GetValue(rIsVisitedVar) = 0;
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
			double n_boundary_nodes = 0;
			for(unsigned int i =0; i<geom.size(); i++)
			{
			  n_boundary_nodes += double( geom[i].FastGetSolutionStepValue(IS_BOUNDARY));
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
/*			if(inode->FastGetSolutionStepValue(IS_BOUNDARY) != 1)
				inode->FastGetSolutionStepValue(rDistanceVar) = - fabs(dist);
			else
				inode->FastGetSolutionStepValue(rDistanceVar) = -dist;*/
		}
				
		KRATOS_CATCH("");
	}
	
	
	//***************************************************************************************
	//***************************************************************************************
	void ExtrapolateVelocities( 	ModelPart& base_model_part,
				    	Variable<double>& rDistanceVar,
	 				Variable< array_1d<double,3> >& rConvVelVar,
			 		Variable< array_1d<double,3> >& rVelVar,
					double extrapolation_distance 
				     )
	{
		KRATOS_TRY
		
// *******
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

			//*************
			// Bucket types
			   typedef Bucket< 3, PointType, PointVector > BucketType;
// 			   typedef Bins< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
// 			   typedef BinsDynamic< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > DynamicBins;
			//*************
			// DynamicBins;
			
			   typedef Tree< KDTreePartition<BucketType> > tree; 		//Kdtree;
// 			   typedef Tree< OCTreePartition<BucketType> > tree; 		//Octree;
// 			   typedef Tree< StaticBins > tree; 		     		//Binstree;
// 			   typedef Tree< KDTreePartition<StaticBins> > tree; 		//KdtreeBins;
// 			   typedef typename KdtreeBins::Partitions SubPartitions;
// 			   typedef Tree< OCTreePartition<StaticBins> > tree; 		//OctreeBins;
/*			
			   typedef Bins< 3, PointType, stdPointVector> stdBins;
			   typedef Tree< Bins<3,PointType,stdPointVector> > tree; 	//stdStaticBins;*/
// *******


				
		//create a list of the nodes on the boundary
		PointVector surface_nodes;
		array_1d<double,3> aux;			
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				inode != base_model_part.NodesEnd();
				inode++)	
		{
			if(inode->FastGetSolutionStepValue(IS_FLUID) == 0)
			{
				noalias(inode->FastGetSolutionStepValue(PRESS_PROJ)) = ZeroVector(3);
				noalias(inode->FastGetSolutionStepValue(rConvVelVar)) = ZeroVector(3);
				inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
				
				//i should NOT touch the boundary!!
				if(inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
					noalias(inode->FastGetSolutionStepValue(rVelVar)) = ZeroVector(3);
				
			}	
			else	//inside the fluid domain copy vel to convection_vel
			{
				if(inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
				{
					inode->FastGetSolutionStepValue(rConvVelVar) = inode->FastGetSolutionStepValue(rVelVar);
 				}
				else //in this case we need to extrapolate the velocities to the boundary
				{
					//make the average of the velocity of the neighbouring fluid nodes
					//if they are not is_structure node
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
				
// 					//remove normal component
// 					const array_1d<double,3>& normal = inode->FastGetSolutionStepValue(NORMAL);
// 					double An2 = inner_prod(normal,normal);
// 						
// 					double temp = inner_prod(aux,normal);
// 					noalias(aux) -= (temp / An2) * normal;
// 					
// 					//save the extrapolated convection velocity
 					inode->FastGetSolutionStepValue(rConvVelVar) = aux;
					
				}
			}
				
                        //select as an extrapolation source the nodes which are inside the domain but not the boundary nodes
                        double dist = inode->FastGetSolutionStepValue(rDistanceVar);
			if( dist <= 0.0 && dist > -0.5*extrapolation_distance)
				surface_nodes.push_back( *(inode.base() ) );
	
			//add boundary nodes to spatial database
			//NOTE: the node on the "advancing front" (boundary+structure) is excluded!!
/*			if(inode->FastGetSolutionStepValue(IS_BOUNDARY) == 1 && inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
				surface_nodes.push_back( *(inode.base() ) );*/
		}
KRATOS_WATCH(surface_nodes.size());
				
		//create a spatial database
		unsigned int bucket_size = 2;
		tree  nodes_tree2(surface_nodes.begin(),surface_nodes.end(),bucket_size);
		array_1d<double,3> aux_vect;	
		//extrapolate velocities from the closesest point
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				inode != base_model_part.NodesEnd();
				   inode++)	
		{
			double dist = inode->FastGetSolutionStepValue(rDistanceVar);
			
			if(dist > 0 && dist <= extrapolation_distance )
			{
// 				if(inode->FastGetSolutionStepValue(IS_BOUNDARY) == 0 )
				//&& inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
				{
					//find closest point from octtree
					double closest_distance = 0.0;
					PointTypePointer p_closest = nodes_tree2.SearchNearestPoint(*inode , closest_distance);

                                        array_1d<double,3>& conv_vel = inode->FastGetSolutionStepValue(rConvVelVar);
					noalias(conv_vel) = (p_closest)->FastGetSolutionStepValue(rVelVar);

// 					//make an average of the velocity of the closest nodes and of its neighbours
// 					noalias(aux_vect) = (p_closest)->FastGetSolutionStepValue(rVelVar);
// 					
// 					WeakPointerVector< Node<3> >& neighb_nodes =(*p_closest)->GetValue(NEIGHBOUR_NODES); 
// 					for( WeakPointerVector< Node<3> >::iterator ineigh =	neighb_nodes.begin(); ineigh != neighb_nodes.end(); ineigh++) 
// 					{
// 						if(ineigh->FastGetSolutionStepValue(IS_BOUNDARY) == 1)
// 							noalias(aux_vect) += ineigh->FastGetSolutionStepValue(rVelVar);	
// 					}
// 					aux_vect /= (1.0 + double(neighb_nodes.size() ) );
// 					
// 					//extrapolate convection velocity from closest point
// 					array_1d<double,3>& conv_vel = inode->FastGetSolutionStepValue(rConvVelVar);
// // 					conv_vel = (*p_closest)->FastGetSolutionStepValue(rVelVar);
// 					noalias(conv_vel) = aux_vect;
					
					//extrapolate velocity EXCEPT on STRUCTURE nodes
					if( inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
					{
						array_1d<double,3>& vel = inode->FastGetSolutionStepValue(rVelVar);
						vel = (p_closest)->FastGetSolutionStepValue(rVelVar);
					}
					
					
					
					//just to debug!!
//					inode->FastGetSolutionStepValue(TEMPERATURE) = (p_closest)->Id();
					
// 					//on the boundary remove normal component
// 					if( inode->FastGetSolutionStepValue(IS_STRUCTURE) == 1 )
// 					{
// 						const array_1d<double,3>& normal = inode->FastGetSolutionStepValue(NORMAL);
// 						double An2 = inner_prod(normal,normal);
// 						
// 						double temp = inner_prod(conv_vel,normal);
// 						noalias(conv_vel) -= (temp / An2) * normal;
// 					}
										
				}
			}
		}		
								
		KRATOS_CATCH("");
	}
	
	
	//***************************************************************************************
	//***************************************************************************************
	//generate a model part containing the extrapolation domain
	//fix the velocity on the boundary
	void ImplicitExtrapolation_PreProcess(ModelPart& base_model_part,
					      ModelPart& extrapolation_model_part,
					      Variable<double>& rDistanceVar,
	   				      Variable< array_1d<double,3> >& rVelVar,
	      				      Variable< array_1d<double,3> >& rConvVelVar,
    					      double extrapolation_distance    				
					     )
	{
		KRATOS_TRY
				
		//clearing the model part
		extrapolation_model_part.Nodes().clear();
		extrapolation_model_part.Elements().clear();
		extrapolation_model_part.Conditions().clear();

// 		for( ModelPart::NodesContainerType::iterator inode = extrapolation_model_part.NodesBegin();
// 				   inode != extrapolation_model_part.NodesEnd();
// 				   inode++)			
// 		{
// 			inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
// 	
// 		}

		//reset velocities of the nodes outside the extrapolation domain (outside the fluid)
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			double dist = inode->FastGetSolutionStepValue(rDistanceVar);
			if( dist > extrapolation_distance)
//			if( dist > 0.0 && inode->FastGetSolutionStepValue(IS_BOUNDARY) != 1.0)
			{
				noalias( inode->FastGetSolutionStepValue(rConvVelVar) ) = ZeroVector(3);
				inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
				noalias( inode->FastGetSolutionStepValue(PRESS_PROJ) ) = ZeroVector(3);
					
					//i should NOT touch the boundary!!
				if(inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
					noalias(inode->FastGetSolutionStepValue(rVelVar)) = ZeroVector(3);
					
			}
			
// if(inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
// 	noalias(inode->FastGetSolutionStepValue(rVelVar)) = ZeroVector(3);

		}
			
		//looping over the elements to identify the nodes of the extrapolation model part
		for( ModelPart::ElementsContainerType::iterator iel = base_model_part.ElementsBegin();
				   iel != base_model_part.ElementsEnd();
				   iel++)
		{
			Geometry< Node<3> >& geom = iel->GetGeometry();
			
			//identify nodes inside the fluid
                        double n_structure = 0;
			double n_test = 0;
			double n_fluid = 0.0;
			for(unsigned int i =0; i<geom.size(); i++)
			{
				double dist = geom[i].FastGetSolutionStepValue(rDistanceVar);
				
				if(dist >= 0.0 && dist <= extrapolation_distance) //if we are OUTside the fluid
				{
					n_test++;
				}
				
				n_fluid +=  static_cast<double>( geom[i].FastGetSolutionStepValue(IS_FLUID));
				n_structure +=  static_cast<double>(geom[i].FastGetSolutionStepValue(IS_STRUCTURE));
			}

                        //we have to avoid elements for which all of the nodes have the velocity fixed
			if( n_test > 0 
                            && n_fluid != double(geom.size() ) 
                            && (n_fluid+n_structure) != double(geom.size() ) 
                        ) // the element is a "fluid" element
			{
				//add the element to the fluid domain
				extrapolation_model_part.Elements().push_back( *(iel.base() ) );
				
				
				
			}		
		}
		
		//temporarily setting to -1 the ids of the nodes of the extrapolation model part
		for( ModelPart::ElementsContainerType::iterator iel = extrapolation_model_part.ElementsBegin();
				   iel != extrapolation_model_part.ElementsEnd();
				   iel++)
		{
			Geometry< Node<3> >& geom = iel->GetGeometry();
// KRATOS_WATCH(iel->Id());
			
			for(unsigned int i =0; i<geom.size(); i++)
				geom[i].FastGetSolutionStepValue(IS_FLUID) = -1.0;
		}
		
		
		
// KRATOS_WATCH("fixing ids");
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			double dist = inode->FastGetSolutionStepValue(rDistanceVar);
			
			
			
			if( inode->FastGetSolutionStepValue(IS_FLUID) == -1.0  )
			{
				//Identifying IS_BOUNDARY = -1 nodes***************************************************
				// if we are at the ''external limit'' of the extrapolation model part and the node is not a structural node set IS_BOUNDARY =-1
				if(dist >= extrapolation_distance && inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1.0 )
					inode->FastGetSolutionStepValue(IS_BOUNDARY) = -1.0;
				

				//If it is a IS_BOUNDARY node FIX VELOCITY 
				if( inode->FastGetSolutionStepValue(IS_BOUNDARY) == 1.0)
				{
					inode->GetValue(PRESSURE) = inode->FastGetSolutionStepValue(PRESSURE);
					noalias( inode->GetValue(PRESS_PROJ) ) = inode->FastGetSolutionStepValue(PRESS_PROJ);

					inode->Fix(VELOCITY_X);
					inode->Fix(VELOCITY_Y);
					inode->Fix(VELOCITY_Z);
					
					inode->Fix(FRACT_VEL_X);
					inode->Fix(FRACT_VEL_Y);
					inode->Fix(FRACT_VEL_Z);
				}
				//else setting to zero velocities
				else 
				{ 	
					noalias( inode->FastGetSolutionStepValue(VELOCITY) )=  ZeroVector(3);
					noalias( inode->FastGetSolutionStepValue(FRACT_VEL) )=  ZeroVector(3);
				}
			}
				
		}	
		
		//add the nodes which were marked as fluid
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			if( inode->FastGetSolutionStepValue(IS_FLUID) == -1.0)
			{
				extrapolation_model_part.Nodes().push_back( *(inode.base() ) );
			}
		}
		
		//reducing the density 
		for( ModelPart::NodesContainerType::iterator inode = extrapolation_model_part.NodesBegin();
				   inode != extrapolation_model_part.NodesEnd();
				   inode++)	
		{
// 			inode->FastGetSolutionStepValue(DENSITY) = 1.0;
			inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
			inode->FastGetSolutionStepValue(PRESS_PROJ) = ZeroVector(3);
			inode->FastGetSolutionStepValue(BODY_FORCE) = ZeroVector(3);

//                         if(inode->FastGetSolutionStepValue(IS_STRUCTURE) == 1)
//                             inode->Fix(PRESSURE);
		}
				
		KRATOS_CATCH("");
	}


	
	//***************************************************************************************
	//***************************************************************************************
	//free the velocity on the boundary
	//copy VELOCITY->CONVECTION_VELOCITY
	void ImplicitExtrapolation_PostProcess(ModelPart& base_model_part,
	     					Variable<double>& rDistanceVar,
     						Variable< array_1d<double,3> >& rConvVelVar,
      						Variable< array_1d<double,3> >& rVelVar,
	    					Vector& body_force
						)	
	{
		KRATOS_TRY
				
		array_1d<double,3> aux;
				
		//free nodes as needed
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
			inode++)	
		{

			//bringing back the density to the original value
//  			inode->FastGetSolutionStepValue(DENSITY) = 1000.0;
//ANTONIA change begin
			if( inode->FastGetSolutionStepValue(IS_FLUID) == 1.0 || inode->FastGetSolutionStepValue(IS_BOUNDARY) == 1.0 )
			{
				noalias(inode->FastGetSolutionStepValue(BODY_FORCE) ) = body_force;
                        }
//ANTONIA change end
			inode->Free(PRESSURE);
		
			
			//resetting IS_FLUID flag
			double dist = inode->FastGetSolutionStepValue(rDistanceVar);
			if( inode->FastGetSolutionStepValue(IS_BOUNDARY) == 1.0 )
				inode->FastGetSolutionStepValue(IS_FLUID) = 1.0;
			else 
			{
				
				if (dist > 0.0)
				{	
					inode->FastGetSolutionStepValue(IS_FLUID) = 0.0;
//Setting to zero pressure of the nodes of the extrapolation_model_part
					inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
					noalias(inode->FastGetSolutionStepValue(PRESS_PROJ)) = ZeroVector(3); //body_force; ???
				}
			}
			
			
			
			if( inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1.0 ) //not on wall
			{
				inode->FastGetSolutionStepValue(rConvVelVar) = inode->FastGetSolutionStepValue(rVelVar);
				
				//free the nodal velocities
				if( inode->FastGetSolutionStepValue(IS_BOUNDARY) == 1.0 )
				{
					inode->FastGetSolutionStepValue(IS_FLUID) = 1.0;
					
					inode->FastGetSolutionStepValue(PRESSURE) = inode->GetValue(PRESSURE);
					noalias( inode->FastGetSolutionStepValue(PRESS_PROJ) ) = inode->GetValue(PRESS_PROJ) ;
					
					inode->Free(VELOCITY_X);
					inode->Free(VELOCITY_Y);
					inode->Free(VELOCITY_Z);
					inode->Free(FRACT_VEL_X);
					inode->Free(FRACT_VEL_Y);
					inode->Free(FRACT_VEL_Z);
				}
			}
			else //on the boundaries
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

				//save the extrapolated convection velocity
				inode->FastGetSolutionStepValue(rConvVelVar) = aux;
			}
				
			//Set pressure and pressure projection to zero in the extrapolation_model_part	
// 			if( inode->FastGetSolutionStepValue(IS_FLUID) == 1.0 && inode->FastGetSolutionstepValue(IS_BOUNDARY != 1.0)
// 			{
// 				inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
// 				inode->FastGetSolutionStepValue(PRESS_PROJ) = 0.0;
// 			}

		}		
				

		KRATOS_CATCH("");
	}
		
	
	//***************************************************************************************
	//***************************************************************************************
	//generate a model part with all of the elements and nodes between min_dist and max_dist
// 	void GenerateModelPart(ModelPart& base_model_part, 
// 			        ModelPart& destination_model_part,
// 	  			Variable<double>& rDistanceVar,
//       				Variable<double>& rIsVisitedVar,
//       				double min_dist, 
// 	  			double max_dist   				
// 				)
// 	{
// 		KRATOS_TRY
// 				
// 		//clearing the model part
// 		destination_model_part.Nodes().clear();
// 		destination_model_part.Elements().clear();
// 		destination_model_part.Conditions().clear();
// 		
// 		//identify nodes to be included
// 		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
// 				   inode != base_model_part.NodesEnd();
// 				   inode++)	
// 		{
// 			double dist = inode->FastGetSolutionStepValue(rDistanceVar);
// // 			if( inode->Id() == 3047 |   inode->Id() == 3138)
// // 			{	KRATOS_WATCH(inode->Id())
// // 				KRATOS_WATCH(inode->FastGetSolutionStepValue(rDistanceVar))
// // 				KRATOS_WATCH(max_dist)
// // 				KRATOS_WATCH(min_dist)}
// 			
// 			if( dist >= min_dist && dist <= max_dist )
// 			{
// 				 
// 				destination_model_part.Nodes().push_back( *(inode.base() ) );
// 				inode->GetValue(rIsVisitedVar) = 1.0;
// 			}
// 			else
// 			{
// 				inode->GetValue(rIsVisitedVar) = 0.0;
// 			}
// 		}
// 			
// 		//looping over the elements to identify elements of interest
// 		for( ModelPart::ElementsContainerType::iterator iel = base_model_part.ElementsBegin();
// 				   iel != base_model_part.ElementsEnd();
// 				   iel++)
// 		{
// 			Geometry< Node<3> >& geom = iel->GetGeometry();
// 			
// 			//identify nodes inside the fluid
// 			double n_test = 0;
// 			for(unsigned int i =0; i<geom.size(); i++)
// 			{
// 				n_test += geom[i].GetValue(rIsVisitedVar);
// 			}
// 			
// 			if( n_test == geom.size()  ) // the element is a "fluid" element
// 			{
// 				//add the element to the fluid domain
// 				destination_model_part.Elements().push_back( *(iel.base() ) );			
// 			}
// 			
// 
// 		/*	if( iel->Id() == 22759)
// 			{	KRATOS_WATCH(iel->Id())
// 				KRATOS_WATCH(n_test)
// 			}*/		
// 
// 	
// 		}
// 				
// 		KRATOS_CATCH("");
// 	}
		
	//***************************************************************************************
	//***************************************************************************************
	//generate a model part with all of the elements and nodes between min_dist and max_dist
	void GenerateModelPart(ModelPart& base_model_part, 
			        ModelPart& destination_model_part,
	  			Variable<double>& rDistanceVar,
      				Variable<double>& rIsVisitedVar,
      				double min_dist, 
	  			double max_dist   				
				)
	{
		KRATOS_TRY
				
		//clearing the EXTRAPOLATION model part
		destination_model_part.Nodes().clear();
		destination_model_part.Elements().clear();
		destination_model_part.Conditions().clear();
		
		//identify nodes to be included
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			inode->GetValue(IS_VISITED) = 0;
		}

		for( ModelPart::ElementsContainerType::iterator iel = base_model_part.ElementsBegin();
		   iel != base_model_part.ElementsEnd();  iel++)
		{
			Geometry< Node<3> >& geom = iel->GetGeometry();
			
			
			//identify nodes inside the extrapolation domain
			double n_test = 0;
			for(unsigned int i =0; i<geom.size(); i++)
			{
				double dist = geom[i].FastGetSolutionStepValue(rDistanceVar);
				if(dist  >= min_dist && dist <= max_dist )
					n_test += 1;
			}
			

			if( n_test == geom.size()  ) // the element is a "fluid" element
			{
				//add the element to the fluid domain
				destination_model_part.Elements().push_back( *(iel.base() ) );	
				for(unsigned int i =0; i<geom.size(); i++)
				{
					geom[i].GetValue(IS_VISITED) = 1 ;
				}		
			}
			
		}

		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			if (inode->GetValue(IS_VISITED) == 1)
				destination_model_part.Nodes().push_back( *(inode.base() ) );
		}

		
				
		KRATOS_CATCH("");
	}

	//***************************************************************************************
	//***************************************************************************************
/*	void ExtrapolateVelocities( 	ModelPart& base_model_part,
					Variable<double>& rDistanceVar,
     Variable< array_1d<double,3> >& rConvVelVar,
     Variable< array_1d<double,3> >& rVelVar,
     double extrapolation_distance 
				  )
	{
		KRATOS_TRY
				
				typedef Node<3> PointType;

		typedef std::vector<PointType::Pointer>           PointVector;
		typedef std::vector<PointType::Pointer>::iterator PointIterator;
		typedef std::vector<double>               DistanceVector;
		typedef std::vector<double>::iterator     DistanceIterator;


		typedef Kratos::KDTreePartition<3, PointType, PointIterator, DistanceIterator> t_partition;
		typedef Kratos::Bucket<3, PointType, PointIterator, DistanceIterator, PointDistance<PointType,3 > > t_leaf;
		typedef Kratos::Tree<3, PointType, t_partition, t_leaf, PointIterator, DistanceIterator> kd_tree;

				
		//create a list of the nodes on the boundary
		PointVector surface_nodes;			
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			if(inode->FastGetSolutionStepValue(IS_FLUID) == 0)
			{
				noalias(inode->FastGetSolutionStepValue(VELOCITY)) = ZeroVector(3);
				inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
			}				
			if(inode->FastGetSolutionStepValue(IS_BOUNDARY) == 1)
				surface_nodes.push_back( *(inode.base() ) );
		}
				
		//create a spatial database
		unsigned int bucket_size = 2;
		kd_tree  nodes_tree2(surface_nodes.begin(),surface_nodes.end(),bucket_size);
				
		//extrapolate velocities from the closesest point
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			double dist = inode->FastGetSolutionStepValue(rDistanceVar);
			
			if(dist > 0 && dist <= extrapolation_distance )
			{
				if(inode->FastGetSolutionStepValue(IS_BOUNDARY) == 0
							       && inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
				{
					//find closest point from octtree
					double closest_distance = 0.0;
					tree::PointerType p_closest = nodes_tree2.SearchNearestPoint(*inode , closest_distance);
					
					//extrapolate velocity from closest point
					
					inode->FastGetSolutionStepValue(VELOCITY) = (p_closest)->FastGetSolutionStepValue(VELOCITY);
					
					
					
				}
			}
// 			else if (inode->FastGetSolutionStepValue(IS_FLUID) != 1)
// 			{
// 				noalias(inode->FastGetSolutionStepValue(VELOCITY)) = ZeroVector(3);
// 			}
		}		
				
				
				
		KRATOS_CATCH("");
	}*/
	
	
	
	//***************************************************************************************
	//***************************************************************************************
	void ExtrapolateVelocitiesByLayer( 	ModelPart& base_model_part,
					Variable<double>& rDistanceVar,
     					Variable< array_1d<double,3> >& rConvVelVar,
     					Variable< array_1d<double,3> >& rVelVar,
     					unsigned int extrapolation_layers
				  )
	{
		KRATOS_TRY
				
		typedef Node<3> PointType;
		typedef PointerVector<PointType >           PointVector;
		typedef PointVector::iterator PointIterator;
		
		//generate a container with the layers to be extrapolated
		std::vector< PointVector > layers(extrapolation_layers);
		
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			inode->GetValue(IS_VISITED) = 0.0;
		}
		
		//generate the convection var inside the fluid and fill the first layer
		array_1d<double,3> aux;	
				

		
		for( ModelPart::NodesContainerType::iterator inode = base_model_part.NodesBegin();
				   inode != base_model_part.NodesEnd();
				   inode++)	
		{
			if(inode->FastGetSolutionStepValue(IS_FLUID) == 0)
			{
				noalias(inode->FastGetSolutionStepValue(PRESS_PROJ)) = ZeroVector(3);
				noalias(inode->FastGetSolutionStepValue(rConvVelVar)) = ZeroVector(3);
				inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
				
				//i should NOT touch the bounary!!
				if(inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
					noalias(inode->FastGetSolutionStepValue(rVelVar)) = ZeroVector(3);
				
			}	
			else	//inside the fluid domain copy vel to convection_vel
			{
				if(inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
				{
					inode->FastGetSolutionStepValue(rConvVelVar) = inode->FastGetSolutionStepValue(rVelVar);
				}
				else //in this case we need to extrapolate the velocities to the boundary
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
				
// 					//remove normal component
// 					const array_1d<double,3>& normal = inode->FastGetSolutionStepValue(NORMAL);
// 					double An2 = inner_prod(normal,normal);
					// 						
// 					double temp = inner_prod(aux,normal);
// 					noalias(aux) -= (temp / An2) * normal;
					// 					
// 					//save the extrapolated convection velocity
					inode->FastGetSolutionStepValue(rConvVelVar) = aux;
					
				}
			}
				
			//NOTE: the node on the "advancing front" (boundary+structure) is excluded!!
			//add boundary nodes to spatial database
			   if(inode->FastGetSolutionStepValue(IS_BOUNDARY) == 1 && inode->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
			{
				layers[0].push_back( *(inode.base() ) );
				inode->GetValue(IS_VISITED) = 1.0;
			}
				
		}
		
		//fill the following layers by neighbour relationships
		//each layer fills the following
		for(unsigned int il = 0; il<extrapolation_layers-1; il++)
		{
			for( PointIterator iii=(layers[il]).begin(); iii!=(layers[il]).end(); iii++)
			{
				WeakPointerVector< Node<3> >& neighb_nodes = iii->GetValue(NEIGHBOUR_NODES); 
				for(WeakPointerVector< Node<3> >::iterator jjj=neighb_nodes.begin(); jjj !=neighb_nodes.end(); jjj++) 
				{ 

					if( jjj->FastGetSolutionStepValue(IS_FLUID) == 0 &&
					jjj->GetValue(IS_VISITED) == 0.0 )
					{
						layers[il+1].push_back( Node<3>::Pointer( *(jjj.base() ) ) );
						jjj->GetValue(IS_VISITED) = double(il+2.0);
					}
				}
			}
		}
		
		for(unsigned int il = 0; il<extrapolation_layers; il++)
		{
			std::cout << "layer = " << il+1<< " nodes " << layers[il].size() <<std::endl;
		}
		
		//now loop over the different layers and obtain the convection velocity by making an average 
		//of the neighbours of lower order
		for(unsigned int il = 1; il<extrapolation_layers; il++)
		{
			for( PointIterator iii=layers[il].begin(); iii!=layers[il].end(); iii++)
			{
				noalias(aux) = ZeroVector(3);
				double avg_number = 0.0;
				
				WeakPointerVector< Node<3> >& neighb_nodes = iii->GetValue(NEIGHBOUR_NODES); 
				for(WeakPointerVector< Node<3> >::iterator i=neighb_nodes.begin(); 				i !=neighb_nodes.end(); i++) 
				{ 
					if(i->GetValue(IS_VISITED) < il+1 && i->GetValue(IS_VISITED) > 0)
					{
						noalias(aux) += i->FastGetSolutionStepValue(rConvVelVar);
						avg_number += 1.0;
					}
				} 
                                if(avg_number != 0.0)
				    aux /= avg_number;
				noalias( iii->FastGetSolutionStepValue(rConvVelVar) ) = aux;
			}
		}		
				
								
		KRATOS_CATCH("");
	}

	//***************************************************************************************
	//***************************************************************************************
	//ensure that each connected component of the extrapolation model part has at least 1 
	//free surface or fixed pressure
	void ApplyMinimumExtrapolationPressureFix( 	ModelPart& extrapolation_model_part
					 )
	{
		KRATOS_TRY
				
		typedef Node<3> PointType;
		typedef PointerVector<PointType >  PointVector;
		typedef PointVector::iterator PointIterator;
		
		//generate a container with the layers to be extrapolated
		PointVector work_array;
		work_array.reserve(extrapolation_model_part.Nodes().size());
		
		//reset the counter
		for( ModelPart::NodesContainerType::iterator inode = extrapolation_model_part.NodesBegin();
				   inode != extrapolation_model_part.NodesEnd();
				   inode++)	
		{
			inode->GetValue(IS_VISITED) = 0.0;
// 			inode->Free(PRESSURE);
		}
		
		double connected_component_index = 1;
// 		PointIterator current_it = work_array.begin();
		unsigned int aux_index = 0;
		unsigned int number_of_connected_components = 0;
		for( ModelPart::NodesContainerType::iterator inode = extrapolation_model_part.NodesBegin();
				   inode != extrapolation_model_part.NodesEnd();
				   inode++)	
		{
			
			if(inode->GetValue(IS_VISITED) == 0.0) //this node is not yet visited - it is in a new 
			{
				number_of_connected_components++;
				unsigned int n_free_boundaries = 0;
				
				//add the node to the work array - this node will be the "seed" of the partition
				inode->GetValue(IS_VISITED) = connected_component_index;
				work_array.push_back( Node<3>::Pointer( *(inode.base()) ) ); 
				
				while(work_array.begin()+aux_index < work_array.end())
				{
					
					if( (work_array.begin()+aux_index)->FastGetSolutionStepValue(IS_BOUNDARY) == -1.0)
					{
						n_free_boundaries++;
/*						KRATOS_WATCH(n_free_boundaries);		*/
					}
					
					//add its neighbours to the work array and mark them
					WeakPointerVector< Node<3> >& neighb_nodes = (work_array.begin()+aux_index)->GetValue(NEIGHBOUR_NODES); 
					for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++) 
					{ 
						if(i->GetValue(IS_VISITED) == 0.0) 
						{
							//add the node to the work array
							i->GetValue(IS_VISITED) = connected_component_index;
							work_array.push_back( Node<3>::Pointer(*(i.base()) ) ); 
						}
						
					}
						
					//update the current counter
					aux_index++;
				}
// 				KRATOS_WATCH("ccc")			
				
				//fix at least one pressure if needed (the pressure will be fixed on the "seed" node
				if( n_free_boundaries==0)
				{
					std::cout << "minimal pressure condition prescribed on node " << inode->Id() << std::endl;
					inode->Fix(PRESSURE);
				}
				
				//increase the component index counter
				connected_component_index++;
			}
			
			
		}
		std::cout << "number of connected components" <<  number_of_connected_components << std::endl;
						
		KRATOS_CATCH("");
	}
	
	private:
		



	};

}  // namespace Kratos.

#endif // KRATOS_LEVEL_SET_UTILITIES_INCLUDED  defined 


