/*
==============================================================================
KratosTestApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2010
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

The  above  copyright  notice  and  this permission  notice  sKRATOS_WATCH(disp);hall  be
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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:31 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_LAGRANGIAN_PFEM2_UTILITIES_INCLUDED )
#define  KRATOS_LAGRANGIAN_PFEM2_UTILITIES_INCLUDED

#define PRESSURE_ON_EULERIAN_MESH
#define USE_FEW_PARTICLES

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
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"


#include <boost/timer.hpp>
#include "utilities/timer.h"

#ifdef _OPENMP
#include "omp.h"
#endif



namespace Kratos
{

template<std::size_t TDim> class LagrangianPFEM2Utilities
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LagrangianPFEM2Utilities<TDim>);
    
   
    //**********************************************************************************************
    //**********************************************************************************************
    ///this function detects the inlet nodes as nodes which have a fixed velocity and a velocity applied which is different from zero
    ///@param rModelPart the model part on which we work
    ///@param rInletNodes is a container which contains the list of inlet nodes
    void DetectInletAndOutlet(ModelPart& rModelPart, ModelPart::NodesContainerType& rInletNodes, ModelPart::NodesContainerType& rOutletNodes)
    {    
 	rInletNodes = ModelPart::NodesContainerType();
	rOutletNodes = ModelPart::NodesContainerType();
        for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
	{
	    const array_1d<double,3>& vel = it->FastGetSolutionStepValue(VELOCITY);
	    double vnorm = norm_2(vel);
	    if(it->IsFixed(VELOCITY_X) == true && vnorm>1e-10)
	    {
	        rInletNodes.push_back( *(it.base() ));
	    }
	    if(it->IsFixed(PRESSURE) == true && it->IsFixed(VELOCITY_X) == false)
	    {
	        rOutletNodes.push_back( *(it.base() ));
	    }
	}
	
    }
    
    //**********************************************************************************************
    //**********************************************************************************************
    ///this function moves the mesh as xn+1 = xn + vn*dt and sets the mesh velocity to vn
    ///@param rModelPart the model part on which we work
    void MoveMesh_ForwardEuler(ModelPart& rModelPart)
    {      
        const double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
	
        for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
	{
	    const array_1d<double,3>& vn = it->FastGetSolutionStepValue(VELOCITY,1);
	    array_1d<double,3>& dn = it->FastGetSolutionStepValue(DISPLACEMENT,1);
	    array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
	    noalias(dn1) = dn;
	    noalias(dn1) += dt*vn;
	    
	    noalias(it->Coordinates()) = it->GetInitialPosition();
	    noalias(it->Coordinates()) += dn1;
	    
	    noalias(it->FastGetSolutionStepValue(MESH_VELOCITY)) = it->FastGetSolutionStepValue(VELOCITY,1);	    
	}
    }
    
    //**********************************************************************************************
    //**********************************************************************************************
    ///this function moves the mesh as xn+1 = xn + vn*dt and sets the mesh velocity to vn
    ///@param rModelPart the model part on which we work
    void MoveMesh_Streamlines(ModelPart& rModelPart, unsigned int substep)
    {      
        const double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
	
	//generate search data structure
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();
	
	//do movement
	array_1d<double, 3 > veulerian;
        array_1d<double, 3 > acc_particle;
        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        const int nparticles = rModelPart.Nodes().size();
//KRATOS_WATCH("551")
        #pragma omp parallel for firstprivate(results,N,veulerian,acc_particle)
        for (int i = 0; i < nparticles; i++)
        {
            unsigned int substep = 0;
            unsigned int subdivisions = 1;
            double small_dt = dt;
	    
	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
	    Node < 3 > ::Pointer pparticle = *(iparticle.base());
	    
	    bool do_move = true;
	    if(iparticle->IsFixed(VELOCITY_X) == true )
	    {
		if(iparticle->FastGetSolutionStepValue(VELOCITY_X) == 0.0) do_move = false;
		else do_move = true;
	    }
	    
	    
	    if( do_move == true  ) //note that we suppose the velocity components to be all fixed
	    {
		
		array_1d<double,3> old_position = pparticle->Coordinates();
		array_1d<double,3> current_position = pparticle->Coordinates();
		noalias(iparticle->GetInitialPosition()) = old_position;
		iparticle->FastGetSolutionStepValue(DISPLACEMENT,1) = ZeroVector(3);
		

		while(substep++ < subdivisions)
		{
		    
		    typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

		    Element::Pointer pelement;

		    bool is_found = SearchStructure.FindPointOnMesh(current_position, N, pelement, result_begin, max_results);
		    (iparticle)->Set(TO_ERASE, true);

		    if (is_found == true)
		    {
			Geometry< Node < 3 > >& geom = pelement->GetGeometry();

			noalias(veulerian) = N[0] * geom[0].FastGetSolutionStepValue(VELOCITY, 1);
			for (unsigned int k = 1; k < geom.size(); k++)
			    noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY, 1);

			//compute adaptive subdivisions
			if(substep == 1)
			{
			    //compute h
			    double h = N[0] * geom[0].FastGetSolutionStepValue(NODAL_H);
			    for (unsigned int k = 1; k < geom.size(); k++)
				h += N[k] * geom[k].FastGetSolutionStepValue(NODAL_H);

			    //compute number of subdivisions needed
			    const unsigned int min_subdivisions = 3;
			    const unsigned int max_subdivisions = 20;
			    double v = norm_2(veulerian);
			    double subdivisions = double(floor(2*dt*v/h));
			    subdivisions = (subdivisions<min_subdivisions) ? min_subdivisions : (subdivisions>max_subdivisions) ? max_subdivisions : subdivisions;

			    //compute subdivisions time step
			    small_dt = dt / subdivisions;

			    //KRATOS_WATCH(subdivisions)

			}

			//move according to the streamline

    //                     array_1d<double, 3 > & disp = (iparticle)->FastGetSolutionStepValue(DISPLACEMENT);

			noalias(current_position) += small_dt*veulerian;

			(pparticle)->Set(TO_ERASE, false);


    //KRATOS_WATCH("619")
		    }
		}
		
		//update the displacement BUT DO NOT OVERWRITE THE POSITION!!                
		iparticle->FastGetSolutionStepValue(DISPLACEMENT) = current_position - iparticle->GetInitialPosition(); 

	    }
            
    
        }
	
	//compute mesh velocity
        for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
	{
// 	    const array_1d<double,3>& vn = it->FastGetSolutionStepValue(VELOCITY,1);
	    array_1d<double,3>& dn = it->FastGetSolutionStepValue(DISPLACEMENT,1);
	    array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
	    array_1d<double,3>& vmesh = it->FastGetSolutionStepValue(MESH_VELOCITY) = dn1;
	    noalias(vmesh) -= dn;
	    vmesh /= dt;
	    
	    //update the position
	    noalias(it->Coordinates()) = it->GetInitialPosition();
	    noalias(it->Coordinates()) += dn1;   
	}
    }
    
    //**********************************************************************************************
    //**********************************************************************************************
    ///this function erases the elements and conditions which have at least one node marked for erase
    ///@param rModelPart the model part on which we work
    void EraseOuterElements(ModelPart& rModelPart)
    {      
	
	int nerased_el = 0;
        for(ModelPart::ElementsContainerType::iterator it = rModelPart.ElementsBegin(); it!=rModelPart.ElementsEnd(); it++)
	{
	    Geometry< Node<3> >& geom = it->GetGeometry();
	    
//	    bool erase_el = false;
	    for(unsigned int i=0; i<geom.size(); i++)
	    {
		if(geom[i].Is(TO_ERASE))
		{
		    it->Set(TO_ERASE,true);
		    nerased_el++;
		    break;
		}
	    }
	}
	
	if(nerased_el > 0)
	{
	    ModelPart::ElementsContainerType temp_elems_container;
	    temp_elems_container.reserve(rModelPart.Elements().size() - nerased_el);

	    temp_elems_container.swap(rModelPart.Elements());

	    for(ModelPart::ElementsContainerType::iterator it = temp_elems_container.begin() ; it != temp_elems_container.end() ; it++)
	    {
		if( static_cast<bool>(it->Is(TO_ERASE)) == false)
		    (rModelPart.Elements()).push_back(*(it.base()));
	    }
	}
    }
    
    //**********************************************************************************************
    //**********************************************************************************************
    ///this function acts on the inlet nodes so to add new nodes with the objective of reproducing a lagrangian inlet
    ///it works by identyfing the nodes on the inlet and making a copy of them,
    ///the nodes on the inlet are then pushed back to their original position while the
    ///new nodes are generated in the position that was occupied by the lagrangian nodes.
    ///@param rModelPart the model part on which we work
    ///@param rInletNodes is a container which contains the list of inlet nodes
    ///@param rNewNodes is a container that will be filled with the list of the newly created nodes
    ///@param add_nodes_to_model_part identifies if the new nodes need to be included in the model part or no
    void ActOnInlet(ModelPart& rModelPart, ModelPart::NodesContainerType& rInletNodes, ModelPart::NodesContainerType& rNewNodes, bool add_nodes_to_model_part)
    {
//  	unsigned int buffer_size = rModelPart.GetBufferSize();
	
	int aux_id = 1;
	rNewNodes.clear(); // = ModelPart::NodesContainerType();
	rNewNodes.reserve(rInletNodes.size());
	
        //create a copy of all of the nodes in the inlet model part
        for(ModelPart::NodesContainerType::iterator it = rInletNodes.begin(); it!=rInletNodes.end(); it++)
	{
	    //create a new node as a copy of one of the inlet nodes
	    Node<3>::Pointer p_new_node( it->Clone() ); //*(*(it.base()))));

	    //assign an arbitrary Id
	    p_new_node->SetId(aux_id++);
	    
	    rNewNodes.push_back(p_new_node);
	    
// 	    KRATOS_WATCH(*p_new_node);
	    
	    //free velocity and pressure
	    p_new_node->Free(VELOCITY_X);
	    p_new_node->Free(VELOCITY_Y);
	    p_new_node->Free(VELOCITY_Z);
	    p_new_node->Free(PRESSURE);
	    
	    p_new_node->FastGetSolutionStepValue(IS_BOUNDARY) = 0.0;
	    p_new_node->FastGetSolutionStepValue(IS_FREE_SURFACE) = 0.0;
	    p_new_node->FastGetSolutionStepValue(IS_STRUCTURE) = 0.0;
	    p_new_node->FastGetSolutionStepValue(IS_FLUID) = 1.0;
	}
	
	//if required renumber the nodes and add the new nodes to the model_part
	if(add_nodes_to_model_part == true)
	{
	    //renumber the nodes in the model part
	    aux_id = 1;
	    for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
		it->SetId(aux_id++);
	    
	    //assign an Id to the new nodes so that it is consecutive
	    for(ModelPart::NodesContainerType::iterator it = rNewNodes.begin(); it!=rNewNodes.end(); it++)
		it->SetId(aux_id++);
	    
	    //add the nodes to the model part
	    for(ModelPart::NodesContainerType::iterator it = rNewNodes.begin(); it!=rNewNodes.end(); it++)
		rModelPart.AddNode(*(it.base()));
	}
	
	//push back the inlet nodes to their original position
	for(ModelPart::NodesContainerType::iterator it = rInletNodes.begin(); it!=rInletNodes.end(); it++)
	{
	    noalias(it->Coordinates()) = it->GetInitialPosition();
	    noalias(it->FastGetSolutionStepValue(DISPLACEMENT)) = ZeroVector(3);
	    noalias(it->FastGetSolutionStepValue(DISPLACEMENT,1)) = ZeroVector(3);
	    noalias(it->FastGetSolutionStepValue(MESH_VELOCITY)) = it->FastGetSolutionStepValue(VELOCITY);	    
	}
	
    }
    
    //**********************************************************************************************
    //**********************************************************************************************
    void ActOnOutlet(ModelPart& rModelPart, ModelPart::NodesContainerType& rOutletNodes)
    {
// 	for(ModelPart::ElementsContainerType::iterator it = rModelPart.ElementsBegin(); it!=rModelPart.ElementsEnd(); it++)
// 	{
// 	    Geometry<Node<3> >& geom = it->GetGeometry();
// 	    
// 	    //count the nodes to be erased
// 	    int nerased = 0;
// 	    for(unsigned int i=0; i<geom.size(); i++)
// 	      if(geom[i].Is(TO_ERASE)== true)
// 		nerased++;
// 	      
// 	    //fix the pressure if needed (don't if the velocity is fixed)
// 	    if(nerased > 0)
// 	    {
// 		for(unsigned int i=0; i<geom.size(); i++)
// 		{
// 		    if(geom[i].IsFixed(VELOCITY_X) == false)
// 		    {
// 			geom[i].Fix(PRESSURE);
// 			geom[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
// 		    }
// 		}	      
// 	    }
// 	}

	for(ModelPart::NodesContainerType::iterator it = rOutletNodes.begin(); it!=rOutletNodes.end(); it++)
	{
	    noalias(it->Coordinates()) = it->GetInitialPosition();
	    noalias(it->FastGetSolutionStepValue(DISPLACEMENT)) = ZeroVector(3);
	    noalias(it->FastGetSolutionStepValue(DISPLACEMENT,1)) = ZeroVector(3);
	    noalias(it->FastGetSolutionStepValue(MESH_VELOCITY)) = it->FastGetSolutionStepValue(VELOCITY);	    
	}
    }
    
    //**********************************************************************************************
    //**********************************************************************************************
    ///this function marks for erase all of the nodes outside of the limits of the bounding box
    ///@param corner1 first corner of the bounding box
    ///@param corner2 second corner of the bounding box   
    ///@param rNodes list of all nodes in the model
    void MarkOuterNodes(const array_1d<double, 3 > & corner1, 
			const array_1d<double, 3 > & corner2,
                        ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY;

        //add a big number to the id of the nodes to be erased
        int n_erased = 0;
        double xmax, xmin, ymax, ymin, zmax, zmin;



        if (corner1[0] > corner2[0])
        {
            xmax = corner1[0];
            xmin = corner2[0];
        }
        else
        {
            xmax = corner2[0];
            xmin = corner1[0];
        }

        if (corner1[1] > corner2[1])
        {
            ymax = corner1[1];
            ymin = corner2[1];
        }
        else
        {
            ymax = corner2[1];
            ymin = corner1[1];
        }

        if (corner1[2] > corner2[2])
        {
            zmax = corner1[2];
            zmin = corner2[2];
        }
        else
        {
            zmax = corner2[2];
            zmin = corner1[2];
        }



        for (ModelPart::NodesContainerType::iterator in = rNodes.begin(); in != rNodes.end(); in++)
        {
            bool erase = false;
            double& x = in->X();
            double& y = in->Y();
            double& z = in->Z();

            if (x < xmin || x > xmax) erase = true;
            else if (y < ymin || y > ymax) erase = true;
            else if (z < zmin || z > zmax) erase = true;

            if (erase == true)
            {
                n_erased += 1;
                in->Set(TO_ERASE, true);
            }
        }

        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************
    void MarkExcessivelyCloseNodes(ModelPart::NodesContainerType& rNodes, const double admissible_distance_factor)
    {
        KRATOS_TRY;

        const double fact2 = admissible_distance_factor*admissible_distance_factor;



        ModelPart::NodeIterator NodesBegin = rNodes.begin();
        int size = rNodes.size();
        ModelPart::NodeIterator in;
        #pragma omp parallel for firstprivate(size,NodesBegin),private(in)
        for (int k = 0; k < size; k++)
        {
            in = NodesBegin + k;

            if (in->FastGetSolutionStepValue(IS_BOUNDARY) == 0) //if it is not a wall node i can erase
            {
                double hnode2 = in->FastGetSolutionStepValue(NODAL_H);
                hnode2 *= hnode2; //take the square

                //loop on neighbours and erase if they are too close
                for (WeakPointerVector< Node < 3 > >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
                        i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
                {
                    if (static_cast<bool> (i->Is(TO_ERASE)) == false) //we can erase the current node only if the neighb is not to be erased
                    {
                        double dx = i->X() - in->X();
                        double dy = i->Y() - in->Y();
                        double dz = i->Z() - in->Z();

                        double dist2 = dx * dx + dy * dy + dz*dz;

                        if (dist2 < fact2 * hnode2)
                            in->Set(TO_ERASE, true);
                    }
                }
            }
        }
        // 	    }

        KRATOS_CATCH("")
    }    


};

} // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_PFEM2_UTILITIES_INCLUDED  defined


