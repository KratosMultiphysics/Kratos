//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti.
//

#if !defined(KRATOS_PARTICLESM_UTILITIES_INCLUDED )
#define  KRATOS_PARTICLESM_UTILITIES_INCLUDED

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
#include "includes/kratos_flags.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "pfem_2_application.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"
//#include "utilities/enrichment_utilities.h"

#include <boost/timer.hpp>
#include "utilities/timer.h"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace Kratos
{
  /*
  template< class T, std::size_t dim >
    class DistanceCalculator1
    {
    public:
      
      double operator()(T const& p1, T const& p2)
      {
	double dist = 0.0;
        for (std::size_t i = 0; i < dim; i++)
	  {
            double tmp = p1[i] - p2[i];
            dist += tmp*tmp;
	  }
        return dist; //square distance because it is easier to work without the square root//
      }
      
    };
  */
  template<std::size_t TDim> class Streamline
    {
    public:
      KRATOS_CLASS_POINTER_DEFINITION(Streamline<TDim>);
      
      
    void MoveMesh_Streamlines(ModelPart& rModelPart, unsigned int substeps)
      {      
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        dt *=0.5; 
	//KRATOS_ERROR(std::logic_error, "element with zero vol found", "");
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();
	
	//do movement
	array_1d<double, 3 > veulerian;
	//double temperature=0.0;
	array_1d<double, 3 > acc_particle;
	array_1d<double, TDim + 1 > N;
	const int max_results = 10000;
	typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
	
	const int nparticles = rModelPart.Nodes().size();
	
#pragma omp parallel for firstprivate(results,N,veulerian,acc_particle)
	for (int i = 0; i < nparticles; i++)
	  {
	    //int substep = 0;
	    int subdivisions = 1;
	    //double temperature=0.0;
	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
	    Node < 3 > ::Pointer pparticle = *(iparticle.base());
	    //small_dt = dt / subdivisions;
	    
	    bool do_move = true;
	    bool first_time=false;
	    iparticle->FastGetSolutionStepValue(DISTANCE)=0.0;
	    
	    //iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY) = iparticle->FastGetSolutionStepValue(VELOCITY,1);  //AUX_VEL
	    
	    if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 1.0) do_move = false;	
	    
	    if( do_move == true  ) //note that we suppose the velocity components to be all fixed
	      {
		array_1d<double,3> old_position = pparticle->Coordinates();
		array_1d<double,3> current_position = pparticle->Coordinates();
		noalias(iparticle->GetInitialPosition()) = old_position;
		iparticle->FastGetSolutionStepValue(DISPLACEMENT,1) = ZeroVector(3);
		subdivisions=10;
		const double small_dt = dt / subdivisions;
		//
		for (int substep = 0; substep < subdivisions; substep++)
		  {
		    typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
		    Element::Pointer pelement;
		    bool is_found = SearchStructure.FindPointOnMesh(current_position, N, pelement, result_begin, max_results);
		    iparticle->Set(TO_ERASE, true);
		    //KRATOS_WATCH(is_found);
		    if (is_found == true)
		      {
			Geometry< Node < 3 > >& geom = pelement->GetGeometry();
			//int nn=0;
			noalias(veulerian) = ZeroVector(3);  
			noalias(acc_particle) = ZeroVector(3);

			for (unsigned int k = 0; k < geom.size(); k++)
			  {
			    noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY);
			    noalias(acc_particle) += N[k] * geom[k].FastGetSolutionStepValue(ANGULAR_ACCELERATION);
			  }
			
	
			//noalias(current_position) += small_dt*veulerian + 0.5 * small_dt * small_dt * acc_particle;
 			noalias(current_position) += small_dt*veulerian + 0.5 * dt * small_dt * acc_particle;
			pparticle->Set(TO_ERASE, false);

			iparticle->FastGetSolutionStepValue(DISTANCE) += small_dt;
		      }
		    else
		      {
			double time1=iparticle->FastGetSolutionStepValue(DISTANCE);
			array_1d<double,3> acc;
			acc[0] =  0.0;
			acc[1] = -10.0;
			acc[2] =  0.0;
			
		      }
		  }//for
		
		//update the displacement BUT DO NOT OVERWRITE THE POSITION!!                
		iparticle->FastGetSolutionStepValue(DISPLACEMENT) = current_position - iparticle->GetInitialPosition(); 
		//KRATOS_WATCH(iparticle->FastGetSolutionStepValue(DISPLACEMENT));
	      }//move
	  }
	//compute mesh velocity
	for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
	  {
	    //array_1d<double,3>& dn = it->FastGetSolutionStepValue(DISPLACEMENT,1);
	    array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
	    noalias(it->Coordinates()) = it->GetInitialPosition();
	    noalias(it->Coordinates()) += dn1; 
	    //noalias(it->Coordinates()) += 0.5 * dt*(it->FastGetSolutionStepValue(VELOCITY)+it->FastGetSolutionStepValue(VELOCITY,1)); 	

	  }
      }
      
void MoveMesh_StreamlinesTn(ModelPart& rModelPart, unsigned int substeps)
      {      
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        dt *=0.5; 
	//KRATOS_ERROR(std::logic_error, "element with zero vol found", "");
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();
	
	//do movement
	array_1d<double, 3 > veulerian;
	//double temperature=0.0;
	array_1d<double, 3 > acc_particle;
	array_1d<double, TDim + 1 > N;
	const int max_results = 10000;
	typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
	
	const int nparticles = rModelPart.Nodes().size();
	
#pragma omp parallel for firstprivate(results,N,veulerian,acc_particle)
	for (int i = 0; i < nparticles; i++)
	  {
	    //int substep = 0;
	    int subdivisions = 1;
	    //double temperature=0.0;
	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
	    Node < 3 > ::Pointer pparticle = *(iparticle.base());
	    //small_dt = dt / subdivisions;
	    
	    bool do_move = true;
	    bool first_time=false;
	    iparticle->FastGetSolutionStepValue(DISTANCE)=0.0;
	    
	    //iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY) = iparticle->FastGetSolutionStepValue(VELOCITY,1);  //AUX_VEL
	    
	    if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 1.0) do_move = false;	
	    
	    if( do_move == true  ) //note that we suppose the velocity components to be all fixed
	      {
		array_1d<double,3> old_position = pparticle->Coordinates();
		array_1d<double,3> current_position = pparticle->Coordinates();
		noalias(iparticle->GetInitialPosition()) = old_position;
		iparticle->FastGetSolutionStepValue(DISPLACEMENT,1) = ZeroVector(3);
		subdivisions=10;
		const double small_dt = dt / subdivisions;
		//
		for (int substep = 0; substep < subdivisions; substep++)
		  {
		    typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
		    Element::Pointer pelement;
		    bool is_found = SearchStructure.FindPointOnMesh(current_position, N, pelement, result_begin, max_results);
		    iparticle->Set(TO_ERASE, true);
		    //KRATOS_WATCH(is_found);
		    if (is_found == true)
		      {
			Geometry< Node < 3 > >& geom = pelement->GetGeometry();
			//int nn=0;
			noalias(veulerian) = ZeroVector(3);  
			noalias(acc_particle) = ZeroVector(3);
			for (unsigned int k = 0; k < geom.size(); k++)
			  {
			    noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY,1);
			    noalias(acc_particle) += N[k] * geom[k].FastGetSolutionStepValue(ANGULAR_ACCELERATION,1);

			  }
			
//			noalias(current_position) += small_dt*veulerian + 0.5 * small_dt * small_dt * acc_particle;
			noalias(current_position) += small_dt*veulerian + 0.5 * dt * small_dt * acc_particle;
			pparticle->Set(TO_ERASE, false);
			iparticle->FastGetSolutionStepValue(DISTANCE) += small_dt;
		      }
		    else
		      {
			double time1=iparticle->FastGetSolutionStepValue(DISTANCE);
			array_1d<double,3> acc;
			acc[0] =  0.0;
			acc[1] = -10.0;
			acc[2] =  0.0;
			
		      }
		  }//for
		
		//update the displacement BUT DO NOT OVERWRITE THE POSITION!!                
		iparticle->FastGetSolutionStepValue(DISPLACEMENT) = current_position - iparticle->GetInitialPosition(); 
		//KRATOS_WATCH(iparticle->FastGetSolutionStepValue(DISPLACEMENT));
	      }//move
	  }
	//compute mesh velocity
	for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
	  {
	    //array_1d<double,3>& dn = it->FastGetSolutionStepValue(DISPLACEMENT,1);
	    array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
	    noalias(it->Coordinates()) = it->GetInitialPosition();
	    noalias(it->Coordinates()) += dn1; 
	    //noalias(it->Coordinates()) += 0.5 * dt*(it->FastGetSolutionStepValue(VELOCITY)+it->FastGetSolutionStepValue(VELOCITY,1)); 	

	  }
      }
      
    private:
      
 
      
    };

} // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED  defined


	
