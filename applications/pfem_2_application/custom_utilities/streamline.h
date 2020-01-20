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

#define pi 3.14159265

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
	    
	    iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY) = iparticle->FastGetSolutionStepValue(VELOCITY);  //AUX_VEL
	    
	    if(iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 1.0) do_move = false;	
	    
	    if( do_move == true  ) //note that we suppose the velocity components to be all fixed
	      {
				array_1d<double,3> old_position = pparticle->Coordinates();
				array_1d<double,3> current_position = pparticle->Coordinates();
				noalias(iparticle->GetInitialPosition()) = old_position;
				iparticle->FastGetSolutionStepValue(DISPLACEMENT,1) = ZeroVector(3);
				subdivisions=20;
				const double small_dt = dt / subdivisions;
	      iparticle->FastGetSolutionStepValue(DISTANCE)=small_dt;
				//
				for (int substep = 0; substep < subdivisions; substep++)
				{
					typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
					Element::Pointer pelement;
					bool is_found = SearchStructure.FindPointOnMesh(current_position, N, pelement, result_begin, max_results);
					iparticle->Set(TO_ERASE, true);
					if (is_found == true)
		      {
					Geometry< Node < 3 > >& geom = pelement->GetGeometry();
					//int nn=0;
					noalias(veulerian) = ZeroVector(3);  
					noalias(acc_particle) = ZeroVector(3);
					for (unsigned int k = 0; k < geom.size(); k++)
					{
					noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY);
					}
					/*
					if(iparticle->FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0) noalias(current_position) += small_dt *0.5*(iparticle->FastGetSolutionStepValue(VELOCITY)+iparticle->FastGetSolutionStepValue(VELOCITY,1));
					else 	
					*/	
					noalias(current_position) += small_dt*veulerian;
					first_time=true;
					pparticle->Set(TO_ERASE, false);
					iparticle->FastGetSolutionStepValue(DISTANCE) += small_dt;
					iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY)=veulerian;
		      }
		    	else
		      {
					double time1=iparticle->FastGetSolutionStepValue(DISTANCE);
					array_1d<double,3> acc;
					acc[0] =  0.0;
					acc[1] = -10.0;
					acc[2] =  0.0;

					noalias(current_position) += small_dt *iparticle->FastGetSolutionStepValue(VELOCITY);
					pparticle->Set(TO_ERASE, false);

		      }//else
		  }//for
				//}//move
				iparticle->FastGetSolutionStepValue(DISPLACEMENT) = current_position - iparticle->GetInitialPosition(); 
	      }//for
	  	}
		//compute mesh velocity
		for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
			{
				array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
				noalias(it->Coordinates()) = it->GetInitialPosition();
				//noalias(it->Coordinates()) += dn1; 
				noalias(it->Coordinates()) += dt * it->FastGetSolutionStepValue(VELOCITY); 
				
			}
    }
      
      void MoveMesh_StreamlinesTn(ModelPart& rModelPart, unsigned int substeps)
      {      
				double deltadt = rModelPart.GetProcessInfo()[DELTA_TIME];
				double dt=0.5 * deltadt; 
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
						
					bool do_move = true;
					bool first_time=false;
					iparticle->FastGetSolutionStepValue(DISTANCE)=0.0;
	    
					iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY) = iparticle->FastGetSolutionStepValue(VELOCITY,1);  //AUX_VEL
							
					if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 1.0) do_move = false;	
	    
					if( do_move == true  ) //note that we suppose the velocity components to be all fixed
						{
							array_1d<double,3> old_position = pparticle->Coordinates();
							array_1d<double,3> current_position = pparticle->Coordinates();
							noalias(iparticle->GetInitialPosition()) = old_position;
							iparticle->FastGetSolutionStepValue(DISPLACEMENT,1) = ZeroVector(3);
							subdivisions=20;
							const double small_dt = dt / subdivisions;
							iparticle->FastGetSolutionStepValue(DISTANCE)=small_dt;
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
									}

									/*
									if(iparticle->FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0) noalias(current_position) += small_dt *iparticle->FastGetSolutionStepValue(VELOCITY);
									else 	
									*/	
									
									noalias(current_position) += small_dt*veulerian;

									pparticle->Set(TO_ERASE, false);
									iparticle->FastGetSolutionStepValue(DISTANCE) += small_dt;
									iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY)=veulerian;
								} //if (is_found == true)
								else
									{
										double time1=iparticle->FastGetSolutionStepValue(DISTANCE);
										array_1d<double,3> acc;
										acc[0] =  0.0;
										acc[1] = -10.0;
										acc[2] =  0.0;

										noalias(current_position) += small_dt *0.5*(iparticle->FastGetSolutionStepValue(VELOCITY,1)+iparticle->FastGetSolutionStepValue(VELOCITY,2));				pparticle->Set(TO_ERASE, false);
		  							}//else
	      					}//for
	  							
							iparticle->FastGetSolutionStepValue(DISPLACEMENT) = current_position - iparticle->GetInitialPosition(); 
      				}//for
					}
    			for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
      				{
							array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
							noalias(it->Coordinates()) = it->GetInitialPosition();
							//noalias(it->Coordinates()) += dn1; 
							noalias(it->Coordinates()) += dt * it->FastGetSolutionStepValue(VELOCITY,1);
      				}
      }

void MoveMesh_RK(ModelPart& rModelPart, unsigned int substeps)
    {      
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        dt *=0.5; 
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();
	
	//do movement
	array_1d<double, 3 > veulerian;
	//double temperature=0.0;
	//array_1d<double, 3 > acc_particle;

        array_1d<double, 3 > v1,v2,v3,v4,vtot,x;

        KRATOS_WATCH("TDim");
        KRATOS_WATCH("TDim");
        KRATOS_WATCH(TDim);
	array_1d<double, TDim + 1 > N;

	const int max_results = 10000;

	typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
	
	const int nparticles = rModelPart.Nodes().size();
	
	#pragma omp parallel for firstprivate(results,N,veulerian,v1,v2,v3,v4,x)

	for (int i = 0; i < nparticles; i++)
	  {


	    array_1d<double,3> current_position;

	    array_1d<double,3> initial_position;		

	    bool is_found=false;	



	   typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;


	    Node < 3 > ::Pointer pparticle = *(iparticle.base());


                initial_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT);

		Element::Pointer pelement;

		//STEP1
		//noalias(current_position) =  initial_position;

		is_found = SearchStructure.FindPointOnMesh(initial_position, N, pelement, result_begin, max_results);

		if(is_found==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
				noalias(v1) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v1) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else noalias(v1) = 0.5 * (pparticle->FastGetSolutionStepValue(VELOCITY)+pparticle->FastGetSolutionStepValue(VELOCITY,1));

		//STEP2
		noalias(x) =  initial_position + 0.5 * dt * v1;
		is_found = SearchStructure.FindPointOnMesh(x, N, pelement, result_begin, max_results);

		if(is_found==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
				noalias(v2) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v2) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else noalias(v2) = 0.5*(pparticle->FastGetSolutionStepValue(VELOCITY)+pparticle->FastGetSolutionStepValue(VELOCITY,1));

		//STEP3
		noalias(x) =  initial_position + 0.5 * dt * v2;
		is_found = SearchStructure.FindPointOnMesh(x, N, pelement, result_begin, max_results);

		if(is_found==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
				noalias(v3) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v3) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else noalias(v3) = 0.5*(pparticle->FastGetSolutionStepValue(VELOCITY)+pparticle->FastGetSolutionStepValue(VELOCITY,1));

		//STEP4

		noalias(x) =  initial_position + dt * v3;
		is_found = SearchStructure.FindPointOnMesh(x, N, pelement, result_begin, max_results);
		if(is_found==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
				noalias(v4) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v4) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else noalias(v4) = 0.5*(pparticle->FastGetSolutionStepValue(VELOCITY)+pparticle->FastGetSolutionStepValue(VELOCITY,1));


                //current_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT);


		if(iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 0.0)
		{
			noalias(x) = initial_position;
	 		noalias(x) += 0.16666666666666666666667*dt*v1;
			noalias(x) += 0.33333333333333333333333*dt*v2;
			noalias(x) += 0.33333333333333333333333*dt*v3;
			noalias(x) += 0.16666666666666666666667*dt*v4;

			pparticle->FastGetSolutionStepValue(DISPLACEMENT) = x - pparticle->GetInitialPosition();
			array_1d<double,3> displacement_aux = x - pparticle->GetInitialPosition();
			pparticle->FastGetSolutionStepValue(DISPLACEMENT) = displacement_aux;
			noalias(pparticle->Coordinates()) = x;

	/*		pparticle->FastGetSolutionStepValue(DISPLACEMENT) += dt * pparticle->FastGetSolutionStepValue(VELOCITY,1)+ dt * dt * (-2.0 * 0.25 + 1.0 ) * 0.5 * pparticle->FastGetSolutionStepValue(ACCELERATION,1) + dt* dt * 0.25 * pparticle->FastGetSolutionStepValue(ACCELERATION);

                        noalias(pparticle->Coordinates()) = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT); 

		
*/



		}//do_move

	}//particles

}


void MoveMesh_FE(ModelPart& rModelPart, unsigned int substeps)
    {      
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        //dt *=0.5; 
	
	//do movement
	array_1d<double, 3 > veulerian;

        array_1d<double, 3 > v1,v2,v3,v4,vtot,x;

	array_1d<double, TDim + 1 > N;

	const int max_results = 10000;

	const int nparticles = rModelPart.Nodes().size();

//	ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

	//ModelPart& model_part=BaseType::GetModelPart();

	
//	#pragma omp parallel for firstprivate(results,N,veulerian,v1,v2,v3,v4,x)



	for (ModelPart::NodeIterator it = rModelPart.NodesBegin();it != rModelPart.NodesEnd(); ++it)
        {

            if(it->FastGetSolutionStepValue(IS_INTERFACE) == 0.0)
		{
		noalias(x) = it->GetInitialPosition() + it->FastGetSolutionStepValue(DISPLACEMENT);
		noalias(x) += dt*it->FastGetSolutionStepValue(VELOCITY);
		it->FastGetSolutionStepValue(DISPLACEMENT) = x - it->GetInitialPosition();
		noalias(it->Coordinates()) = x;
	        }
	}

}


	void Force( ModelPart& model_part)
	{
	  
        KRATOS_TRY

		//ModelPart& model_part=BaseType::GetModelPart();
		double time = model_part.GetProcessInfo()[TIME];

		const double dt = model_part.GetProcessInfo()[DELTA_TIME];

		double time_dtn = time-dt;

		//T^{n+1}
		for(ModelPart::NodesContainerType::iterator in = model_part.NodesBegin() ; in != model_part.NodesEnd() ; ++in)
		{
		double nu=0.001;

		double fx=(in->X()*in->X()) * ( (1.0-in->X())*(1.0-in->X())); 

		double dfx= 2.0*in->X()*( (1.0-in->X()) *(1.0-in->X())) - 2.0*(in->X() * in->X() )*(1.0-in->X()); 

		double d2fx = 2.0*( (1.0-in->X())*(1.0-in->X()) ) - 8.0*in->X()*(1.0-in->X()) + 2.0*(in->X() * in->X()); 

		double fy = (in->Y()*in->Y()) * ( (1.0-in->Y())* (1.0-in->Y()) ); 

		double dfy = 2.0*in->Y()*( (1.0-in->Y())*(1.0-in->Y()) ) - 2.0*(in->Y()*in->Y())*(1.0-in->Y()); 

		double d2fy = 2.0*( (1.0-in->Y())*(1.0-in->Y()) ) - 8.0*in->Y()*(1.0-in->Y()) + 2.0*(in->Y()*in->Y());

		double d3fx = -12.0*(1.0-in->X()) + 12.0*in->X(); 
 
		double d3fy =-12.0*(1.0-in->Y()) + 12.0*in->Y(); 

		double ht = cos( pi * time * 4.0) * exp(-time);

		double dht = - 4.0 *pi * sin( pi * time * 4.0) * exp(-time) - cos( pi * time * 4.0) * exp(-time); 

		double ux =  100.0 * ht * fx * dfy; 
		double uy = -100.0 * ht * dfx * fy; 

		double duxdx = 100.0 * ht * dfx * dfy;
		double duxdy = 100.0 * ht * fx * d2fy; 
		double duydx = -100.0 * ht * d2fx * fy; 
		double duydy = -100.0 * ht * dfx * dfy;

		double d2uxdx2 = 100.0 * ht * ( d2fx * dfy );
		double d2uxdy2 = 100.0 * ht * ( fx * d3fy );
		double d2uydx2 = -100.0 * ht * ( d3fx * fy );
		double d2uydy2 = -100.0 * ht * ( dfx * d2fy ); ///

		double d2uxdxdy = 100.0 * ht * dfx * d2fy;
                double d2uydxdy = -100.0 * ht * d2fx * dfy;

		double dpx = 200.0*in->X();
		double dpy = 0.0;
 		//1. dynamic term
		double Fx = 100.0*dht*fx*dfy;
		double Fy = -100.0*dht*fy*dfx;
		// 2. convective term
    		//SHOULD BE ZERO FOR LAGRANGIAN MODELS
		Fx += ux * duxdx + uy * duxdy;
		Fy += ux * duydx + uy * duydy;
		//3. viscous term
		Fx -= nu * 2.0 * ( d2uxdx2 + 0.5*(d2uydxdy + d2uxdy2));
		Fy -= nu * 2.0 * ( d2uydy2 + 0.5*(d2uxdxdy + d2uydx2));
    
		//4. pressure gradient
		Fx += dpx;
		Fy += dpy;

		in->FastGetSolutionStepValue(BODY_FORCE_X)=Fx;
		in->FastGetSolutionStepValue(BODY_FORCE_Y)=Fy;
		}

		//T^{n}
/*		for(ModelPart::NodesContainerType::iterator in = model_part.NodesBegin() ; in != model_part.NodesEnd() ; ++in)
		{
		double nu=0.001;

		double fx=(in->X()*in->X()) * ( (1.0-in->X())*(1.0-in->X())); 

		double dfx= 2.0*in->X()*( (1.0-in->X()) *(1.0-in->X())) - 2.0*(in->X() * in->X() )*(1.0-in->X()); 

		double d2fx = 2.0*( (1.0-in->X())*(1.0-in->X()) ) - 8.0*in->X()*(1.0-in->X()) + 2.0*(in->X() * in->X()); 

		double fy = (in->Y()*in->Y()) * ( (1.0-in->Y())* (1.0-in->Y()) ); 

		double dfy = 2.0*in->Y()*( (1.0-in->Y())*(1.0-in->Y()) ) - 2.0*(in->Y()*in->Y())*(1.0-in->Y()); 

		double d2fy = 2.0*( (1.0-in->Y())*(1.0-in->Y()) ) - 8.0*in->Y()*(1.0-in->Y()) + 2.0*(in->Y()*in->Y());

		double d3fx = -12.0*(1.0-in->X()) + 12.0*in->X(); 
 
		double d3fy =-12.0*(1.0-in->Y()) + 12.0*in->Y(); 

		double ht = cos( pi * time_dtn * 4.0) * exp(-time_dtn);

		double dht = - 4.0 *pi * sin( pi * time_dtn * 4.0) * exp(-time_dtn) - cos( pi * time_dtn * 4.0) * exp(-time_dtn); 

		double ux =  100.0 * ht * fx * dfy; 
		double uy = -100.0 * ht * dfx * fy; 

		double duxdx = 100.0 * ht * dfx * dfy;
		double duxdy = 100.0 * ht * fx * d2fy; 
		double duydx = -100.0 * ht * d2fx * fy; 
		double duydy = -100.0 * ht * dfx * dfy;

		double d2uxdx2 = 100.0 * ht * ( d2fx * dfy );
		double d2uxdy2 = 100.0 * ht * ( fx * d3fy );
		double d2uydx2 = -100.0 * ht * ( d3fx * fy );
		double d2uydy2 = -100.0 * ht * ( dfx * d2fy ); ///

		double d2uxdxdy = 100.0 * ht * dfx * d2fy;
                double d2uydxdy = -100.0 * ht * d2fx * dfy;

		double dpx = 200.0*in->X();
		double dpy = 0.0;
 		//1. dynamic term
		double Fx = 100.0*dht*fx*dfy;
		double Fy = -100.0*dht*fy*dfx;
		// 2. convective term
    		//SHOULD BE ZERO FOR LAGRANGIAN MODELS
		Fx += ux * duxdx + uy * duxdy;
		Fy += ux * duydx + uy * duydy;
		//3. viscous term
		Fx -= nu * 2.0 * ( d2uxdx2 + 0.5*(d2uydxdy + d2uxdy2));
		Fy -= nu * 2.0 * ( d2uydy2 + 0.5*(d2uxdxdy + d2uydx2));
    
		//4. pressure gradient
		Fx += dpx;
		Fy += dpy;

		in->FastGetSolutionStepValue(BODY_FORCE_X)+=Fx;
		in->FastGetSolutionStepValue(BODY_FORCE_Y)+=Fy;
		}

		//*0.5
		for(ModelPart::NodesContainerType::iterator in = model_part.NodesBegin() ; in != model_part.NodesEnd() ; ++in)
		{
		in->FastGetSolutionStepValue(BODY_FORCE_X)*=0.5;
		in->FastGetSolutionStepValue(BODY_FORCE_Y)*=0.5;
		}

*/

	KRATOS_CATCH("")
	}


/*
void MoveMesh_RK(ModelPart& rModelPart, unsigned int substeps)
    {      
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        dt *=0.5; 
	//KRATOS_ERROR(std::logic_error, "element with zero vol found", "");
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();
	
	//do movement
	array_1d<double, 3 > veulerian;
	//double temperature=0.0;
	//array_1d<double, 3 > acc_particle;

	array_1d<double, 3 > coefficients_RK;

	array_1d<double, 2 > coefficients_RK_aux;

        coefficients_RK[0]=0.5; coefficients_RK[1]=0.5; coefficients_RK[2]=1.0; coefficients_RK[3]=0.166666667;

	coefficients_RK_aux[0]=2.0; coefficients_RK_aux[1]=2.0; coefficients_RK_aux[2]=1.0; 

	array_1d<double, TDim + 1 > N;

	const int max_results = 10000;

	typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
	
	const int nparticles = rModelPart.Nodes().size();
	
	#pragma omp parallel for firstprivate(results,N,veulerian)

	for (int i = 0; i < nparticles; i++)
	  {
	    int subdivisions;

	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;

	    Node < 3 > ::Pointer pparticle = *(iparticle.base());
	    
	    bool do_move = true;

	    bool first_time=false;

	    iparticle->FastGetSolutionStepValue(IS_POROUS)=false;	

	    iparticle->FastGetSolutionStepValue(MESH_VELOCITY)=iparticle->FastGetSolutionStepValue(VELOCITY); 
	    
	    iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY) = iparticle->FastGetSolutionStepValue(VELOCITY);  //AUX_VEL
	    
	    if(iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 1.0) do_move = false;	


	    array_1d<double,3> current_position;

	    array_1d<double,3> initial_position;

	    bool not_found=false;	
	    if( do_move == true  ) //note that we suppose the velocity components to be all fixed
	      {

                //array_1d<double,3> current_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT);

	        initial_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT);


		subdivisions=3;

	for (int substep = 0; substep < subdivisions; substep++)
		{


			noalias(current_position) =  initial_position + coefficients_RK[substep] * dt * iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY);

			typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
		

			Element::Pointer pelement;

			bool is_found = SearchStructure.FindPointOnMesh(current_position, N, pelement, result_begin, max_results);

			//iparticle->Set(TO_ERASE, true);

			if (is_found == true)
		      		{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();

				noalias(veulerian) = ZeroVector(3);  

				for (unsigned int k = 0; k < geom.size(); k++)
					{
						noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY);
					}

					iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY)=veulerian; //save velocity

					first_time=true;


					iparticle->FastGetSolutionStepValue(MESH_VELOCITY) += coefficients_RK_aux[substep]*iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY);


					//pparticle->Set(TO_ERASE, false);
		      		}


		    	else
		      		{
					iparticle->FastGetSolutionStepValue(IS_POROUS)=true;
					not_found=true;

		      		}//else
		  }//for



		//iparticle->FastGetSolutionStepValue(DISPLACEMENT) = current_position - iparticle->GetInitialPosition(); 


	      }//for
	  	}
		//compute mesh velocity

		for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
			{


                               if(it->FastGetSolutionStepValue(IS_INTERFACE)==false) {






				
				if(it->FastGetSolutionStepValue(IS_POROUS)==true)
					{


						array_1d<double,3> current_position = it->GetInitialPosition() + it->FastGetSolutionStepValue(DISPLACEMENT);

						noalias(current_position) += dt * it->FastGetSolutionStepValue(VELOCITY);

		 				it->FastGetSolutionStepValue(DISPLACEMENT) = current_position - it->GetInitialPosition();

						noalias(it->Coordinates()) = current_position;


					}
				else
					{
						array_1d<double,3> current_position = it->GetInitialPosition() + it->FastGetSolutionStepValue(DISPLACEMENT);

						noalias(current_position) += 0.166666667 * dt * it->FastGetSolutionStepValue(MESH_VELOCITY);

		 				it->FastGetSolutionStepValue(DISPLACEMENT) = current_position - it->GetInitialPosition();

						noalias(it->Coordinates()) = current_position;
					}
					

				}
			}







    }*/
      
        /*for (ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); it++)
        {
            array_1d<double,3> delta_disp = it->FastGetSolutionStepValue(DISPLACEMENT);
            noalias(delta_disp) -= it->FastGetSolutionStepValue(DISPLACEMENT,1);

            double norm_delta_disp = norm_2(delta_disp);

            array_1d<double,3> v_old = it->FastGetSolutionStepValue(VELOCITY,1);
            double norm_v = norm_2(v_old);

            if(norm_delta_disp*3.0 < norm_v*dt*2.0 )
                it->Set(TO_ERASE, true);
            //if(norm_delta_disp* (0.333333333333333*0.001) >  norm_v*dt*2.0 )
            //    it->Set(TO_ERASE, true);
        }*/



void MoveMesh_RKTn(ModelPart& rModelPart, unsigned int substeps)
    {      
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        dt *=0.5; 
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();
	
	//do movement
	array_1d<double, 3 > veulerian;
	//double temperature=0.0;
	//array_1d<double, 3 > acc_particle;

        array_1d<double, 3 > v1,v2,v3,v4,vtot,x;


	array_1d<double, TDim + 1 > N;

	const int max_results = 10000;

	typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
	
	const int nparticles = rModelPart.Nodes().size();
	
	#pragma omp parallel for firstprivate(results,N,veulerian,v1,v2,v3,v4,x)

	for (int i = 0; i < nparticles; i++)
	  {


	    array_1d<double,3> current_position;

	    array_1d<double,3> initial_position;		

	    bool is_found=false;	



	   typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;


	    Node < 3 > ::Pointer pparticle = *(iparticle.base());



  	        //ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;

		


                initial_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT,1);

		Element::Pointer pelement;

		//STEP1
		//noalias(current_position) =  initial_position;

		is_found = SearchStructure.FindPointOnMesh(initial_position, N, pelement, result_begin, max_results);

		if(is_found==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
				noalias(v1) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v1) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else noalias(v1) = pparticle->FastGetSolutionStepValue(VELOCITY);

		//STEP2
		noalias(x) =  initial_position + 0.5 * dt * v1;
		is_found = SearchStructure.FindPointOnMesh(x, N, pelement, result_begin, max_results);

		if(is_found==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
				noalias(v2) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v2) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else noalias(v2) = pparticle->FastGetSolutionStepValue(VELOCITY);

		//STEP3
		noalias(x) =  initial_position + 0.5 * dt * v2;
		is_found = SearchStructure.FindPointOnMesh(x, N, pelement, result_begin, max_results);

		if(is_found==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
				noalias(v3) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v3) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else noalias(v3) = pparticle->FastGetSolutionStepValue(VELOCITY);

		//STEP4

		noalias(x) =  initial_position + dt * v3;
		is_found = SearchStructure.FindPointOnMesh(x, N, pelement, result_begin, max_results);
		if(is_found==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
				noalias(v4) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v4) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else noalias(v4) = pparticle->FastGetSolutionStepValue(VELOCITY);


                //current_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT);


		if(iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 0.0)
		{
			noalias(x) = initial_position;
	 		/*noalias(x) += 0.16666666666666666666667*dt*v1;
			noalias(x) += 0.33333333333333333333333*dt*v2;
			noalias(x) += 0.33333333333333333333333*dt*v3;
			noalias(x) += 0.16666666666666666666667*dt*v4;*/

                        noalias(x) += dt * iparticle->FastGetSolutionStepValue(VELOCITY,1);

			pparticle->FastGetSolutionStepValue(DISPLACEMENT) = x - pparticle->GetInitialPosition();

			noalias(pparticle->Coordinates()) = x;


/*		
			pparticle->FastGetSolutionStepValue(DISPLACEMENT) = pparticle->FastGetSolutionStepValue(DISPLACEMENT,1) + dt * pparticle->FastGetSolutionStepValue(VELOCITY,1)+ dt * dt * (-2.0 * 0.25 + 1.0 ) * 0.5 * pparticle->FastGetSolutionStepValue(ACCELERATION,1) + dt* dt * 0.25 * pparticle->FastGetSolutionStepValue(ACCELERATION);

                        noalias(pparticle->Coordinates()) = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT); 

*/

		}//do_move

	}//particles

}


/*void MoveMesh_RKTn(ModelPart& rModelPart, unsigned int substeps)
    {      
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        dt *=0.5; 
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();
	
	//do movement
	array_1d<double, 3 > veulerian;
	//double temperature=0.0;
	//array_1d<double, 3 > acc_particle;

	array_1d<double, 3 > coefficients_RK;

	array_1d<double, 2 > coefficients_RK_aux;

        coefficients_RK[0]=0.5; coefficients_RK[1]=0.5; coefficients_RK[2]=1.0; coefficients_RK[3]=0.166666667;

	coefficients_RK_aux[0]=2.0; coefficients_RK_aux[1]=2.0; coefficients_RK_aux[2]=1.0; 

	array_1d<double, TDim + 1 > N;

	const int max_results = 10000;

	typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
	
	const int nparticles = rModelPart.Nodes().size();
	
	#pragma omp parallel for firstprivate(results,N,veulerian)

	for (int i = 0; i < nparticles; i++)
	  {
	    int subdivisions;

	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;

	    Node < 3 > ::Pointer pparticle = *(iparticle.base());
	    
	    bool do_move = true;

	    bool first_time=false;

	    iparticle->FastGetSolutionStepValue(IS_POROUS)=false;	

	    iparticle->FastGetSolutionStepValue(MESH_VELOCITY)=iparticle->FastGetSolutionStepValue(VELOCITY,1); 
	    
	    iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY) = iparticle->FastGetSolutionStepValue(VELOCITY,1);  //AUX_VEL
	    
	    if(iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 1.0) do_move = false;

	    array_1d<double,3> current_position;

	    array_1d<double,3> initial_position;		

	    bool not_found=false;	
	    if( do_move == true  ) //note that we suppose the velocity components to be all fixed
	      {

                initial_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT,1);



		subdivisions=3;

	for (int substep = 0; substep < subdivisions; substep++)
		{


			noalias(current_position) =  initial_position + coefficients_RK[substep] * dt * iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY);

			typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
		

			Element::Pointer pelement;

			bool is_found = SearchStructure.FindPointOnMesh(current_position, N, pelement, result_begin, max_results);


			if (is_found == true)
		      		{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();

				noalias(veulerian) = ZeroVector(3);  

				for (unsigned int k = 0; k < geom.size(); k++)
					{
						noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY,1);
					}

					iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY)=veulerian; 

					first_time=true;


					iparticle->FastGetSolutionStepValue(MESH_VELOCITY) += coefficients_RK_aux[substep]*iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY);


					
		      		}


		    	else
		      		{
					iparticle->FastGetSolutionStepValue(IS_POROUS)=true;
					not_found=true;

		      		}//else
		  }//for

	      }//for
	  	}
		//compute mesh velocity

		for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
			{


                               if(it->FastGetSolutionStepValue(IS_INTERFACE)==false) {


				
				if(it->FastGetSolutionStepValue(IS_POROUS)==true)
				{


				array_1d<double,3> current_position = it->GetInitialPosition() + it->FastGetSolutionStepValue(DISPLACEMENT,1);

				noalias(current_position) += dt * it->FastGetSolutionStepValue(VELOCITY,1);

 				it->FastGetSolutionStepValue(DISPLACEMENT) = current_position - it->GetInitialPosition();

				noalias(it->Coordinates()) = current_position;


				}
				else
				{
				array_1d<double,3> current_position = it->GetInitialPosition() + it->FastGetSolutionStepValue(DISPLACEMENT,1);

				noalias(current_position) += 0.166666667 * dt * it->FastGetSolutionStepValue(MESH_VELOCITY);

 				it->FastGetSolutionStepValue(DISPLACEMENT) = current_position - it->GetInitialPosition();

				noalias(it->Coordinates()) = current_position;
				}
			}
		}
    }

*/

/*
void MoveMesh_RKTn(ModelPart& rModelPart, unsigned int substeps)
    {      
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];

        dt *=0.5; 

	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);

	SearchStructure.UpdateSearchDatabase();
	
	//do movement
	array_1d<double, 3 > veulerian;

	array_1d<double, 3 > coefficients_RK;

	array_1d<double, 2 > coefficients_RK_aux;

        coefficients_RK[0]=0.5; coefficients_RK[1]=0.5; coefficients_RK[2]=1.0; coefficients_RK[3]=0.166666667;

	coefficients_RK_aux[0]=2.0; coefficients_RK_aux[1]=2.0; coefficients_RK_aux[2]=1.0; 

	array_1d<double, TDim + 1 > N;

	const int max_results = 10000;

	typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
	
	const int nparticles = rModelPart.Nodes().size();
	
	#pragma omp parallel for firstprivate(results,N,veulerian)

	for (int i = 0; i < nparticles; i++)
	  {
	    int subdivisions;

	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;

	    Node < 3 > ::Pointer pparticle = *(iparticle.base());
	    
	    bool do_move = true;

	    bool first_time=false;

	    bool not_found=false;

	    iparticle->FastGetSolutionStepValue(IS_POROUS)=false;	

	    iparticle->FastGetSolutionStepValue(MESH_VELOCITY)=iparticle->FastGetSolutionStepValue(VELOCITY,1); 
	    
	    iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY) = iparticle->FastGetSolutionStepValue(VELOCITY,1);  //AUX_VEL
	    
	    if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 1.0) do_move = false;	
	    
	    if( do_move == true  ) //note that we suppose the velocity components to be all fixed
	      {

		array_1d<double,3> current_position = pparticle->Coordinates();


		iparticle->FastGetSolutionStepValue(DISPLACEMENT,1) = ZeroVector(3);

		subdivisions=3;


		for (int substep = 0; substep < subdivisions; substep++)
		{


			noalias(current_position) =  pparticle->Coordinates() + coefficients_RK[substep] * dt * iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY);

			typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

			Element::Pointer pelement;

			bool is_found = SearchStructure.FindPointOnMesh(current_position, N, pelement, result_begin, max_results);

			iparticle->Set(TO_ERASE, true);

			if (is_found == true)
		      		{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();

				noalias(veulerian) = ZeroVector(3);  

				for (unsigned int k = 0; k < geom.size(); k++)
					{
						noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY,1);
					}

					iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY)=veulerian; //save velocity

					first_time=true;


					iparticle->FastGetSolutionStepValue(MESH_VELOCITY) += coefficients_RK_aux[substep]*iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY);


					pparticle->Set(TO_ERASE, false);
		      		}
		    	else
		      		{
					
					not_found=true;

					iparticle->FastGetSolutionStepValue(IS_POROUS)=true;


		      		}//else
		  }//for



		//iparticle->FastGetSolutionStepValue(DISPLACEMENT) = current_position - iparticle->GetInitialPosition(); 


	      }//for
	  	}
		//compute mesh velocity
		for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
			{

				if(it->FastGetSolutionStepValue(IS_POROUS)==true)  noalias(it->Coordinates()) +=  dt * it->FastGetSolutionStepValue(VELOCITY);
				else
				noalias(it->Coordinates()) += 0.166666667 * dt * it->FastGetSolutionStepValue(MESH_VELOCITY);
				
			}
    }
      
*/

      double Calculate_Vol(ModelPart & rLagrangianModelPart)
      {
        KRATOS_TRY


	  //particles
	  for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin(); node_it != rLagrangianModelPart.NodesEnd(); node_it++)
	    {
	      //if( node_it->GetValue(NEIGHBOUR_ELEMENTS).size() != 0) (node_it)->FastGetSolutionStepValue(K0) = 0.0;
	      (node_it)->FastGetSolutionStepValue(K0) = 0.0;	
	    }

	for (ModelPart::ElementsContainerType::iterator el_it = rLagrangianModelPart.ElementsBegin();el_it != rLagrangianModelPart.ElementsEnd(); el_it++)
	  {
            Geometry<Node < 3 > >& geom = el_it->GetGeometry();
	    int nd=2; 
	    double area=0.0;
	    double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3; 	
            x0 = geom[0].X();
            y0 = geom[0].Y();
	    if(nd==3) z0 = geom[0].Z();
            x1 = geom[1].X();
            y1 = geom[1].Y();
            if(nd==3) z1 = geom[1].Z();
	    x2 = geom[2].X();
            y2 = geom[2].Y();
	    if(nd==3) z2 = geom[2].Z();
	    if(nd==3)
	    {
	    x3 = geom[3].X();
            y3 = geom[3].Y();
	    z3 = geom[3].Z();
	    }
	    if(nd==2)	
	    area = CalculateVol(x0, y0, x1, y1, x2, y2);
	    else
	    area=CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

	    geom[0].FastGetSolutionStepValue(K0) += area * 0.3333333333;
            geom[1].FastGetSolutionStepValue(K0) += area * 0.3333333333;
            geom[2].FastGetSolutionStepValue(K0) += area * 0.3333333333;
	    if(nd==3)	
	    geom[3].FastGetSolutionStepValue(K0) += area * 0.3333333333;
	    	
	  }

	double sum=0.0;

	for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin(); node_it != rLagrangianModelPart.NodesEnd(); node_it++)
	  {

	    sum +=(node_it)->FastGetSolutionStepValue(K0) ;
	  }

	return sum;

        KRATOS_CATCH("")
	  }

      inline double CalculateVol(const double x0, const double y0,
				 const double x1, const double y1,
				 const double x2, const double y2
				 )
      {
	return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
      }


      inline double CalculateVol(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2,const double x3, const double y3, const double z3 )
      {
	double x10 = x1 - x0;
	double y10 = y1 - y0;
	double z10 = z1 - z0;

	double x20 = x2 - x0;
	double y20 = y2 - y0;
	double z20 = z2 - z0;

	double x30 = x3 - x0;
	double y30 = y3 - y0;
	double z30 = z3 - z0;

	double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
	return detJ * 0.1666666666666666666667;
      }
      
    private:
      
};

} // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED  defined


	
