// KRATOS
// _____   __               __  __      _ _   _
//|  __ \ / _|             |  \/  |    | | | (_)
//| |__) | |_ ___ _ __ ___ | \  / | ___| | |_ _ _ __   __ _
//|  ___/|  _/ _ \ '_ ` _ \| |\/| |/ _ \ | __| | '_ \ / _` |
//| |    | ||  __/ | | | | | |  | |  __/ | |_| | | | | (_| |
//|_|    |_| \___|_| |_| |_|_|  |_|\___|_|\__|_|_| |_|\__, |
//                                                     __/ |
//                                                    |___/ APPLICATION
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//


#if !defined(KRATOS_PARTICLES2M_UTILITIES_INCLUDED )
#define  KRATOS_PARTICLES2M_UTILITIES_INCLUDED

//#define PRESSURE_ON_EULERIAN_MESH
//#define USE_FEW_PARTICLES

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
#include "pfem_melting_application.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"
#include "includes/deprecated_variables.h"
#include "utilities/variable_utils.h"
//#include "utilities/enrichment_utilities.h"
//#
#include <boost/timer.hpp>
#include "utilities/timer.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define pi 3.14159265

namespace Kratos
{

  template<std::size_t TDim> class Streamline
    {
    public:
      KRATOS_CLASS_POINTER_DEFINITION(Streamline<TDim>);



      void SubSteppingElementbasedSI(ModelPart& rModelPart, unsigned int substeps)
      {
				double deltadt = rModelPart.GetProcessInfo()[DELTA_TIME];
				double dt=1.0 * deltadt;

				BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
				SearchStructure.UpdateSearchDatabase();

				//do movement
				array_1d<double, 3 > veulerian;
				//double temperature=0.0;
				array_1d<double, 3 > acc_particle;
				///array_1d<double, TDim + 1 > N;
				Vector N; N.resize(TDim + 1);
				const int max_results = 10000;
				typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

				const int nparticles = rModelPart.Nodes().size();

				#pragma omp parallel for firstprivate(results,N,veulerian,acc_particle)
				for (int i = 0; i < nparticles; i++)
					{
					//int substep = 0;
					int subdivisions = 10;
					//double temperature=0.0;
					ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
					Node < 3 > ::Pointer pparticle = *(iparticle.base());

					bool do_move = true;
					//bool first_time=false;
					//iparticle->FastGetSolutionStepValue(DISTANCE)=0.0;

					//iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY) = iparticle->FastGetSolutionStepValue(VELOCITY);  //AUX_VEL

					if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 1.0) do_move = false;

					if( do_move == true  ) //note that we suppose the velocity components to be all fixed
						{
							array_1d<double,3> old_position = pparticle->Coordinates();
							array_1d<double,3> current_position = pparticle->Coordinates();
							noalias(iparticle->GetInitialPosition()) = old_position;
							iparticle->FastGetSolutionStepValue(DISPLACEMENT,1) = ZeroVector(3);
							const double small_dt = dt / subdivisions;
							//iparticle->FastGetSolutionStepValue(DISTANCE)=small_dt;
							//
							for (int substep = 0; substep < subdivisions; substep++)
								{
								typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
								Element::Pointer pelement;
								bool is_found = SearchStructure.FindPointOnMesh(current_position, N, pelement, result_begin, max_results);
								iparticle->Set(TO_ERASE, true);
								//KRATOS_WATCH(is_found);
								is_found = false;
								if (is_found == true)
									{
									Geometry< Node < 3 > >& geom = pelement->GetGeometry();
									//int nn=0;
									noalias(veulerian) = ZeroVector(3);
									//noalias(acc_particle) = ZeroVector(3);
									for (unsigned int k = 0; k < geom.size(); k++)
									{
										noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY);
									}


									noalias(current_position) += small_dt*veulerian;

									pparticle->Set(TO_ERASE, false);

									//iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY)=veulerian;
								} //if (is_found == true)
								else
									{
										//double time1=iparticle->FastGetSolutionStepValue(DISTANCE);
										/*array_1d<double,3> acc;
										acc[0] =  0.0;
										acc[1] = -10.0;
										acc[2] =  0.0;*/

										noalias(current_position) += small_dt * iparticle->FastGetSolutionStepValue(VELOCITY);
										pparticle->Set(TO_ERASE, false);
										iparticle->FastGetSolutionStepValue(IS_INTERFACE) = 1.0;
		  							}//else
	      					}//for

							iparticle->FastGetSolutionStepValue(DISPLACEMENT) = current_position - iparticle->GetInitialPosition();
      				}//for
					}
    			for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
      				{
							array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
							noalias(it->Coordinates()) = it->GetInitialPosition();
							noalias(it->Coordinates()) += dn1;
      				}
      }

void RungeKutta4ElementbasedSI(ModelPart& rModelPart, unsigned int substeps)
    {
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        dt *=0.5;
        //dt *=1.0;
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();

	//do movement
	array_1d<double, 3 > veulerian;
	//double temperature=0.0;
	//array_1d<double, 3 > acc_particle;

        array_1d<double, 3 > v1,v2,v3,v4,vtot,x;

	//array_1d<double, TDim + 1 > N;
	Vector N; N.resize(TDim + 1);

	const int max_results = 10000;

	typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

	const int nparticles = rModelPart.Nodes().size();
	
	bool element_is_active = true;

	#pragma omp parallel for firstprivate(results,N,veulerian,v1,v2,v3,v4,x)

	for (int i = 0; i < nparticles; i++)
	  {
	    array_1d<double,3> current_position;

	    array_1d<double,3> initial_position;

	    //bool is_found=false;
	    bool is_found1=false;
	    bool is_found2=false;
	    bool is_found3=false;
	    bool is_found4=false;

	   typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;

	    Node < 3 > ::Pointer pparticle = *(iparticle.base());

	    if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0) {
            initial_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT);

   	    //if(iparticle->FastGetSolutionStepValue(IS_FLUID) == 1.0) {

		Element::Pointer pelement;

		//STEP1
		//noalias(current_position) =  initial_position;

		is_found1 = SearchStructure.FindPointOnMesh(initial_position, N, pelement, result_begin, max_results);

		if(is_found1==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
				//if (pelement.Is(STRUCTURE)==true) double pepe=0.0;
				
				//element_is_active = true;
                               //if ((pelement)->IsDefined(STRUCTURE)==true)
                               element_is_active = (pelement)->Is(STRUCTURE);
                               //if(element_is_active==true) iparticle->Set(TO_ERASE, true);                            
                            
				noalias(v1) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v1) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else {
		noalias(v1) = pparticle->FastGetSolutionStepValue(VELOCITY);
		//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");
		}
		//STEP2
		noalias(x) =  initial_position + 0.5 * dt * v1;
		is_found2 = SearchStructure.FindPointOnMesh(x, N, pelement, result_begin, max_results);

		if(is_found2==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
				
				element_is_active = (pelement)->Is(STRUCTURE);
                               //if(element_is_active==true) iparticle->Set(TO_ERASE, true);      
                               
				noalias(v2) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v2) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else {
		noalias(v2) = pparticle->FastGetSolutionStepValue(VELOCITY);
		}

		//STEP3
		noalias(x) =  initial_position + 0.5 * dt * v2;
		is_found3 = SearchStructure.FindPointOnMesh(x, N, pelement, result_begin, max_results);

		if(is_found3==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
				
				element_is_active = (pelement)->Is(STRUCTURE);
                               //if(element_is_active==true) iparticle->Set(TO_ERASE, true);      
                               
				noalias(v3) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v3) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else {
		noalias(v3) = pparticle->FastGetSolutionStepValue(VELOCITY);
		//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");
		}
		//STEP4

		noalias(x) =  initial_position + dt * v3;
		is_found4 = SearchStructure.FindPointOnMesh(x, N, pelement, result_begin, max_results);
		if(is_found4==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
				
				element_is_active = (pelement)->Is(STRUCTURE);
                               //if(element_is_active==true) iparticle->Set(TO_ERASE, true);      
                               
				noalias(v4) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v4) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else {
		noalias(v4) = pparticle->FastGetSolutionStepValue(VELOCITY);
		//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");
		}

}

	       if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
		{

		if (iparticle->FastGetSolutionStepValue(IS_SOLID) == true) is_found1=false; //hypo=true
		
		if(is_found1==true)
			if(is_found2==true)
				if(is_found3==true)
					if(is_found4==true)
					{
					noalias(x) = initial_position;
			 		noalias(x) += 0.16666666666666666666667*dt*v1;
					noalias(x) += 0.33333333333333333333333*dt*v2;
					noalias(x) += 0.33333333333333333333333*dt*v3;
					noalias(x) += 0.16666666666666666666667*dt*v4;
					pparticle->FastGetSolutionStepValue(DISPLACEMENT) = x - pparticle->GetInitialPosition();

					}
		if(is_found1==false || is_found2==false || is_found3==false || is_found4==false)
		{
					noalias(x) = initial_position;
					noalias(x) += dt*pparticle->FastGetSolutionStepValue(VELOCITY);
					pparticle->FastGetSolutionStepValue(DISPLACEMENT) = x - pparticle->GetInitialPosition();
		}


}


	}//particles



		for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
			{
				array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
				noalias(it->Coordinates()) = it->GetInitialPosition();
				noalias(it->Coordinates()) += dn1;

			}


}



void RungeKutta4KernelbasedSI(ModelPart& rModelPart, unsigned int substeps)
    {
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        dt *=0.5;
        //dt *=1.0;
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();

	//do movement
	array_1d<double, 3 > veulerian;
	//double temperature=0.0;
	//array_1d<double, 3 > acc_particle;

        array_1d<double, 3 > v1,v2,v3,v4,vtot,x;


	const int nparticles = rModelPart.Nodes().size();
//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");

        //new code!
	typedef Node < 3 > PointType;
        typedef Node < 3 > ::Pointer PointTypePointer;
        typedef std::vector<PointType::Pointer> PointVector;
        typedef std::vector<PointType::Pointer>::iterator PointIterator;
        typedef std::vector<double> DistanceVector;
        typedef std::vector<double>::iterator DistanceIterator;
	//creating an auxiliary list for the new nodes
        PointVector list_of_nodes;
 	typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
	typedef Tree< KDTreePartition<BucketType> > tree; //Kdtree;




        for (ModelPart::NodesContainerType::iterator node_it = rModelPart.NodesBegin();
                node_it != rModelPart.NodesEnd(); ++node_it)
        {
            PointTypePointer pnode = *(node_it.base());

            //putting the nodes of the destination_model part in an auxiliary list
            list_of_nodes.push_back(pnode);
        }

	//create a spatial database with the list of new nodes
        unsigned int bucket_size = 20;
        tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);



	//#pragma omp parallel for firstprivate(results,N,veulerian,v1,v2,v3,v4,x)
	#pragma omp parallel for firstprivate(v1,v2,v3,v4,x)
	for (int i = 0; i < nparticles; i++)
	  {


	    array_1d<double,3> current_position;

	    array_1d<double,3> initial_position;

	    //bool is_found=false;
	    bool is_found1=false;
	    bool is_found2=false;
	    bool is_found3=false;
	    bool is_found4=false;


            double sigma = 0.0;
            if (TDim == 2)
                sigma = 10.0 / (7.0 * 3.1415926);
            else
                sigma = 1.0 / 3.1415926;



            unsigned int MaximumNumberOfResults = 10000;
            PointVector Results(MaximumNumberOfResults);
            DistanceVector SquaredResultsDistances(MaximumNumberOfResults);



	   //typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;

	   double radius = 1.5 * iparticle->FastGetSolutionStepValue(NODAL_H);


	    Node < 3 > ::Pointer pparticle = *(iparticle.base());

		if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0/*iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 0.0 or iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 1.0*/) {
		//if(iparticle->FastGetSolutionStepValue(IS_FLUID) == 1.0) {		
                initial_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT);

		bool flag_particle=false;
		
		if(iparticle->FastGetSolutionStepValue(IS_SOLID)==true)flag_particle=true;
		else flag_particle=false;
		
		
		Element::Pointer pelement;

		Node < 3 > work_point(0, 0.0, 0.0, 0.0);

		work_point.X() = initial_position(0);
		work_point.Y() = initial_position(1);
		work_point.Z() = initial_position(2);


		int number_of_points_in_radius;
	        number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),SquaredResultsDistances.begin(), MaximumNumberOfResults);

		if (number_of_points_in_radius > 0)
			{
				is_found1=true;
				double tot_weight = 0.0;
				noalias(v1)=ZeroVector(3);

				for (int k = 0; k < number_of_points_in_radius; k++)
					{

						double distance = sqrt(*(SquaredResultsDistances.begin() + k));
				    		double weight = SPHCubicKernel(sigma, distance, radius);

   		    				PointIterator it_found = Results.begin() + k;

                                               if(flag_particle==true) //solid
						     if((*it_found)->FastGetSolutionStepValue(IS_SOLID)==true)
							{
				    			tot_weight += weight;


							v1 += weight * (*it_found)->FastGetSolutionStepValue(VELOCITY);
							}
						else
							{
							tot_weight += weight;


							v1 += weight * (*it_found)->FastGetSolutionStepValue(VELOCITY);
							
							}	
							
							
							
							
					}
				v1 /= tot_weight;

			}
		else {

		noalias(v1) = pparticle->FastGetSolutionStepValue(VELOCITY);
		//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");
		}
		//STEP2
		noalias(x) =  initial_position + 0.5 * dt * v1;


                work_point.X() = x(0);
		work_point.Y() = x(1);
		work_point.Z() = x(2);


	        number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),SquaredResultsDistances.begin(), MaximumNumberOfResults);

	        if (number_of_points_in_radius > 0)
			{
				is_found2=true;
				double tot_weight = 0.0;
				noalias(v2)=ZeroVector(3);

				for (int k = 0; k < number_of_points_in_radius; k++)
					{

						double distance = sqrt(*(SquaredResultsDistances.begin() + k));
				    		double weight = SPHCubicKernel(sigma, distance, radius);

   		    				PointIterator it_found = Results.begin() + k;

                                               if(flag_particle==true) //solid
						     if((*it_found)->FastGetSolutionStepValue(IS_SOLID)==true)
							{
				    			tot_weight += weight;


							v2 += weight * (*it_found)->FastGetSolutionStepValue(VELOCITY);
							}
						else
							{
							tot_weight += weight;


							v2 += weight * (*it_found)->FastGetSolutionStepValue(VELOCITY);
							
							}	
					}
				v2 /= tot_weight;

			}
		else {

		//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");
		noalias(v2) = pparticle->FastGetSolutionStepValue(VELOCITY);
		}

		//STEP3
		noalias(x) =  initial_position + 0.5 * dt * v2;

		work_point.X() = x(0);
		work_point.Y() = x(1);
		work_point.Z() = x(2);


		number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),SquaredResultsDistances.begin(), MaximumNumberOfResults);


		//is_found2 = SearchStructure.FindPointOnMesh(x, N, pelement, result_begin, max_results);

	        if (number_of_points_in_radius > 0)
			{
				double tot_weight = 0.0;
				noalias(v3)=ZeroVector(3);
				is_found3 = true;
				for (int k = 0; k < number_of_points_in_radius; k++)
					{

						double distance = sqrt(*(SquaredResultsDistances.begin() + k));
				    		double weight = SPHCubicKernel(sigma, distance, radius);

   		    				PointIterator it_found = Results.begin() + k;

                                               if(flag_particle==true) //solid
						     if((*it_found)->FastGetSolutionStepValue(IS_SOLID)==true)
							{
				    			tot_weight += weight;


							v3 += weight * (*it_found)->FastGetSolutionStepValue(VELOCITY);
							}
						else
							{
							tot_weight += weight;


							v3 += weight * (*it_found)->FastGetSolutionStepValue(VELOCITY);
							
							}	
					}
				v3 /= tot_weight;

			}

		else {
		noalias(v3) = pparticle->FastGetSolutionStepValue(VELOCITY);
		}
		//STEP4

		noalias(x) =  initial_position + dt * v3;

		work_point.X() = x(0);
		work_point.Y() = x(1);
		work_point.Z() = x(2);

		number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),SquaredResultsDistances.begin(), MaximumNumberOfResults);


	        if (number_of_points_in_radius > 0)
			{
				is_found4 =true;
				double tot_weight = 0.0;
				noalias(v4)=ZeroVector(3);

				for (int k = 0; k < number_of_points_in_radius; k++)
					{

						double distance = sqrt(*(SquaredResultsDistances.begin() + k));
				    		double weight = SPHCubicKernel(sigma, distance, radius);

   		    				PointIterator it_found = Results.begin() + k;

                                               if(flag_particle==true) //solid
						     if((*it_found)->FastGetSolutionStepValue(IS_SOLID)==true)
							{
				    			tot_weight += weight;


							v4 += weight * (*it_found)->FastGetSolutionStepValue(VELOCITY);
							}
						else
							{
							tot_weight += weight;


							v4 += weight * (*it_found)->FastGetSolutionStepValue(VELOCITY);
							
							}	
					}
				v4 /= tot_weight;

			}
		else {
		noalias(v4) = pparticle->FastGetSolutionStepValue(VELOCITY);
		}

		}

	       if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0/*iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 0.0 or iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 1.0 */)
	       //if(iparticle->FastGetSolutionStepValue(IS_FLUID) == 1.0/*iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 0.0 or iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 1.0 */)

		{

		//if (iparticle->FastGetSolutionStepValue(IS_SOLID) == true) is_found1=false; //hypo=true
		//is_found1=false;
		if(is_found1==true)
			if(is_found2==true)
				if(is_found3==true)
					if(is_found4==true)
					{

					noalias(x) = initial_position;
			 		noalias(x) += 0.16666666666666666666667*dt*v1;
					noalias(x) += 0.33333333333333333333333*dt*v2;
					noalias(x) += 0.33333333333333333333333*dt*v3;
					noalias(x) += 0.16666666666666666666667*dt*v4;
					pparticle->FastGetSolutionStepValue(DISPLACEMENT) = x - pparticle->GetInitialPosition();

					}
			if(is_found1==false || is_found2==false || is_found3==false || is_found4==false)
					{
					noalias(x) = initial_position;
					noalias(x) += dt*pparticle->FastGetSolutionStepValue(VELOCITY);
					pparticle->FastGetSolutionStepValue(DISPLACEMENT) = x - pparticle->GetInitialPosition();
					}


		}

	//do_move

	}//particles

		for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
			{
				array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
				noalias(it->Coordinates()) = it->GetInitialPosition();
				noalias(it->Coordinates()) += dn1;

			}
}

void MovingParticlesN(ModelPart& rModelPart, unsigned int substeps)
    {
     KRATOS_TRY
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        //dt *=0.5;
        KRATOS_WATCH("jhgjhgjhgjhgjhgjhg")
        KRATOS_WATCH("jhgjhgjhgjhgjhgjhg")
        
        //KRATOS_THROW_ERROR(std::logic_error,"not fffffffffffffffffff","");
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();

	//do movement
	array_1d<double, 3 > veulerian;
	//double temperature=0.0;
	//array_1d<double, 3 > acc_particle;

        array_1d<double, 3 > v1,v2,v3,v4,vtot,x;

	//array_1d<double, TDim + 1 > N;

	//const int max_results = 10000;

	//typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

	const int nparticles = rModelPart.Nodes().size();
//KRATOS_THROW_ERROR(std::logic_error,"not fffffffffffffffffff","");

	//#pragma omp parallel for firstprivate(results,N,veulerian,v1,v2,v3,v4,x)
	#pragma omp parallel for firstprivate(veulerian,v1,v2,v3,v4,x)
	for (int i = 0; i < nparticles; i++)
	  {
	    array_1d<double,3> initial_position;

	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;

	    Node < 3 > ::Pointer pparticle = *(iparticle.base());

		if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0) {
                initial_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT,1);

 		noalias(x) = initial_position;
		noalias(x) += dt*pparticle->FastGetSolutionStepValue(VELOCITY,1);
		pparticle->FastGetSolutionStepValue(DISPLACEMENT) = x - pparticle->GetInitialPosition();
		KRATOS_WATCH(pparticle->GetInitialPosition())
		}


		}


		for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
			{
			if(it->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0) {
				array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
				//KRATOS_WATCH(it->FastGetSolutionStepValue(DISPLACEMENT,1))
				//KRATOS_WATCH(it->FastGetSolutionStepValue(DISPLACEMENT))
				noalias(it->Coordinates()) = it->GetInitialPosition();
				noalias(it->Coordinates()) += dn1;
				//KRATOS_THROW_ERROR(std::logic_error,"not fffffffffffffffffff","");
				}
			}

		for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
			{
			if(it->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0) {
				it->FastGetSolutionStepValue(DISPLACEMENT,1)=it->FastGetSolutionStepValue(DISPLACEMENT);
				KRATOS_WATCH(it->FastGetSolutionStepValue(DISPLACEMENT,1))
				KRATOS_WATCH(it->FastGetSolutionStepValue(DISPLACEMENT))
				//KRATOS_ERROR<<"error: this part is emptyyyy";
				//KRATOS_THROW_ERROR(std::logic_error,"not fffffffffffffffffff","");
				}
			}
			
KRATOS_CATCH("")
}
void MovingParticles(ModelPart& rModelPart, unsigned int substeps)
    {
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        //dt *=0.5;
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();
        //KRATOS_THROW_ERROR(std::logic_error,"not fffffffffffffffffff","");
	//do movement
	array_1d<double, 3 > veulerian;
	//double temperature=0.0;
	//array_1d<double, 3 > acc_particle;

        array_1d<double, 3 > v1,v2,v3,v4,vtot,x;

	//array_1d<double, TDim + 1 > N;

	//const int max_results = 10000;

	//typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

	const int nparticles = rModelPart.Nodes().size();

	//#pragma omp parallel for firstprivate(results,N,veulerian,v1,v2,v3,v4,x)
	#pragma omp parallel for firstprivate(veulerian,v1,v2,v3,v4,x)
/*	for (int i = 0; i < nparticles; i++)
	  {
	    array_1d<double,3> initial_position;

	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;

	    Node < 3 > ::Pointer pparticle = *(iparticle.base());

		if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0) {
                initial_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT);

		//HABRIA QUE GUARDAR EL DISPLACEMENT()
		//pparticle->FastGetSolutionStepValue(DISPLACEMENT,1)=pparticle->FastGetSolutionStepValue(DISPLACEMENT);
 		noalias(x) = initial_position;
		noalias(x) += dt*pparticle->FastGetSolutionStepValue(VELOCITY);
		pparticle->FastGetSolutionStepValue(DISPLACEMENT) = x - pparticle->GetInitialPosition();
		}


		}

*/
		for(ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
			{
			if(it->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0) {
				array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
				noalias(it->Coordinates()) = it->GetInitialPosition();
				noalias(it->Coordinates()) += dn1;

			}
}

}


    bool CheckInvertElement(ModelPart& ThisModelPart, int domain_size, double mesh_element_size)
    {
        KRATOS_TRY

        //set to zero the nodal area
        bool inverted=false;
	double vol=0.0;
	array_1d<double,TDim+1> N;
        //double dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];
	//defining work arrays
        //PointerVector< Element > elements_to_solve;
        std::vector<GlobalPointersVector<Element>> nodal_neighbours;
        //elements_to_solve.reserve(ThisModelPart.Elements().size());

        VariableUtils().SetFlag(TO_ERASE, false, ThisModelPart.Elements());

        VariableUtils().SetFlag(TO_ERASE, false, ThisModelPart.Nodes());
        //VariableUtils().SetFlag(INTERFACE, 0, ThisModelPart.Nodes());
        //rDestinationModelPart.RemoveNodesFromAllLevels(TO_ERASE);
        //rDestinationModelPart.RemoveElementsFromAllLevels(TO_ERASE);

        if(domain_size == 2)
        {
	    KRATOS_ERROR<<"error: this part is emptyyyy";
        }
        else
        {

	      //reference volumen
             double vol_r=sqrt(2)* pow(mesh_element_size, 3)	/12.0;
	     for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
                in!=ThisModelPart.NodesEnd(); in++)
		{
		    in->FastGetSolutionStepValue(RADIATIVE_INTENSITY) = 0.00;
		}

	      //bool erase_node=false;
             for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin(); i!=ThisModelPart.ElementsEnd(); i++)
		    {
		        //calculating shape functions values
		        Geometry< Node<3> >& geom = i->GetGeometry();


		        //counting number of structural nodes
		        vol = geom.Volume();

//	     	        if(vol <= 0 or vol<= 0.00000000000001)
     	     	        //if(vol <= 0 or vol<= 0.0001*vol_r)
     	     	        //if(vol <= 0 or vol<= 0.001*vol_r)

     	     	        /////

			//erase_node=false;
                    	//for (unsigned int k = 0; k < geom.size(); k++)
			//{
            		//array_1d<double,3> delta_disp = geom[k].FastGetSolutionStepValue(DISPLACEMENT);
            		//noalias(delta_disp) -= geom[k].FastGetSolutionStepValue(DISPLACEMENT,1);

            		//double norm_delta_disp = norm_2(delta_disp);

            		//array_1d<double,3> v_old = geom[k].FastGetSolutionStepValue(VELOCITY,1);
            		//double norm_v = norm_2(v_old);

            		//if(norm_delta_disp*3.0 < norm_v*dt ) erase_node=true;
                       //if(norm_delta_disp* (0.333333333333333*0.001) >  norm_v*dt ) erase_node=true;
                       //if(norm_v >1.5) erase_node=true;
        		//}



     	     	        /////
     	     	        if(vol<= 0.01*vol_r || vol> 2.5*vol_r )//or erase_node==true)  //0.001

	     	        {
	     	        	inverted=true;

	     	        	i->Set(TO_ERASE, true);

	     	        	geom[0].FastGetSolutionStepValue(RADIATIVE_INTENSITY)=1.0;
				geom[1].FastGetSolutionStepValue(RADIATIVE_INTENSITY)=1.0;
				geom[2].FastGetSolutionStepValue(RADIATIVE_INTENSITY)=1.0;
				geom[3].FastGetSolutionStepValue(RADIATIVE_INTENSITY)=1.0;

				auto& r_nodes_current_element = i->GetGeometry();

    				auto& pNode0 = r_nodes_current_element[0];
    				auto& pNode1 = r_nodes_current_element[1];
    				auto& pNode2 = r_nodes_current_element[2];
    				auto& pNode3 = r_nodes_current_element[3];

				// Neighbour elements of each node of the current element
				GlobalPointersVector<Element>& r_neigh_node_0 = pNode0.GetValue(NEIGHBOUR_ELEMENTS);
				GlobalPointersVector<Element>& r_neigh_node_1 = pNode1.GetValue(NEIGHBOUR_ELEMENTS);
				GlobalPointersVector<Element>& r_neigh_node_2 = pNode2.GetValue(NEIGHBOUR_ELEMENTS);
				GlobalPointersVector<Element>& r_neigh_node_3 = pNode3.GetValue(NEIGHBOUR_ELEMENTS);


				nodal_neighbours.push_back(r_neigh_node_0);
				nodal_neighbours.push_back(r_neigh_node_1);
				nodal_neighbours.push_back(r_neigh_node_2);
				nodal_neighbours.push_back(r_neigh_node_3);

				bool is_inside_0 = false;
		    		bool is_inside_1 = false;
		    		bool is_inside_2 = false;
		    		bool is_inside_3 = false;
		    		bool all_neigh_0 =false;
		    		bool all_neigh_1 =false;
		    		bool all_neigh_2 =false;
		    		bool all_neigh_3 =false;

				//first node
           			 for (unsigned int neigh_elem = 0; neigh_elem < nodal_neighbours.size(); neigh_elem++) { //loop for nodes
            			// Nodes of the neigh element
            				for (unsigned int elem = 0; elem < nodal_neighbours[neigh_elem].size(); elem++) //loop for elements
            				{
            				Element& parent = nodal_neighbours[neigh_elem][elem];

            				/*auto& Node0 = parent.GetGeometry()[0];
    					auto& Node1 = parent.GetGeometry()[1];
    					auto& Node2 = parent.GetGeometry()[2];
    					auto& Node3 = parent.GetGeometry()[3];*/
    					is_inside_0 = false;
		    			is_inside_1 = false;
		    			is_inside_2 = false;
		    			is_inside_3 = false;

            				is_inside_0 = CalculatePosition(parent.GetGeometry(),pNode0[0],pNode0[1],pNode0[2],N); //observe if node is inside the neighbours!
            				is_inside_1 = CalculatePosition(parent.GetGeometry(),pNode1[0],pNode1[1],pNode1[2],N);
            				is_inside_2 = CalculatePosition(parent.GetGeometry(),pNode2[0],pNode2[1],pNode2[2],N);
            				is_inside_3 = CalculatePosition(parent.GetGeometry(),pNode3[0],pNode3[1],pNode3[2],N);
            				/*if(is_inside_0 == true) parent.Set(TO_ERASE, true);
            				if(is_inside_1 == true) parent.Set(TO_ERASE, true);
            				if(is_inside_2 == true) parent.Set(TO_ERASE, true);
            				if(is_inside_3 == true) parent.Set(TO_ERASE, true);*/

            				if(is_inside_0 == true) {all_neigh_0=true; pNode0.Set(TO_ERASE, true);}
            				if(is_inside_1 == true) {all_neigh_1=true;pNode1.Set(TO_ERASE, true);}
            				if(is_inside_2 == true) {all_neigh_2=true;pNode2.Set(TO_ERASE, true);}
            				if(is_inside_3 == true) {all_neigh_3=true;pNode3.Set(TO_ERASE, true);}

            				//KRATOS_WATCH("--------------------------->")
            				//KRATOS_WATCH("--------------------------->")
            				//KRATOS_WATCH("--------------------------->")
          				//KRATOS_WATCH("NODOSSSS PADRESSSSSSSSSSSSSSSSSSS")
            				//KRATOS_WATCH(pNode0)
            				//KRATOS_WATCH(pNode1)
            				//KRATOS_WATCH(pNode2)
            				//KRATOS_WATCH(pNode3)
            				//KRATOS_WATCH("ELEMENTOOOOOOOOOO VECINOOOOOOOOOOOO")
            				//KRATOS_WATCH(parent.GetGeometry())
                			}

                		}
                		if(all_neigh_0==true)
                		{
                		for (unsigned int elem = 0; elem < nodal_neighbours[0].size(); elem++) //loop for elements
            				{
            					Element& parent = nodal_neighbours[0][elem];
            					parent.Set(TO_ERASE, true);
            				}
                		}

                		if(all_neigh_1==true)
                		{
                		for (unsigned int elem = 0; elem < nodal_neighbours[1].size(); elem++) //loop for elements
            				{
            					Element& parent = nodal_neighbours[1][elem];
            					parent.Set(TO_ERASE, true);
            				}
                		}
                		if(all_neigh_2==true)
                		{
                		for (unsigned int elem = 0; elem < nodal_neighbours[2].size(); elem++) //loop for elements
            				{
            					Element& parent = nodal_neighbours[2][elem];
            					parent.Set(TO_ERASE, true);
            				}
                		}
                		if(all_neigh_3==true)
                		{
                		for (unsigned int elem = 0; elem < nodal_neighbours[3].size(); elem++) //loop for elements
            				{
            					Element& parent = nodal_neighbours[3][elem];
            					parent.Set(TO_ERASE, true);
            				}
                		}
				//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");


			nodal_neighbours.clear();
	     	        }
	     	        //KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");
		    }
	}

	  ThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);
          ThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);


	return inverted;

        KRATOS_CATCH("")
    }


    double CalculateVolume(ModelPart& ThisModelPart, int domain_size)
    {
        KRATOS_TRY

        //set to zero the nodal area
        //bool inverted=false;
	double vol=0.0;
	double totvol=0.0;
	array_1d<double,TDim+1> N;


        if(domain_size == 2)
        {
	    KRATOS_ERROR<<"error: this part is emptyyyy";
        }
        else
        {

            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
		    {
		        //calculating shape functions values
		        Geometry< Node<3> >& geom = i->GetGeometry();

		        //counting number of structural nodes
		        vol = geom.Volume();
                       totvol += vol;
		    }
	}



	return totvol;

        KRATOS_CATCH("")
    }


    void Aux(ModelPart& ThisModelPart, int domain_size)
    {
        KRATOS_TRY

           Matrix dummy=ZeroMatrix(2,2);
           double Area;
           array_1d<double, 3> N;

           BoundedMatrix<double, 3, 2> msDN_Dx;

           BoundedMatrix<double,3,3> Sigma=ZeroMatrix(2,2);

	    BoundedMatrix<double,2,2> CauchyStress=ZeroMatrix(2,2);


             for(ModelPart::ElementsContainerType::iterator im = ThisModelPart.ElementsBegin() ; im != ThisModelPart.ElementsEnd() ; ++im)
                 {
                 dummy=ZeroMatrix(2,2);
                 im->SetValue(CAUCHY_STRESS_TENSOR, dummy);
}
 
        KRATOS_CATCH("")
    }
    private:
      inline double SPHCubicKernel(const double sigma, const double r, const double hmax)
    {


        double h_half = 0.5 * hmax;
        const double s = r / h_half;
        const double coeff = sigma / pow(h_half, static_cast<int>(TDim));

        if (s <= 1.0)
            return coeff * (1.0 - 1.5 * s * s + 0.75 * s * s * s);
        else if (s <= 2.0)
            return 0.25 * coeff * pow(2.0 - s, 3);
        else
            return 0.0;
    }

 inline bool CalculatePosition(
        Geometry<Node < 3 > >&geom,
        const double xc,
         const double yc,
         const double zc,
        array_1d<double, 4 > & N
        )
    {

        const double x0 = geom[0].X();
        const double y0 = geom[0].Y();
        const double z0 = geom[0].Z();
        const double x1 = geom[1].X();
        const double y1 = geom[1].Y();
        const double z1 = geom[1].Z();
        const double x2 = geom[2].X();
        const double y2 = geom[2].Y();
        const double z2 = geom[2].Z();
        const double x3 = geom[3].X();
        const double y3 = geom[3].Y();
        const double z3 = geom[3].Z();

        const double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

        double inv_vol = 0.0;
        if (vol < 0.0000000000001)
        {
            //The interpolated node will not be inside an elemente with zero area
            return false;
        }
        else
        {
            inv_vol = 1.0 / vol;
        }

        N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
        N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;
        N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
        N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;


        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >= 0.0 &&
                N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <= 1.0)
            //if the xc yc zc is inside the tetrahedron return true
            return true;

        return false;
    }
    /**
     * This method computes the volume of a tetrahedra
     */
    inline double CalculateVol(const double x0, const double y0, const double z0,
                               const double x1, const double y1, const double z1,
                               const double x2, const double y2, const double z2,
                               const double x3, const double y3, const double z3
                              )
    {
        const double x10 = x1 - x0;
        const double y10 = y1 - y0;
        const double z10 = z1 - z0;

        const double x20 = x2 - x0;
        const double y20 = y2 - y0;
        const double z20 = z2 - z0;

        const double x30 = x3 - x0;
        const double y30 = y3 - y0;
        const double z30 = z3 - z0;

        const double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
        return detJ * 0.1666666666666666666667;
    }



};

} // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED  defined



