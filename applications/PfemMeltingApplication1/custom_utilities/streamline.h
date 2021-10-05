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
//#include "utilities/enrichment_utilities.h"

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
				double dt=0.5 * deltadt;

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
					//bool first_time=false;
					iparticle->FastGetSolutionStepValue(DISTANCE)=0.0;

					iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY) = iparticle->FastGetSolutionStepValue(VELOCITY);  //AUX_VEL

					if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 1.0) do_move = false;

					if( do_move == true  ) //note that we suppose the velocity components to be all fixed
						{
							array_1d<double,3> old_position = pparticle->Coordinates();
							array_1d<double,3> current_position = pparticle->Coordinates();
							noalias(iparticle->GetInitialPosition()) = old_position;
							iparticle->FastGetSolutionStepValue(DISPLACEMENT,1) = ZeroVector(3);
							subdivisions=4;
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
										noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY);
									}


									noalias(current_position) += small_dt*veulerian;

									pparticle->Set(TO_ERASE, false);
									iparticle->FastGetSolutionStepValue(DISTANCE) += small_dt;
									iparticle->FastGetSolutionStepValue(EMBEDDED_VELOCITY)=veulerian;
								} //if (is_found == true)
								else
									{
										//double time1=iparticle->FastGetSolutionStepValue(DISTANCE);
										array_1d<double,3> acc;
										acc[0] =  0.0;
										acc[1] = -10.0;
										acc[2] =  0.0;

										noalias(current_position) += small_dt * iparticle->FastGetSolutionStepValue(VELOCITY);									pparticle->Set(TO_ERASE, false);
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
        //dt *=0.5; 
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

	    //bool is_found=false;
	    bool is_found1=false;
	    bool is_found2=false;
	    bool is_found3=false;
	    bool is_found4=false;


	   typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;


	    Node < 3 > ::Pointer pparticle = *(iparticle.base());

if(iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 0.0) {
                initial_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT);

		Element::Pointer pelement;

		//STEP1
		//noalias(current_position) =  initial_position;

		is_found1 = SearchStructure.FindPointOnMesh(initial_position, N, pelement, result_begin, max_results);

		if(is_found1==true)
			{
				Geometry< Node < 3 > >& geom = pelement->GetGeometry();
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
				noalias(v4) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                		for (unsigned int k = 1; k < geom.size(); k++)
                    		noalias(v4) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
			}
		else {
		noalias(v4) = pparticle->FastGetSolutionStepValue(VELOCITY);
		//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");
		}

}

	       if(iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 0.0)
		{

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
	BinBasedFastPointLocator<TDim> SearchStructure(rModelPart);
	SearchStructure.UpdateSearchDatabase();

	//do movement
	array_1d<double, 3 > veulerian;
	//double temperature=0.0;
	//array_1d<double, 3 > acc_particle;

        array_1d<double, 3 > v1,v2,v3,v4,vtot,x;


	const int nparticles = rModelPart.Nodes().size();


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

if(iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 0.0) {
                initial_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT);

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

				    		tot_weight += weight;

				    		PointIterator it_found = Results.begin() + k;

						v1 += weight * (*it_found)->FastGetSolutionStepValue(VELOCITY);

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

				    		tot_weight += weight;

				    		PointIterator it_found = Results.begin() + k;

						v2 += weight * (*it_found)->FastGetSolutionStepValue(VELOCITY);

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

				    		tot_weight += weight;

				    		PointIterator it_found = Results.begin() + k;

						v3 += weight * (*it_found)->FastGetSolutionStepValue(VELOCITY);

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

				    		tot_weight += weight;

				    		PointIterator it_found = Results.begin() + k;

						v4 += weight * (*it_found)->FastGetSolutionStepValue(VELOCITY);

					}
				v4 /= tot_weight;

			}
		else {
		noalias(v4) = pparticle->FastGetSolutionStepValue(VELOCITY);
		}

		}

	       if(iparticle->FastGetSolutionStepValue(IS_INTERFACE) == 0.0)
		{


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

    bool CheckInvertElement(ModelPart& ThisModelPart, int domain_size)
    {
        KRATOS_TRY

        //set to zero the nodal area
        bool inverted=false;
	double vol=0.0;
        
        if(domain_size == 2)
        {
	    KRATOS_ERROR<<"error: this part is emptyyyy";
        }
        else 
        {

	     for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
                in!=ThisModelPart.NodesEnd(); in++)
		{
		    in->FastGetSolutionStepValue(RADIATIVE_INTENSITY) = 0.00;
		}
            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
		    {
		        //calculating shape functions values
		        Geometry< Node<3> >& geom = i->GetGeometry();
		        //counting number of structural nodes
		        vol = GeometryUtils::CalculateVolume3D(geom);

			geom[0].FastGetSolutionStepValue(RADIATIVE_INTENSITY)+=0.25*vol;
			geom[1].FastGetSolutionStepValue(RADIATIVE_INTENSITY)+=0.25*vol;
			geom[2].FastGetSolutionStepValue(RADIATIVE_INTENSITY)+=0.25*vol;
			geom[3].FastGetSolutionStepValue(RADIATIVE_INTENSITY)+=0.25*vol;

	     	        if(vol <= 0)  inverted=true;	
		    }
	}
	return inverted;

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
};

} // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED  defined



