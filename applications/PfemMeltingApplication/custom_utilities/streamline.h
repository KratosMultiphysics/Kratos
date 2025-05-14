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
//#include "processes/node_erase_process.h"
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



void MovingParticlesN(ModelPart& rModelPart, unsigned int substeps)
    {
     KRATOS_TRY
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        //dt *=0.5;
        //KRATOS_WATCH("jhgjhgjhgjhgjhgjhg")
        //KRATOS_WATCH("jhgjhgjhgjhgjhgjhg")
        
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
        ///KRATOS_THROW_ERROR(std::logic_error,"not fffffffffffffffffff","");

	//#pragma omp parallel for firstprivate(results,N,veulerian,v1,v2,v3,v4,x)
	#pragma omp parallel for firstprivate(veulerian,v1,v2,v3,v4,x)
	for (int i = 0; i < nparticles; i++)
	  {
	    array_1d<double,3> initial_position;

	    ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;

	    Node ::Pointer pparticle = *(iparticle.base());

		if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0) {
                initial_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT,1);

 		noalias(x) = initial_position;
		noalias(x) += dt*pparticle->FastGetSolutionStepValue(VELOCITY,1);
		pparticle->FastGetSolutionStepValue(DISPLACEMENT) = x - pparticle->GetInitialPosition();
		//KRATOS_WATCH(pparticle->GetInitialPosition())
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
				//it->FastGetSolutionStepValue(DISPLACEMENT,1)=it->FastGetSolutionStepValue(DISPLACEMENT);
				array_1d<double,3>& dn1 = it->FastGetSolutionStepValue(DISPLACEMENT);
				array_1d<double,3>& dn2 = it->FastGetSolutionStepValue(DISPLACEMENT,1);
				dn2=dn1;
				//KRATOS_WATCH(it->FastGetSolutionStepValue(DISPLACEMENT,1))
				//KRATOS_WATCH(it->FastGetSolutionStepValue(DISPLACEMENT))
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

void RungeKutta4ElementbasedSI(ModelPart& rModelPart, unsigned int substeps)
    {
	double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        //dt *=0.5;
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

	    Node ::Pointer pparticle = *(iparticle.base());

	    if(iparticle->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0) {
            initial_position = pparticle->GetInitialPosition() + pparticle->FastGetSolutionStepValue(DISPLACEMENT);
            
            //KRATOS_WATCH(initial_position)	
            //KRATOS_WATCH(pparticle->FastGetSolutionStepValue(DISPLACEMENT))	
            
   	    //if(iparticle->FastGetSolutionStepValue(IS_FLUID) == 1.0) {

		Element::Pointer pelement;

		//STEP1
		//noalias(current_position) =  initial_position;

		is_found1 = SearchStructure.FindPointOnMesh(initial_position, N, pelement, result_begin, max_results);

		if(is_found1==true)
			{
				Geometry< Node >& geom = pelement->GetGeometry();
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
				Geometry< Node >& geom = pelement->GetGeometry();
				
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
				Geometry< Node >& geom = pelement->GetGeometry();
				
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
				Geometry< Node >& geom = pelement->GetGeometry();
				
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
				
				//KRATOS_WATCH(it->GetInitialPosition())	
            			//KRATOS_WATCH(it->Coordinates())
            			//KRATOS_WATCH(it->FastGetSolutionStepValue(DISPLACEMENT))

			}


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
        Geometry<Node >&geom,
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



