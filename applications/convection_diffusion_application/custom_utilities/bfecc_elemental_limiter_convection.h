// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

#if !defined(KRATOS_BFECC_LIMITER_CONVECTION_INCLUDED )
#define  KRATOS_BFECC_LIMITER_CONVECTION_INCLUDED

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
#include "includes/variables.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"
#include "processes/find_nodal_neighbours_process.h"


#include <boost/timer.hpp>
#include "utilities/timer.h"

#ifdef _OPENMP
#include "omp.h"
#endif



namespace Kratos
{

template<std::size_t TDim> class BFECCLimiterConvection
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BFECCLimiterConvection<TDim>);

    BFECCLimiterConvection(typename BinBasedFastPointLocator<TDim>::Pointer pSearchStructure)
        : mpSearchStructure(pSearchStructure)
    {
    }

    ~BFECCLimiterConvection()
    {
    }

    //**********************************************************************************************
    //**********************************************************************************************
    void BFECCconvect(ModelPart& rModelPart, const Variable< double >& rVar, const Variable<array_1d<double,3> >& conv_var, const double substeps)
    {
        KRATOS_TRY
        const double dt = rModelPart.GetProcessInfo()[DELTA_TIME];

        //do movement
        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        const int nparticles = rModelPart.Nodes().size();
         
        PointerVector< Element > elem_backward( rModelPart.Nodes().size());
        std::vector< array_1d<double,TDim+1> > Ns( rModelPart.Nodes().size());
        std::vector< bool > found( rModelPart.Nodes().size());
		std::vector< bool > foundf( rModelPart.Nodes().size());
        
        //FIRST LOOP: estimate rVar(n+1) 
        #pragma omp parallel for firstprivate(results,N)
        for (int i = 0; i < nparticles; i++)
        {
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
            
            ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
            
            Element::Pointer pelement;
            array_1d<double,3> bckPos = iparticle->Coordinates();
            const array_1d<double,3>& vel = iparticle->FastGetSolutionStepValue(conv_var);
            bool is_found = ConvectBySubstepping(dt,bckPos,vel, N, pelement, result_begin, max_results, -1.0, substeps);
            found[i] = is_found;
            
            if(is_found) {
                //save position backwards
                elem_backward(i) = pelement;
                Ns[i] = N;
                
                Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                double phi1 = N[0] * ( geom[0].FastGetSolutionStepValue(rVar,1));
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi1 += N[k] * ( geom[k].FastGetSolutionStepValue(rVar,1) );
                }
                
                iparticle->FastGetSolutionStepValue(rVar) = phi1;
            }
        }
        
        //now obtain the value AT TIME STEP N by taking it from N+1
        #pragma omp parallel for firstprivate(results,N)
        for (int i = 0; i < nparticles; i++)
        {
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
            
            ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
            
            Element::Pointer pelement;
            array_1d<double,3> fwdPos = iparticle->Coordinates();
            const array_1d<double,3>& vel = iparticle->FastGetSolutionStepValue(conv_var,1);
            bool is_found = ConvectBySubstepping(dt,fwdPos,vel, N, pelement, result_begin, max_results, 1.0, substeps);
            foundf[i] = is_found;

            if(is_found) {
                Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                double phi_old = N[0] * ( geom[0].FastGetSolutionStepValue(rVar));
                
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi_old  += N[k] * ( geom[k].FastGetSolutionStepValue(rVar) );
                }
                
                //Computing error 1 and modified solution at time N to be interpolated again
				iparticle->GetValue(BFECC_ERROR_1) = 0.5*iparticle->FastGetSolutionStepValue(rVar,1) - 0.5*phi_old;//computing error1 as e1 = 0.5*(rVar(n) - phi_old)
                iparticle->GetValue(rVar) = iparticle->FastGetSolutionStepValue(rVar,1) + iparticle->GetValue(BFECC_ERROR_1);//rVar(n)+e1
			}
        }
		//Backward with modified solution
         #pragma omp parallel for 
        for (int i = 0; i < nparticles; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
            bool is_found = found[i];
            if(is_found) {
                array_1d<double,TDim+1> N = Ns[i];
                Geometry< Node < 3 > >& geom = elem_backward[i].GetGeometry();
                double phi1 = N[0] * ( geom[0].GetValue(rVar));
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi1 += N[k] * ( geom[k].GetValue(rVar) );
                }   
                iparticle->FastGetSolutionStepValue(rVar) = phi1;
            }
//             else
//                 std::cout << "it should find it" << std::endl;
        }



		// computing error 2 with forward of phi1
		int nelements = rModelPart.NumberOfElements();
		for(int i = 0 ; i < nelements; i++)
		{
			typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

			ModelPart::ElementsContainerType::iterator i_element = rModelPart.ElementsBegin() + i;
			Element::GeometryType& element_geometry = i_element->GetGeometry();	

			Element::Pointer pelement;
            array_1d<double,3> fwdPos = i_element->GetGeometry().Center();
            array_1d<double,3> vel = ZeroVector(3);

			for(unsigned int j = 0 ; j < element_geometry.size(); j++){
				for(int k = 0 ; k < 3; k++){
					vel[k] += element_geometry[j].GetSolutionStepValue(conv_var)[k] / element_geometry.size();
				}
			}

            bool is_found = ConvectBySubstepping(dt,fwdPos,vel, N, pelement, result_begin, max_results, 1.0, substeps);//seeking forwards
			double e1 = 0.00f;
			
			for(unsigned int j = 0 ; j < element_geometry.size(); j++){ 
				e1 += element_geometry[j].GetValue(BFECC_ERROR_1);
			}
			e1 /= element_geometry.size(); 

			double e2 = e1;
            if(is_found) {                
				//Forward with
                Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                double phi2 = N[0] * ( geom[0].FastGetSolutionStepValue(rVar));
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi2 += N[k] * ( geom[k].FastGetSolutionStepValue(rVar) );
                }
				double solution_in_center = 0;
				//Computing error2 as e2 = rVar(n)-(phi2+e1)
				for(unsigned int j = 0 ; j < element_geometry.size(); j++){
					solution_in_center += i_element->GetGeometry()[j].FastGetSolutionStepValue(rVar,1);
				}
				solution_in_center /= element_geometry.size(); 

				e2 = solution_in_center - (phi2 + e1);

				if(std::abs(e2) > std::abs(e1)){
					for(unsigned int j = 0 ; j < element_geometry.size(); j++){
						element_geometry[j].GetValue(BFECC_ERROR) = minmod(e1,element_geometry[j].GetValue(BFECC_ERROR_1));
					}
				}
			}
			
		}

		#pragma omp parallel for
        for (int i = 0; i < nparticles; i++){
			ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
			bool is_found = foundf[i];
            if(is_found) {
				iparticle->GetValue(rVar) = iparticle->FastGetSolutionStepValue(rVar,1) + iparticle->GetValue(BFECC_ERROR);
			}
		}
		//Backward with modified solution
		#pragma omp parallel for 
        for (int i = 0; i < nparticles; i++){
            ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
            bool is_found = found[i];
            if(is_found) {
                array_1d<double,TDim+1> N = Ns[i];
                Geometry< Node < 3 > >& geom = elem_backward[i].GetGeometry();
                double phi1 = N[0] * ( geom[0].GetValue(rVar));
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi1 += N[k] * ( geom[k].GetValue(rVar) );
                }
                iparticle->FastGetSolutionStepValue(rVar) = phi1;
            }
        }

        KRATOS_CATCH("")
    }

	double minmod(double x, double y) {
		double f;
		if(x > 0.0f && y > 0.0f)
			f = std::min(x,y); 
		else if(x < 0.0f && y < 0.0f)
			f = std::max(x,y); 
		else
			f = 0;
		return f;
	}
	
    bool ConvectBySubstepping(
                const double dt,
                 array_1d<double,3>& position, //IT WILL BE MODIFIED
                 const array_1d<double,3>& initial_velocity, 
                 array_1d<double,TDim+1>& N, 
                 Element::Pointer& pelement, 
                 typename BinBasedFastPointLocator<TDim>::ResultIteratorType& result_begin,
                 const unsigned int max_results,
                 const double velocity_sign,
                 const double subdivisions)
    {
        bool is_found = false;
        array_1d<double,3> veulerian;
        const double small_dt = dt/subdivisions;
        
        
        if(velocity_sign > 0.0) //going from the past to the future
        {
            noalias(position) += small_dt*initial_velocity;
            unsigned int substep=0;
            while(substep++ < subdivisions)
            {
                is_found = mpSearchStructure->FindPointOnMesh(position, N, pelement, result_begin, max_results);

                if (is_found == true)
                {
                    Geometry< Node < 3 > >& geom = pelement->GetGeometry();

                    const double new_step_factor = static_cast<double>(substep)/subdivisions;
                    const double old_step_factor = (1.0 - new_step_factor);

                    noalias(veulerian) = N[0] * ( new_step_factor*geom[0].FastGetSolutionStepValue(VELOCITY) + old_step_factor*geom[0].FastGetSolutionStepValue(VELOCITY,1));
                    for (unsigned int k = 1; k < geom.size(); k++)
                        noalias(veulerian) += N[k] * ( new_step_factor*geom[k].FastGetSolutionStepValue(VELOCITY) + old_step_factor*geom[k].FastGetSolutionStepValue(VELOCITY,1) );
                    
                    noalias(position) += small_dt*veulerian;

 
                }    
                else
                    break;
            }
        }
        else //going from the future to the past
        {
            noalias(position) -= small_dt*initial_velocity;
            unsigned int substep=0;
            while(substep++ < subdivisions)
            {
                is_found = mpSearchStructure->FindPointOnMesh(position, N, pelement, result_begin, max_results);

                if (is_found == true)
                {
                    Geometry< Node < 3 > >& geom = pelement->GetGeometry();

                    //this factors get inverted from the other case
                   const double old_step_factor = static_cast<double>(substep)/subdivisions;
                   const double new_step_factor = (1.0 - old_step_factor);
 
                    noalias(veulerian) = N[0] * ( new_step_factor*geom[0].FastGetSolutionStepValue(VELOCITY) + old_step_factor*geom[0].FastGetSolutionStepValue(VELOCITY,1));
                    for (unsigned int k = 1; k < geom.size(); k++)
                        noalias(veulerian) += N[k] * ( new_step_factor*geom[k].FastGetSolutionStepValue(VELOCITY) + old_step_factor*geom[k].FastGetSolutionStepValue(VELOCITY,1) );
                    
                    noalias(position) -= small_dt*veulerian;

 
                }         
             else 
                 break;
            }
        }
        
                return is_found;
        
    }
    
private:
    typename BinBasedFastPointLocator<TDim>::Pointer mpSearchStructure;



};

} // namespace Kratos.

#endif // KRATOS_BFECC_LIMITER_CONVECTION_INCLUDED  defined


