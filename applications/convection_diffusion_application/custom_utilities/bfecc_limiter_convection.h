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
				iparticle->GetValue(ERROR_1) = 0.5*iparticle->FastGetSolutionStepValue(rVar,1) - 0.5*phi_old;//computing error1 as e1 = 0.5*(rVar(n) - phi_old)
                iparticle->GetValue(rVar) = iparticle->FastGetSolutionStepValue(rVar,1) + iparticle->GetValue(ERROR_1);//rVar(n)+e1
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
		//computing error 2 with forward of phi1
		#pragma omp parallel for firstprivate(results,N)
        for (int i = 0; i < nparticles; i++)
        {
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
            
            ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
            
            Element::Pointer pelement;
            array_1d<double,3> fwdPos = iparticle->Coordinates();
            const array_1d<double,3>& vel = iparticle->FastGetSolutionStepValue(conv_var);
            bool is_found = ConvectBySubstepping(dt,fwdPos,vel, N, pelement, result_begin, max_results, 1.0, substeps);//seeking forwards
            found[i] = is_found;
            
            if(is_found) {                
				//Forward with
                Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                double phi2 = N[0] * ( geom[0].FastGetSolutionStepValue(rVar));
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi2 += N[k] * ( geom[k].FastGetSolutionStepValue(rVar) );
                }       
				//Computing error2 as e2 = rVar(n)-(phi2+e1)
			    iparticle->GetValue(ERROR_2) = iparticle->FastGetSolutionStepValue(rVar,1) - (phi2 + iparticle->GetValue(ERROR_1));
			}
			
        }
		//Limiter
		FindNodalNeighboursProcess adj_process(rModelPart, 10, 10);
		adj_process.ClearNeighbours();
		adj_process.Execute();

		#pragma omp parallel for
        for (int i = 0; i < nparticles; i++){
            ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
			iparticle->GetValue(ERROR) = iparticle->GetValue(ERROR_1);
		}
		
		#pragma omp parallel for
        for (int i = 0; i < nparticles; i++){
			ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
			WeakPointerVector<Node<3> > & adj_iparticle = iparticle->GetValue(NEIGHBOUR_NODES);

			//double errorA = std::abs(iparticle->GetValue(ERROR_1));
			//double errorB = std::abs(iparticle->GetValue(ERROR_2));
			//double diff = errorB - errorA;
						
			//if( diff > std::numeric_limits<double>::epsilon() ){
			if(std::abs(iparticle->GetValue(ERROR_2)) > std::abs(iparticle->GetValue(ERROR_1))){
				for(int j = 0 ; j < iparticle->GetValue(NEIGHBOUR_NODES).size(); j++){
					adj_iparticle[j].GetValue(ERROR) = minmod(iparticle->GetValue(ERROR_1),adj_iparticle[j].GetValue(ERROR));
				}
			}
		}
		#pragma omp parallel for
        for (int i = 0; i < nparticles; i++){
			ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
			bool is_found = foundf[i];
            if(is_found) {
				iparticle->GetValue(rVar) = iparticle->FastGetSolutionStepValue(rVar,1) + iparticle->GetValue(ERROR);
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


