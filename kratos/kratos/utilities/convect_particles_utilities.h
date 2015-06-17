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


#if !defined(KRATOS_CONVECT_PARTICLES_UTILITIES_INCLUDED )
#define  KRATOS_CONVECT_PARTICLES_UTILITIES_INCLUDED

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


#include <boost/timer.hpp>
#include "utilities/timer.h"

#ifdef _OPENMP
#include "omp.h"
#endif



namespace Kratos
{

template<std::size_t TDim> class ParticleConvectUtily
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ParticleConvectUtily<TDim>);

    ParticleConvectUtily(typename BinBasedFastPointLocator<TDim>::Pointer pSearchStructure)
        : mpSearchStructure(pSearchStructure)
    {
    }

    ~ParticleConvectUtily()
    {
    }



    //**********************************************************************************************
    //**********************************************************************************************
    ///this function moves all the nodes contained in rModelPart from their position at time tn to the one at time
    ///tn+1 by following the trajectories. This is done by performing "subdivions" forward euler steps within each time step
    ///@param rModelPart the model part on which we work
    ///@param subdivisions number of forward euler substeps used in advancing in time
    void MoveParticles_Substepping(ModelPart& rModelPart, unsigned int subdivisions)
    {
        KRATOS_TRY
        const double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        const double small_dt = dt/ static_cast<double>(subdivisions);

        //do movement
        array_1d<double, 3 > veulerian;
        array_1d<double, 3 > acc_particle;
        array_1d<double, TDim + 1 > N;
        const int max_results = rModelPart.Nodes().size();

        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        const int nparticles = rModelPart.Nodes().size();

        #pragma omp parallel for firstprivate(results,N,veulerian,acc_particle)
        for (int i = 0; i < nparticles; i++)
        {
            unsigned int substep = 0;

            ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
            Node < 3 > ::Pointer pparticle = *(iparticle.base());

            array_1d<double,3> current_position = iparticle->GetInitialPosition() + iparticle->FastGetSolutionStepValue(DISPLACEMENT,1);

            Element::Pointer pelement;
            bool is_found = false;


            while(substep++ < subdivisions)
            {

                typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

                is_found = false;


                if(substep > 1 ) //first check if it falls within the same element
                {
                    is_found = mpSearchStructure->CalculatePosition(pelement->GetGeometry(), current_position[0], current_position[1], current_position[2], N);

                    if(is_found == false)
                        is_found = mpSearchStructure->FindPointOnMesh(current_position, N, pelement, result_begin, max_results);
                }
                else //if not found use the search structure
                {
                    is_found = mpSearchStructure->FindPointOnMesh(current_position, N, pelement, result_begin, max_results);
                }

                (iparticle)->Set(TO_ERASE, true);

                if (is_found == true)
                {
                    Geometry< Node < 3 > >& geom = pelement->GetGeometry();

                    const double new_step_factor = static_cast<double>(substep)/subdivisions;
                    const double old_step_factor = 1.0 - new_step_factor;

                    noalias(veulerian) = N[0] * ( new_step_factor*geom[0].FastGetSolutionStepValue(VELOCITY) + old_step_factor*geom[0].FastGetSolutionStepValue(VELOCITY,1));
                    for (unsigned int k = 1; k < geom.size(); k++)
                        noalias(veulerian) += N[k] * ( new_step_factor*geom[k].FastGetSolutionStepValue(VELOCITY) + old_step_factor*geom[k].FastGetSolutionStepValue(VELOCITY,1) );
                    noalias(current_position) += small_dt*veulerian;

                    (iparticle)->Set(TO_ERASE, false);

                }
                else
                    break;


            }

            if (is_found == true)
            {


                iparticle->FastGetSolutionStepValue(DISPLACEMENT) = current_position - iparticle->GetInitialPosition();
                noalias(pparticle->Coordinates()) = current_position;
            }
        }

        KRATOS_CATCH("")

    }


    //**********************************************************************************************
    //**********************************************************************************************
    ///this function moves the mesh as xn+1 = xn + vn*dt and sets the mesh velocity to vn
    ///@param rModelPart the model part on which we work
    void MoveParticles_RK4(ModelPart& rModelPart)
    {
        KRATOS_TRY
        const double dt = rModelPart.GetProcessInfo()[DELTA_TIME];

        //do movement
        array_1d<double, 3 > v1,v2,v3,v4,vtot,x;
        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        const int nparticles = rModelPart.Nodes().size();

        #pragma omp parallel for firstprivate(results,N,v1,v2,v3,v4,vtot,x)
        for (int i = 0; i < nparticles; i++)
        {
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

            ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
            Node < 3 > ::Pointer pparticle = *(iparticle.base());

            array_1d<double,3> initial_position = iparticle->GetInitialPosition() + iparticle->FastGetSolutionStepValue(DISPLACEMENT,1);

            Element::Pointer pelement;
            bool is_found = false;
            //STEP1
            {
                is_found = mpSearchStructure->FindPointOnMesh(initial_position, N, pelement, result_begin, max_results);
                if( is_found == false) goto end_of_particle;
                Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                noalias(v1) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY,1));
                for (unsigned int k = 1; k < geom.size(); k++)
                    noalias(v1) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY,1) );
            }

            //STEP2
//             if(is_found == true)
            {
                noalias(x) = initial_position + (0.5*dt)*v1;
                is_found = mpSearchStructure->FindPointOnMesh(x, N, pelement, result_begin, max_results);
                if( is_found == false) goto end_of_particle;
                Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                const double new_step_factor = 0.5;
                const double old_step_factor = 0.5;

                noalias(v2) = N[0] * ( new_step_factor*geom[0].FastGetSolutionStepValue(VELOCITY) + old_step_factor*geom[0].FastGetSolutionStepValue(VELOCITY,1));
                for (unsigned int k = 1; k < geom.size(); k++)
                    noalias(v2) += N[k] * ( new_step_factor*geom[k].FastGetSolutionStepValue(VELOCITY) + old_step_factor*geom[k].FastGetSolutionStepValue(VELOCITY,1) );
            }

            //STEP3
//             if(is_found == true)
            {
                const array_1d<double,3> x = initial_position + (0.5*dt)*v2;
                is_found = mpSearchStructure->FindPointOnMesh(x, N, pelement, result_begin, max_results);
                if( is_found == false) goto end_of_particle;
                Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                const double new_step_factor = 0.5; //as the step before
                const double old_step_factor = 0.5;

                noalias(v3) = N[0] * ( new_step_factor*geom[0].FastGetSolutionStepValue(VELOCITY) + old_step_factor*geom[0].FastGetSolutionStepValue(VELOCITY,1));
                for (unsigned int k = 1; k < geom.size(); k++)
                    noalias(v3) += N[k] * ( new_step_factor*geom[k].FastGetSolutionStepValue(VELOCITY) + old_step_factor*geom[k].FastGetSolutionStepValue(VELOCITY,1) );
            }

            //STEP4
//             if(is_found == true)
            {
                const array_1d<double,3> x = initial_position + (dt)*v3;
                is_found = mpSearchStructure->FindPointOnMesh(x, N, pelement, result_begin, max_results);
                if( is_found == false) goto end_of_particle;
                Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                noalias(v4) = N[0] * ( geom[0].FastGetSolutionStepValue(VELOCITY));
                for (unsigned int k = 1; k < geom.size(); k++)
                    noalias(v4) += N[k] * ( geom[k].FastGetSolutionStepValue(VELOCITY) );
            }


            (iparticle)->Set(TO_ERASE, false);
            //finalize step
            noalias(x) = initial_position;
            noalias(x) += 0.16666666666666666666667*dt*v1;
            noalias(x) += 0.33333333333333333333333*dt*v2;
            noalias(x) += 0.33333333333333333333333*dt*v3;
            noalias(x) += 0.16666666666666666666667*dt*v4;

            iparticle->FastGetSolutionStepValue(DISPLACEMENT) = x - iparticle->GetInitialPosition();
            noalias(pparticle->Coordinates()) = x;
            
            end_of_particle:  (iparticle)->Set(TO_ERASE, true);
        }

        KRATOS_CATCH("")

    }

//**********************************************************************************************
//**********************************************************************************************
///this function erases the elements and conditions which have at least one node marked for erase
///@param rModelPart the model part on which we work
    void EraseOuterElements(ModelPart& rModelPart)
    {
        KRATOS_TRY
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
                if( it->IsNot(TO_ERASE) )
                    (rModelPart.Elements()).push_back(*(it.base()));
            }
        }
        KRATOS_CATCH("")
    }

private:
    typename BinBasedFastPointLocator<TDim>::Pointer mpSearchStructure;



};

} // namespace Kratos.

#endif // KRATOS_CONVECT_PARTICLES_UTILITIES_INCLUDED  defined


