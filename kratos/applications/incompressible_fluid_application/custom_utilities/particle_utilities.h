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


#if !defined(KRATOS_PARTICLES_UTILITIES_INCLUDED )
#define  KRATOS_PARTICLES_UTILITIES_INCLUDED

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

template<std::size_t TDim> class ParticleUtils
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ParticleUtils<TDim>);

    //**********************************************************************************************
    //**********************************************************************************************

    void StreamlineMove(array_1d<double, 3 > & body_force, const double density, const double speficit_heat, const double dt, const double subdivisions, ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, bool use_eulerian_velocity, BinBasedFastPointLocator<TDim>& node_locator)
    {
        KRATOS_TRY


        Timer time;
        Timer::Start("StreamlineMove");


        if (rEulerianModelPart.NodesBegin()->SolutionStepsDataHas(FORCE) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----FORCE---- variable!!!!!! ERROR", "");


        //reset particle position to the beginning of the step
        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                node_it != rLagrangianModelPart.NodesEnd(); ++node_it)
        {
            Node < 3 > ::Pointer pnode = *(node_it.base());

            pnode->Set(TO_ERASE, true);
            node_it->GetValue(IS_VISITED) = 0;

        }
        //            KRATOS_WATCH("539")
        array_1d<double, 3 > veulerian;
        array_1d<double, 3 > acc_particle;
        array_1d<double, 3 > acc_particle1;
        array_1d<double, TDim + 1 > N;
        //double G;
        const int max_results = 1000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        //const double small_dt = dt / subdivisions;

        const int nparticles = rLagrangianModelPart.Nodes().size();

        #pragma omp parallel for firstprivate(results,N,veulerian,acc_particle)
        for (int i = 0; i < nparticles; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;
            double subdivisions = (iparticle)->FastGetSolutionStepValue(LAMBDA);
            subdivisions=20.0;
            const double small_dt = dt / subdivisions;
            //KRATOS_WATCH(subdivisions);
            //KRATOS_WATCH(small_dt);
            for (unsigned int substep = 0; substep < subdivisions; substep++)
            {
                //ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;

                Node < 3 > ::Pointer pparticle = *(iparticle.base());
                typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

                Element::Pointer pelement;


                bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);


                if (is_found == true)
                {
                    (pparticle)->GetValue(IS_VISITED) = 1;

                    Geometry< Node < 3 > >& geom = pelement->GetGeometry();

                    //move according to the streamline
                    noalias(veulerian) = N[0] * geom[0].FastGetSolutionStepValue(VELOCITY, 1);
                    for (unsigned int k = 1; k < geom.size(); k++)
                        noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY, 1);

                    array_1d<double, 3 > & disp = (iparticle)->FastGetSolutionStepValue(DISPLACEMENT);

                    noalias(disp) += small_dt*veulerian;

                    (pparticle)->Set(TO_ERASE, false);

                    noalias(acc_particle) = ZeroVector(3);
                    noalias(acc_particle1) = ZeroVector(3);
                    //G = 0.0;
                    for (unsigned int k = 0; k < geom.size(); k++)
                    {
                        noalias(acc_particle) += N[k] * geom[k].FastGetSolutionStepValue(FORCE) /** density_inverse*/;

                        noalias(acc_particle1) += N[k] * geom[k].FastGetSolutionStepValue(ANGULAR_ACCELERATION) /** density_inverse*/;

                    }


                    array_1d<double, 3 > & vel_particle = (pparticle)->FastGetSolutionStepValue(VELOCITY);

                    //double density_inverse_1 = 1.0 / (pparticle)->FastGetSolutionStepValue(DENSITY);

                    //double density_inverse_1 =1.0;

                    noalias(vel_particle) += small_dt*acc_particle;// * density_inverse_1;
                    //noalias(vel_particle) += small_dt*acc_particle1;//* density_inverse_1 ;


                    //update position
                    noalias(iparticle->Coordinates()) = iparticle->GetInitialPosition();
                    noalias(iparticle->Coordinates()) += iparticle->FastGetSolutionStepValue(DISPLACEMENT);
                    (iparticle)->GetValue(IS_VISITED) = 0;

                }

            }

        }

        //erase nodes whose velocity is far inconsistent with the displacement increment (typically nodes that get stuck to the wall)
        for (ModelPart::NodesContainerType::iterator it = rLagrangianModelPart.NodesBegin();
                it != rLagrangianModelPart.NodesEnd(); it++)
        {
            array_1d<double,3> delta_disp = it->FastGetSolutionStepValue(DISPLACEMENT);
            noalias(delta_disp) -= it->FastGetSolutionStepValue(DISPLACEMENT,1);

            double norm_delta_disp = norm_2(delta_disp);

            array_1d<double,3> v_old = it->FastGetSolutionStepValue(VELOCITY,1);
            double norm_v = norm_2(v_old);

            if(norm_delta_disp*3.0 < norm_v*dt )
                it->Set(TO_ERASE, true);
            if(norm_delta_disp* (0.333333333333333*0.001) >  norm_v*dt )
                it->Set(TO_ERASE, true);
        }

        //perform the erase
        //NodeEraseProcess(rLagrangianModelPart).Execute();

        Timer::Stop("StreamlineMove");
        KRATOS_WATCH(time)

        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************

    void aa(array_1d<double, 3 > & body_force, const double density, const double dt, const double subdivisions, ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, bool use_eulerian_velocity, BinBasedFastPointLocator<TDim>& node_locator)
    {
        KRATOS_TRY

        Timer time;
        Timer::Start("Actualizacion");

        if (rEulerianModelPart.NodesBegin()->SolutionStepsDataHas(FORCE) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----FORCE---- variable!!!!!! ERROR", "");

        //should be done outside!!!
        //            BinBasedFastPointLocator<TDim> node_locator(rEulerianModelPart);
        //            node_locator.UpdateSearchDatabase();

        //   double density_inverse = 1.0 / density;

        //reset particle position to the beginning of the step
        /*for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                node_it != rLagrangianModelPart.NodesEnd(); ++node_it)
        {
            Node < 3 > ::Pointer pnode = *(node_it.base());

        }*/
        //            KRATOS_WATCH("539")
        array_1d<double, 3 > veulerian;
        array_1d<double, 3 > acc_particle;
        array_1d<double, TDim + 1 > N;
        const int max_results = 1000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        //   const double small_dt = dt / subdivisions;

        const int nparticles = rLagrangianModelPart.Nodes().size();

        #pragma omp parallel for firstprivate(results,N,veulerian,acc_particle)
        for (int i = 0; i < nparticles; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;

            Node < 3 > ::Pointer pparticle = *(iparticle.base());
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

            Element::Pointer pelement;

            bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);

            if (is_found == true)
            {
                //(pparticle)->GetValue(IS_VISITED) = 1;

                Geometry< Node < 3 > >& geom = pelement->GetGeometry();


                noalias(acc_particle) = ZeroVector(3);

                //double courant = pelement->GetValue(POISSON_RATIO);

                //int Ndt_min = 3;
                //int Ndt_max = 20;
                //int K = 5;
                //int  Ndt = floor(courant *K);
                //int nNdt = (Ndt<Ndt_min) ? Ndt_min : (Ndt>Ndt_max) ? Ndt_max : Ndt;

                //double & numeros_de_pasos = (pparticle)->FastGetSolutionStepValue(LAMBDA);
                //numeros_de_pasos =nNdt;


                for (unsigned int k = 0; k < geom.size(); k++)
                {
                    acc_particle += N[k] * geom[k].FastGetSolutionStepValue(FORCE) /** density_inverse*/;

                }


                //double density_inverse = 1.0 / (pparticle)->FastGetSolutionStepValue(DENSITY);
                array_1d<double, 3 > & vel_particle = (pparticle)->FastGetSolutionStepValue(VELOCITY);


                noalias(vel_particle) += acc_particle;// * density_inverse;

            }

        }


        Timer::Stop("Actualizacion");
        KRATOS_WATCH(time)
        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************
    //function to seed a list of new nodes

    void Reseed(ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart)
    {
        KRATOS_TRY;

        unsigned int id = (rEulerianModelPart.Nodes().end() - 1)->Id() + 1;
        rLagrangianModelPart.Nodes().clear();

        /*for (ModelPart::NodesContainerType::iterator node_it = rEulerianModelPart.NodesBegin();
                  node_it != rEulerianModelPart.NodesEnd(); node_it++)
          {
              int node_id = id++;
              double x = node_it->X();
              double y = node_it->Y();
              double z = node_it->Z();
              Node < 3 > ::Pointer pnode = rLagrangianModelPart.CreateNewNode(node_id, x, y, z);

              pnode->FastGetSolutionStepValue(VELOCITY) = node_it->FastGetSolutionStepValue(VELOCITY);
          }*/

#ifdef USE_FEW_PARTICLES
        boost::numeric::ublas::bounded_matrix<double, TDim + 2, TDim + 1 > pos;
        boost::numeric::ublas::bounded_matrix<double, TDim + 2, TDim + 1 > N;
#else
        boost::numeric::ublas::bounded_matrix<double, 16, 3 > pos;
        boost::numeric::ublas::bounded_matrix<double, 16, 3 > N;
#endif
        for (ModelPart::ElementsContainerType::iterator el_it = rEulerianModelPart.ElementsBegin();
                el_it != rEulerianModelPart.ElementsEnd(); el_it++)
        {
            Geometry<Node < 3 > >& geom = el_it->GetGeometry();
            ComputeGaussPointPositions(geom, pos, N);

            for (unsigned int i = 0; i < pos.size1(); i++)
            {
                int node_id = id++;
                Node < 3 > ::Pointer pnode = rLagrangianModelPart.CreateNewNode(node_id, pos(i, 0), pos(i, 1), pos(i, 2));

                array_1d<double, 3 > & vel = pnode->FastGetSolutionStepValue(VELOCITY);
                noalias(vel) = ZeroVector(3);
                for (unsigned int j = 0; j < TDim + 1; j++)
                    noalias(vel) += N(i, j) * geom[j].FastGetSolutionStepValue(VELOCITY);
            }
        }

        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                node_it != rLagrangianModelPart.NodesEnd(); node_it++)
        {
            node_it->FastGetSolutionStepValue(VELOCITY, 1) = node_it->FastGetSolutionStepValue(VELOCITY);
        }

        KRATOS_CATCH("");
    }

    void ReseedEmptyElements(ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, BinBasedFastPointLocator<TDim>& node_locator)
    {
        KRATOS_TRY;

        //KRATOS_THROW_ERROR(std::logic_error, "Add  ----FORCE---- variable!!!!!! ERROR", "");
        Timer time;
        Timer::Start("Puntos");

        //generate a tree with the position of the lagrangian nodes
//         typedef Node < 3 > PointType;
//         typedef Node < 3 > ::Pointer PointTypePointer;
        typedef Node<3> NodeType;

        //unsigned int min_number_of_particles = 4  ;

//	  KRATOS_THROW_ERROR(std::logic_error, "Add  ----FORCE---- variable!!!!!! ERROR", "");
        int id;
        if (rLagrangianModelPart.Nodes().size() != 0)
            id = (rLagrangianModelPart.NodesEnd() - 1)->Id();
        else
            id = 1;

        for (ModelPart::ElementsContainerType::iterator el_it = rEulerianModelPart.ElementsBegin();
                el_it != rEulerianModelPart.ElementsEnd(); el_it++)
        {
            el_it->SetValue(YOUNG_MODULUS,0.0);
        }

        for (ModelPart::NodesContainerType::iterator pparticle = rLagrangianModelPart.NodesBegin(); pparticle != rLagrangianModelPart.NodesEnd(); pparticle++)
        {
            pparticle->Set(TO_ERASE,false);
        }


        int last_id=id;
        //count particles that fall within an element
        array_1d<double, TDim + 1 > N;
        const int max_results = 1000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        int nparticles = rLagrangianModelPart.Nodes().size();

        //count particles within an element
        #pragma omp parallel for firstprivate(results,N)
        for (int i = 0; i < nparticles; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;

            Node < 3 > ::Pointer pparticle = *(iparticle.base());
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

            Element::Pointer pelement;

            bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);

            if (is_found == true)
            {
                double& counter = pelement->GetValue(YOUNG_MODULUS);
                #pragma omp atomic
                counter += 1.0;
            }

        }

        //erase particles within elements for which reseeding is needed
        #pragma omp parallel for firstprivate(results,N)
        for (int i = 0; i < nparticles; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;

            Node < 3 > ::Pointer pparticle = *(iparticle.base());
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

            Element::Pointer pelement;

            bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);

            if (is_found == true)
            {
                double& counter = pelement->GetValue(YOUNG_MODULUS);
                /*KRATOS_WATCH("cantidad de particulas");
                KRATOS_WATCH("cantidad de particulas");
                KRATOS_WATCH("cantidad de particulas");
                KRATOS_WATCH(counter);*/
                double densi_ty = pelement->GetValue(DENSITY);
                if(densi_ty!=500.0 )
                {
                    if(counter  > 10.0)   //14
                    {
                        /* KRATOS_WATCH("BOOOOOOOOOOOOOOOOOOOOOOORRRRRRRRRRRRRRAAAAAAAAAAAAAARRRRRRRRRRRRRRRRRRR");
                         KRATOS_WATCH("BOOOOOOOOOOOOOOOOOOOOOOORRRRRRRRRRRRRRAAAAAAAAAAAAAARRRRRRRRRRRRRRRRRRR");
                         KRATOS_WATCH("BOOOOOOOOOOOOOOOOOOOOOOORRRRRRRRRRRRRRAAAAAAAAAAAAAARRRRRRRRRRRRRRRRRRR");	*/
                        pparticle->Set(TO_ERASE,true);
                        #pragma omp atomic
                        counter -=1.0;
                    }
                }
                else if (densi_ty==500.0)
                {
                    if(counter  > 20.0)
                    {
                        pparticle->Set(TO_ERASE,true);
                        #pragma omp atomic
                        counter -=1.0;
                    }

                }

            }

        }
        //perform the erase
        NodeEraseProcess(rLagrangianModelPart).Execute();



        int number_of_threads = OpenMPUtils::GetNumThreads();
        KRATOS_WATCH(number_of_threads);
        vector<unsigned int> elem_partition;

        int number_of_rows=rEulerianModelPart.Elements().size();
        KRATOS_WATCH(number_of_threads);
        //KRATOS_THROW_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");

        elem_partition.resize(number_of_threads + 1);
        int elem_partition_size = number_of_rows / number_of_threads;
        elem_partition[0] = 0;
        elem_partition[number_of_threads] = number_of_rows;
        KRATOS_WATCH(elem_partition_size);
        for ( int i = 1; i < number_of_threads; i++)
            elem_partition[i] = elem_partition[i - 1] + elem_partition_size;

//         typedef Node < 3 > PointType;
        std::vector<ModelPart::NodesContainerType> aux;// aux;
        aux.resize(number_of_threads);


        #pragma omp parallel firstprivate(elem_partition)
        {
            int k = OpenMPUtils::ThisThread();
            ModelPart::ElementsContainerType::iterator it_begin = rEulerianModelPart.ElementsBegin() +  elem_partition[k];
            ModelPart::ElementsContainerType::iterator it_end = rEulerianModelPart.ElementsBegin() + elem_partition[k+1] ;
            //KRATOS_WATCH(min_number_of_particles);
            //ModelPart::NodesContainerType local_list=aux[k];
            PointerVectorSet<Node<3>, IndexedObject> & list=aux[k];
            //KRATOS_WATCH(k);

            boost::numeric::ublas::bounded_matrix<double, 4, 3 > pos;
            boost::numeric::ublas::bounded_matrix<double, 4, 3 > Nnew;



            int local_id=1;
            for (ModelPart::ElementsContainerType::iterator el_it = it_begin; el_it != it_end; el_it++)
            {

                if (el_it->GetValue(YOUNG_MODULUS) < 4.0)
//		if (el_it->GetValue(YOUNG_MODULUS) < 1.0)
                {
                    //KRATOS_THROW_ERROR(std::logic_error, "Add  ----FORCE---- variable!!!!!! ERROR", "");
                    Geometry< Node<3> >& geom = el_it->GetGeometry();

                    //if(el_it->GetValue(DENSITY)==1.0 or el_it->GetValue(DENSITY)==1000.0 ){

                    //KRATOS_WATCH("INSERTO PUNTOS");
                    //KRATOS_WATCH("INSERTO PUNTOS");
                    //KRATOS_WATCH("INSERTO PUNTOS");
//		    double dens=el_it->GetValue(DENSITY);
                    ComputeGaussPointPositions(geom, pos, Nnew);

                    for (unsigned int i = 0; i < pos.size1(); i++)
                    {
                        int node_id = local_id++;

                        Node < 3 > ::Pointer pnode = NodeType::Pointer(new NodeType(node_id, pos(i,0), pos(i,1), pos(i,2)));

                        pnode->SetSolutionStepVariablesList(&(rEulerianModelPart.GetNodalSolutionStepVariablesList()));
                        pnode->SetBufferSize(3);

                        array_1d<double, 3 > & vel = pnode->FastGetSolutionStepValue(VELOCITY);

                        array_1d<double, 3 > & disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
                        noalias(disp) = ZeroVector(3);

                        noalias(vel) = ZeroVector(3);

                        for (unsigned int j = 0; j < 3; j++)
                        {

                            noalias(vel) += Nnew(i, j) * geom[j].FastGetSolutionStepValue(VELOCITY);

                        }
                        pnode->GetSolutionStepValue(FLAG_VARIABLE)=1.0;
                        //pnode->FastGetSolutionStepValue(DENSITY)=dens;
                        (list).push_back(pnode);
                    }


                }

            }
        }

        for(int k=0; k<number_of_threads; k++)
        {
            PointerVectorSet<Node<3>, IndexedObject> & list=aux[k];
            for(PointerVectorSet<Node<3>, IndexedObject>::ptr_iterator it =  list.ptr_begin(); it!=list.ptr_end(); it++)
            {
                //rLagrangianModelPart.AddNode(*it);
                rLagrangianModelPart.Nodes().push_back(*it);
                //KRATOS_WATCH(k);
                int node_id=last_id++;//last_id++;
                (*it)->SetId(node_id);
                //KRATOS_WATCH(rLagrangianModelPart.Nodes().size());
            }
        }


        Timer::Stop("Puntos");
        KRATOS_WATCH(time);

        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************
    //**********************************************************************************************
    void Back(array_1d<double, 3 > & body_force, const double density, const double dt, const double subdivisions, ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, bool use_eulerian_velocity, BinBasedFastPointLocator<TDim>& node_locator)
    {

        KRATOS_TRY

        if (rEulerianModelPart.NodesBegin()->SolutionStepsDataHas(FORCE) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----FORCE---- variable!!!!!! ERROR", "");

        //should be done outside!!!
        //            BinBasedFastPointLocator<TDim> node_locator(rEulerianModelPart);
        //            node_locator.UpdateSearchDatabase();

//	  double density_inverse = 1.0 / density;

        //reset particle position to the beginning of the step
        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin(); node_it != rLagrangianModelPart.NodesEnd(); ++node_it)
        {
            Node < 3 > ::Pointer pnode = *(node_it.base());

            //pnode->Set(TO_ERASE, true);
            node_it->GetValue(IS_VISITED) = 0;

        }
        //            KRATOS_WATCH("539")
        array_1d<double, 3 > veulerian;
        array_1d<double, 3 > acc_particle;
        array_1d<double, TDim + 1 > N;
        const int max_results = 1000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        //subdivisions=10.0;
        const double small_dt = dt / subdivisions;

        const int nparticles = rLagrangianModelPart.Nodes().size();

        //#pragma omp parallel for firstprivate(results,N,veulerian,acc_particle)
        for (int i = 0; i < nparticles; i++)
        {
            for (unsigned int substep = 0; substep < subdivisions; substep++)
            {
                ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;

                Node < 3 > ::Pointer pparticle = *(iparticle.base());

                if((iparticle)->FastGetSolutionStepValue(FLAG_VARIABLE)==1.0)
                {

                    typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

                    Element::Pointer pelement;

                    bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);
                    if(is_found==false)
                    {
                        pparticle->Set(TO_ERASE,true);
                    }
                    if (is_found == true)
                    {
                        (pparticle)->GetValue(IS_VISITED) = 1;

                        Geometry< Node < 3 > >& geom = pelement->GetGeometry();

                        //move according to the streamline
                        noalias(veulerian) = N[0] * geom[0].FastGetSolutionStepValue(VELOCITY);
                        for (unsigned int k = 1; k < geom.size(); k++)
                            noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY);
                        //KRATOS_WATCH(veulerian);
                        array_1d<double, 3 > & disp = (iparticle)->FastGetSolutionStepValue(DISPLACEMENT);

                        noalias(disp) -= small_dt*veulerian;
                        //KRATOS_WATCH(disp);

                        ///////////////
                        array_1d<double, 3 > & vel_particle = (pparticle)->FastGetSolutionStepValue(VELOCITY);

                        noalias(vel_particle) = veulerian;
                        ////////////////////


                        //(pparticle)->Set(TO_ERASE, false);
                        //double pp=(iparticle)->X();
                        //double pp1=(iparticle)->Y();
                        //double pp2=(iparticle)->Z();
                        //KRATOS_WATCH(pp);// + old_disp[0];
                        //KRATOS_WATCH(pp1);// + old_disp[1];
                        //KRATOS_WATCH(pp2);// + old_disp[2];

                        //update position
                        noalias(iparticle->Coordinates()) = iparticle->GetInitialPosition();
                        noalias(iparticle->Coordinates()) += iparticle->FastGetSolutionStepValue(DISPLACEMENT);
                        (iparticle)->GetValue(IS_VISITED) = 0;

                        //double pp_1=(iparticle)->X();
                        //double pp1_1=(iparticle)->Y();
                        //double pp2_1=(iparticle)->Z();
                        //KRATOS_WATCH(pp_1);// + old_disp[0];
                        //KRATOS_WATCH(pp1_1);// + old_disp[1];
                        //KRATOS_WATCH(pp2_1);// + old_disp[2];


                        //KRATOS_WATCH((iparticle)->X());// + old_disp[0];
                        //KRATOS_WATCH((iparticle)->Y());// + old_disp[1];
                        //KRATOS_WATCH((iparticle)->Z());// + old_disp[2];


                    }
                }
            }



        }

        // NodeEraseProcess(rLagrangianModelPart).Execute();

        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    void Back1(array_1d<double, 3 > & body_force, const double density, const double dt, const double subdivisions, ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, bool use_eulerian_velocity, BinBasedFastPointLocator<TDim>& node_locator)
    {

        KRATOS_TRY

        if (rEulerianModelPart.NodesBegin()->SolutionStepsDataHas(FORCE) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----FORCE---- variable!!!!!! ERROR", "");

        //should be done outside!!!
        //            BinBasedFastPointLocator<TDim> node_locator(rEulerianModelPart);
        //            node_locator.UpdateSearchDatabase();

//	  double density_inverse = 1.0 / density;

        //reset particle position to the beginning of the step
        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin(); node_it != rLagrangianModelPart.NodesEnd(); ++node_it)
        {
            Node < 3 > ::Pointer pnode = *(node_it.base());

            //pnode->Set(TO_ERASE, true);
            node_it->GetValue(IS_VISITED) = 0;

        }
        //            KRATOS_WATCH("539")
        array_1d<double, 3 > veulerian;
        array_1d<double, 3 > acc_particle;
        array_1d<double, TDim + 1 > N;
        const int max_results = 1000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        //subdivisions=10.0;
        const double small_dt = dt / subdivisions;

        const int nparticles = rLagrangianModelPart.Nodes().size();

//#pragma omp parallel for firstprivate(results,N,veulerian,acc_particle)
        for (int i = 0; i < nparticles; i++)
        {
            for (unsigned int substep = 0; substep < subdivisions; substep++)
            {
                ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;

                Node < 3 > ::Pointer pparticle = *(iparticle.base());

                if((iparticle)->FastGetSolutionStepValue(FLAG_VARIABLE)==1.0)
                {

                    typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

                    Element::Pointer pelement;

                    bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);
                    if(is_found==false)
                    {
                        pparticle->Set(TO_ERASE,true);
                    }
                    if (is_found == true)
                    {
                        (pparticle)->GetValue(IS_VISITED) = 1;

                        Geometry< Node < 3 > >& geom = pelement->GetGeometry();

                        //move according to the streamline
                        noalias(veulerian) = N[0] * geom[0].FastGetSolutionStepValue(VELOCITY);
                        for (unsigned int k = 1; k < geom.size(); k++)
                            noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY);
                        //KRATOS_WATCH(veulerian);
                        array_1d<double, 3 > & disp = (iparticle)->FastGetSolutionStepValue(DISPLACEMENT);




                        //KRATOS_WATCH((iparticle)->X());// + old_disp[0];
                        //KRATOS_WATCH((iparticle)->Y());// + old_disp[1];
                        //KRATOS_WATCH((iparticle)->Z());// + old_disp[2];


                        noalias(disp) += small_dt*veulerian;
                        //KRATOS_WATCH(disp);

                        ///////////////
                        array_1d<double, 3 > & vel_particle = (pparticle)->FastGetSolutionStepValue(VELOCITY);

                        noalias(vel_particle) = veulerian;
                        ////////////////////


                        //(pparticle)->Set(TO_ERASE, false);
                        //double pp=(iparticle)->X();
                        //double pp1=(iparticle)->Y();
                        //double pp2=(iparticle)->Z();
                        //KRATOS_WATCH(pp);// + old_disp[0];
                        //KRATOS_WATCH(pp1);// + old_disp[1];
                        //KRATOS_WATCH(pp2);// + old_disp[2];


                        //update position
                        noalias(iparticle->Coordinates()) = iparticle->GetInitialPosition();
                        noalias(iparticle->Coordinates()) += iparticle->FastGetSolutionStepValue(DISPLACEMENT);
                        (iparticle)->GetValue(IS_VISITED) = 0;

                        //double pp_1=(iparticle)->X();
                        //double pp1_1=(iparticle)->Y();
                        //double pp2_1=(iparticle)->Z();
                        //KRATOS_WATCH(pp_1);// + old_disp[0];
                        //KRATOS_WATCH(pp1_1);// + old_disp[1];
                        //KRATOS_WATCH(pp2_1);// + old_disp[2];
                        //KRATOS_WATCH((iparticle)->X());// + old_disp[0];
                        //KRATOS_WATCH((iparticle)->Y());// + old_disp[1];
                        //KRATOS_WATCH((iparticle)->Z());// + old_disp[2];


                    }
                }
            }

        }


        // NodeEraseProcess(rLagrangianModelPart).Execute();
        KRATOS_CATCH("")
    }

    //**********************************************************************************************
    //**********************************************************************************************
    void StreamlineMove2(array_1d<double, 3 > & body_force, const double density, const double speficit_heat, const double dt, const double subdivisions, ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, bool use_eulerian_velocity, BinBasedFastPointLocator<TDim>& node_locator)
    {
        KRATOS_TRY


        Timer time;
        Timer::Start("StreamlineMove");


        if (rEulerianModelPart.NodesBegin()->SolutionStepsDataHas(FORCE) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----FORCE---- variable!!!!!! ERROR", "");


        //reset particle position to the beginning of the step
        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                node_it != rLagrangianModelPart.NodesEnd(); ++node_it)
        {
            Node < 3 > ::Pointer pnode = *(node_it.base());
            pnode->Set(TO_ERASE, true);
            node_it->GetValue(IS_VISITED) = 0;

        }
        //            KRATOS_WATCH("539")
        array_1d<double, 3 > veulerian;
        array_1d<double, 3 > acc_particle;
        array_1d<double, 3 > acc_particle1;
        array_1d<double, TDim + 1 > N;
        //double G;
        const int max_results = 1000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        //const double small_dt = dt / subdivisions;

        const int nparticles = rLagrangianModelPart.Nodes().size();

        #pragma omp parallel for firstprivate(results,N,veulerian,acc_particle)
        for (int i = 0; i < nparticles; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;
            double subdivisions = (iparticle)->FastGetSolutionStepValue(LAMBDA);
            subdivisions=10.0;
            const double small_dt = dt / subdivisions;

            for (unsigned int substep = 0; substep < subdivisions; substep++)
            {

                Node < 3 > ::Pointer pparticle = *(iparticle.base());

                if((iparticle)->FastGetSolutionStepValue(FLAG_VARIABLE)==1.0)
                {

                    typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

                    Element::Pointer pelement;


                    bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);


                    if (is_found == true)
                    {
                        (pparticle)->GetValue(IS_VISITED) = 1;

                        Geometry< Node < 3 > >& geom = pelement->GetGeometry();

                        //move according to the streamline
                        noalias(veulerian) = N[0] * geom[0].FastGetSolutionStepValue(VELOCITY, 1);
                        for (unsigned int k = 1; k < geom.size(); k++)
                            noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY, 1);

                        array_1d<double, 3 > & disp = (iparticle)->FastGetSolutionStepValue(DISPLACEMENT);

                        noalias(disp) += small_dt*veulerian;

                        (pparticle)->Set(TO_ERASE, false);

                        noalias(acc_particle) = ZeroVector(3);
                        noalias(acc_particle1) = ZeroVector(3);
                        //G = 0.0;
                        for (unsigned int k = 0; k < geom.size(); k++)
                        {
                            noalias(acc_particle) += N[k] * geom[k].FastGetSolutionStepValue(FORCE) /** density_inverse*/;

                            //  noalias(acc_particle1) += N[k] * geom[k].FastGetSolutionStepValue(ANGULAR_ACCELERATION) /** density_inverse*/;

                        }


                        array_1d<double, 3 > & vel_particle = (pparticle)->FastGetSolutionStepValue(VELOCITY);

//			  double density_inverse_1 = 1.0 / (pparticle)->FastGetSolutionStepValue(DENSITY);

                        //double density_inverse_1 =1.0;

                        noalias(vel_particle) += small_dt*acc_particle;// * density_inverse_1;
                        //noalias(vel_particle) += small_dt*acc_particle1;//* density_inverse_1 ;


                        //update position
                        noalias(iparticle->Coordinates()) = iparticle->GetInitialPosition();
                        noalias(iparticle->Coordinates()) += iparticle->FastGetSolutionStepValue(DISPLACEMENT);
                        (iparticle)->GetValue(IS_VISITED) = 0;
                    }
                }

            }

        }

        //erase nodes whose velocity is far inconsistent with the displacement increment (typically nodes that get stuck to the wall)


        Timer::Stop("StreamlineMove");
        KRATOS_WATCH(time)

        KRATOS_CATCH("")
    }
    //**********************************************************************************************

    //function to seed a list of new nodes

    void Density(ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, BinBasedFastPointLocator<TDim>& node_locator)
    {
        KRATOS_TRY;

        Timer time;
        Timer::Start("Puntos");

        //generate a tree with the position of the lagrangian nodes
//         typedef Node < 3 > PointType;
//         typedef Node < 3 > ::Pointer PointTypePointer;
//         typedef Node<3> NodeType;

        //unsigned int min_number_of_particles = 4  ;

//        int id;
//        if (rLagrangianModelPart.Nodes().size() != 0)
//            id = (rLagrangianModelPart.NodesEnd() - 1)->Id();
//        else
//            id = 1;

        for (ModelPart::ElementsContainerType::iterator el_it = rEulerianModelPart.ElementsBegin();
                el_it != rEulerianModelPart.ElementsEnd(); el_it++)
        {
            //                el_it->SetValue(DENSITY,0.0);
            el_it->SetValue(YOUNG_MODULUS,0.0);
            el_it->SetValue(POISSON_RATIO,0.0);
        }

        //int last_id=id;
        //count particles that fall within an element
        array_1d<double, TDim + 1 > N;
        const int max_results = 1000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        int nparticles = rLagrangianModelPart.Nodes().size();
        //double density=0.0;
        //int number=0.0;

        //count particles within an element
        //#pragma omp parallel for firstprivate(results,N)
        for (int i = 0; i < nparticles; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;

            Node < 3 > ::Pointer pparticle = *(iparticle.base());
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

            Element::Pointer pelement;

            bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);

            if (is_found == true)
            {
                double density = (pparticle)->FastGetSolutionStepValue(DENSITY);
                double& counter_min = pelement->GetValue(YOUNG_MODULUS);
                double& counter_max = pelement->GetValue(POISSON_RATIO);

                if(density==1.0)
                {
                    counter_min +=1.0;
                }

                if(density==1000.0)
                {
                    counter_max +=1.0;
                }
            }
        }


        for (ModelPart::ElementsContainerType::iterator el_it = rEulerianModelPart.ElementsBegin(); el_it != rEulerianModelPart.ElementsEnd(); el_it++)
        {
            double& dens = el_it->GetValue(DENSITY);
            double counter_min = el_it->GetValue(YOUNG_MODULUS);
            double counter_max = el_it->GetValue(POISSON_RATIO);
            if(counter_min==0.0 && counter_max>0.0)
            {
                dens=1000.0;
            }

            if(counter_min>0.0 && counter_max==0.0)
            {
                dens=1.0;
            }
            if(counter_min>0.0 && counter_max>0.0)
            {
                dens=500.0;
            }

            if(counter_min==0.0 &&  counter_max==0.0)
            {
                //dens=1.0;
                KRATOS_WATCH("SIN PARTICULAS");
                KRATOS_WATCH("SIN PARTICULAS");
                KRATOS_WATCH("SIN PARTICULAS");
                KRATOS_WATCH("SIN PARTICULAS");
                KRATOS_WATCH(dens);
                KRATOS_THROW_ERROR(std::logic_error,  "NEGATIVE VALUE OF Time step estimated" , "");
            }

            //KRATOS_WATCH(dens);

        }



        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    void Density1(ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, BinBasedFastPointLocator<TDim>& node_locator)
    {
        KRATOS_TRY;

        Timer time;
        Timer::Start("Puntos");

        //generate a tree with the position of the lagrangian nodes
//         typedef Node < 3 > PointType;
//         typedef Node < 3 > ::Pointer PointTypePointer;
//         typedef Node<3> NodeType;

        //unsigned int min_number_of_particles = 4  ;

//        int id;

//        if (rLagrangianModelPart.Nodes().size() != 0)
//            id = (rLagrangianModelPart.NodesEnd() - 1)->Id();
//        else
//            id = 1;

        //int last_id=id;
        //count particles that fall within an element
        array_1d<double, TDim + 1 > N;
        const int max_results = 1000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        int nparticles = rLagrangianModelPart.Nodes().size();

        //count particles within an element
        //#pragma omp parallel for firstprivate(results,N)
        for (int i = 0; i < nparticles; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;

            Node < 3 > ::Pointer pparticle = *(iparticle.base());

            if((iparticle)->FastGetSolutionStepValue(FLAG_VARIABLE)==1.0)
            {

                typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

                Element::Pointer pelement;

                bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);

                if (is_found == true)
                {

                    Geometry< Node < 3 > >& geom = pelement->GetGeometry();

                    double& density = (pparticle)->FastGetSolutionStepValue(DENSITY);

                    double density_element = pelement->GetValue(DENSITY);
                    double density_node;


                    if(N(0)>0.5)
                    {
                        density_node= geom[0].FastGetSolutionStepValue(DENSITY);
                        if(density_node<600.0) density=1.0;
                        else
                            density=1000.0;

                    }
                    else if(N(1)>0.5)
                    {
                        density_node= geom[1].FastGetSolutionStepValue(DENSITY);
                        if(density_node<600.0) density=1.0;
                        else
                            density=1000.0;

                    }
                    else if(N(2)>0.5)
                    {

                        density_node= geom[2].FastGetSolutionStepValue(DENSITY);
                        if(density_node<600.0) density=1.0;
                        else
                            density=1000.0;

                    }
                    else
                    {
                        if(density_element>600.0)
                        {
                            density=1000.0;
                        }
                        else if(density_element<=600.0)
                        {
                            density=1.0;
                        }

                    }

                    /*if( density_element ==1.00 )density=1.0;
                    else if(density_element >1.00 and density_element<500){density=1.0; pparticle->Set(TO_ERASE, true);}
                    else if( density_element >=500.00 and density_element<1000) {density=1000.0;  pparticle->Set(TO_ERASE, true);}
                    else if(density_element=1000.0) density=1000.0;*/
                }
            }
        }

        //NodeEraseProcess(rLagrangianModelPart).Execute();
        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************
    void EstimateTime(ModelPart& rEulerianModelPart,const double max_dt)
    {
        KRATOS_TRY
        // KRATOS_THROW_ERROR(std::logic_error,  "NEGATIVE VALUE OF Time step estimated" , "");
        //initializee dt with max dt
        //initialize dt with incredible value
        double /*dt, glob_min_dt,*/ dummy;
// 	    double h, nu;
        array_1d<double,3> N = ZeroVector(3);
        array_1d<double,3> aux = ZeroVector(3); //dimension = number of nodes
        array_1d<double,3> vel = ZeroVector(3); //dimension = number of nodes
        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
        array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension

        //initialize it with given value
//        glob_min_dt=max_dt;


//        dt=0.0;
        for(ModelPart::ElementsContainerType::iterator im = rEulerianModelPart.ElementsBegin() ; im !=rEulerianModelPart.ElementsEnd() ; ++im)
        {
            GeometryUtils::CalculateGeometryData(im->GetGeometry(),DN_DX,N,dummy);

            double h = sqrt(2.00*dummy);

            array_1d<double,3> const& v = im->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            ms_vel_gauss[0] = v[0];
            ms_vel_gauss[1] = v[1];

            //direction of the height is stored in the auxilliary vector
            for (unsigned int i=1; i<3; i++)
            {
                array_1d<double,3> const& vi = im->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
                ms_vel_gauss[0] += vi[0];
                ms_vel_gauss[1] += vi[1];
            }
            ms_vel_gauss *=0.3333;

            double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
            norm_u = sqrt(norm_u);

            double courant= norm_u * max_dt / h;


            double& counter = im->GetValue(POISSON_RATIO);
            counter = courant;

        }


        KRATOS_CATCH("");
    }
    //**********************************************************************************************
    //**********************************************************************************************


    void VisualizationModelPart(ModelPart& rCompleteModelPart, ModelPart& rEulerianModelPart, ModelPart & rLagrangianModelPart)
    {
        KRATOS_TRY;

        rCompleteModelPart.Elements() = rEulerianModelPart.Elements();
        rCompleteModelPart.Nodes() = rEulerianModelPart.Nodes();


        unsigned int id;
        if(rEulerianModelPart.Nodes().size()!= 0)
            id = (rEulerianModelPart.Nodes().end() - 1)->Id() + 1;
        else
            id = 1;

        //preallocate the memory needed
        int tot_nodes = rEulerianModelPart.Nodes().size() + rLagrangianModelPart.Nodes().size();
        rCompleteModelPart.Nodes().reserve( tot_nodes );

        //note that here we renumber the nodes
        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                node_it != rLagrangianModelPart.NodesEnd(); node_it++)
        {
			node_it->SetId(id++);
            rCompleteModelPart.AddNode(*(node_it.base()));
        }

        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************

    void TransferToEulerianMesh(ModelPart& rEulerianModelPart, ModelPart & rLagrangianModelPart)
    {
        KRATOS_TRY

        //defintions for spatial search
        typedef Node < 3 > PointType;
        typedef Node < 3 > ::Pointer PointTypePointer;
        typedef std::vector<PointType::Pointer> PointVector;
        typedef std::vector<PointType::Pointer>::iterator PointIterator;
        typedef std::vector<double> DistanceVector;
        typedef std::vector<double>::iterator DistanceIterator;

        //creating an auxiliary list for the new nodes
        PointVector list_of_nodes;

        //*************
        // Bucket types
        typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
        // 			   typedef Bins< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
        // 			   typedef BinsDynamic< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > DynamicBins;
        //*************
        // DynamicBins;

        typedef Tree< KDTreePartition<BucketType> > tree; //Kdtree;
        // 			   typedef Tree< OCTreePartition<BucketType> > tree; 		//Octree;
        // 			   typedef Tree< StaticBins > tree; 		     		//Binstree;
        // 			   typedef Tree< KDTreePartition<StaticBins> > tree; 		//KdtreeBins;
        // 			   typedef typename KdtreeBins::Partitions SubPartitions;
        // 			   typedef Tree< OCTreePartition<StaticBins> > tree; 		//OctreeBins;
        /*
                                   typedef Bins< TDim, PointType, stdPointVector> stdBins;
                                   typedef Tree< Bins<TDim,PointType,stdPointVector> > tree; 	//stdStaticBins;*/

        //starting calculating time of construction of the kdtree
        boost::timer kdtree_construction;

        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                node_it != rLagrangianModelPart.NodesEnd(); ++node_it)
        {
            PointTypePointer pnode = *(node_it.base());

            //putting the nodes of the destination_model part in an auxiliary list
            list_of_nodes.push_back(pnode);
        }

        std::cout << "kdt constructin time " << kdtree_construction.elapsed() << std::endl;

        //create a spatial database with the list of new nodes
        unsigned int bucket_size = 20;
        tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);

        //work arrays
        Node < 3 > work_point(0, 0.0, 0.0, 0.0);
        unsigned int MaximumNumberOfResults = 10000;
        PointVector Results(MaximumNumberOfResults);
        DistanceVector SquaredResultsDistances(MaximumNumberOfResults);


        if (rEulerianModelPart.NodesBegin()->SolutionStepsDataHas(NODAL_H) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");

        double sigma = 0.0;
        if (TDim == 2)
            sigma = 10.0 / (7.0 * 3.1415926);
        else
            sigma = 1.0 / 3.1415926;

        for (ModelPart::NodesContainerType::iterator node_it = rEulerianModelPart.NodesBegin();
                node_it != rEulerianModelPart.NodesEnd(); node_it++)
        {
            work_point.X() = node_it->X();
            work_point.Y() = node_it->Y();
            work_point.Z() = node_it->Z();

            double radius = 0.6 * node_it->FastGetSolutionStepValue(NODAL_H);

            //find all of the new nodes within the radius
            int number_of_points_in_radius;

            //look between the new nodes which of them is inside the radius of the circumscribed cyrcle
            number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),
                                         SquaredResultsDistances.begin(), MaximumNumberOfResults);


            //PRUEBA
            array_1d<double, 3 > & vel1 = (node_it)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > original_vel1 = vel1;

            if (number_of_points_in_radius > 0)
            {

                array_1d<double, 3 > & vel = (node_it)->FastGetSolutionStepValue(VELOCITY);
                //double& temperature = (node_it)->FastGetSolutionStepValue(TEMPERATURE);

                array_1d<double, 3 > original_vel = vel;
                //double original_temperature = temperature;

                noalias(vel) = ZeroVector(3);
                //temperature = 0.0;

                double tot_weight = 0.0;

                for (int k = 0; k < number_of_points_in_radius; k++)
                {
                    //                           double weight = 1.0;
                    double distance = sqrt(*(SquaredResultsDistances.begin() + k));
                    double weight = SPHCubicKernel(sigma, distance, radius);

                    tot_weight += weight;

                    PointIterator it_found = Results.begin() + k;

                    //                            array_1d<double,3> aux = (*it_found)->Coordinates()-node_it->Coordinates();
                    //                            KRATOS_WATCH(norm_2(aux));
                    //                            KRATOS_WATCH( *(SquaredResultsDistances.begin()+k) );

                    const array_1d<double, 3 > particle_velocity = (*it_found)->FastGetSolutionStepValue(VELOCITY);
                    //const double particle_temperature = (*it_found)->FastGetSolutionStepValue(TEMPERATURE);

                    noalias(vel) += weight * particle_velocity;
                    //temperature += weight * particle_temperature;
                }

                vel /= tot_weight;
                //temperature /= tot_weight;

                if (node_it->IsFixed(VELOCITY_X))
                {
                    noalias(vel) = original_vel;
                }
                /*if (node_it->IsFixed(TEMPERATURE))
                    temperature = original_temperature;*/
            }
            else
            {
                if (node_it->IsFixed(VELOCITY_X))
                {
                    noalias(vel1) = original_vel1;
                }
            }
        }
        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************

    void TransferToEulerianMeshShapeBased(ModelPart& rEulerianModelPart, ModelPart & rLagrangianModelPart, BinBasedFastPointLocator<TDim>& node_locator)
    {
        KRATOS_TRY

        Timer time;
        Timer::Start("Interpolacion");

        if (rEulerianModelPart.NodesBegin()->SolutionStepsDataHas(FORCE) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----FORCE---- variable!!!!!! ERROR", "");
        if (rEulerianModelPart.NodesBegin()->SolutionStepsDataHas(TEMPERATURE) == false)
            KRATOS_THROW_ERROR(std::logic_error, "Add  ----TEMPERATURE---- variable!!!!!! ERROR", "");

        //defintions for spatial search
//         typedef Node < 3 > PointType;
//         typedef Node < 3 > ::Pointer PointTypePointer;


        for (ModelPart::NodesContainerType::iterator node_it = rEulerianModelPart.NodesBegin();
                node_it != rEulerianModelPart.NodesEnd(); node_it++)
        {
            (node_it)->GetValue(POISSON_RATIO) = 0.0;
            (node_it)->FastGetSolutionStepValue(DENSITY) = 0.0;
            //(node_it)->FastGetSolutionStepValue(FLAG_VARIABLE) = 0.0;
            if (node_it->IsFixed(VELOCITY_X) == false)
            {
                (node_it)->FastGetSolutionStepValue(VELOCITY) = ZeroVector(3);
                (node_it)->GetValue(YOUNG_MODULUS) = 0.0;
            }
        }

        array_1d<double, TDim + 1 > N;
        const int max_results = 1000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        const int nparticles = rLagrangianModelPart.Nodes().size();

        #pragma omp parallel for firstprivate(results,N)
        for (int i = 0; i < nparticles; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;

            Node < 3 > ::Pointer pparticle = *(iparticle.base());
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

            Element::Pointer pelement;

            bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);

            if (is_found == true)
            {
                Geometry<Node<3> >& geom = pelement->GetGeometry();
                const array_1d<double, 3 > & vel_particle = (iparticle)->FastGetSolutionStepValue(VELOCITY);
                double density_particle = (iparticle)->FastGetSolutionStepValue(DENSITY);
                //double flag_variable_particle = (iparticle)->FastGetSolutionStepValue(FLAG_VARIABLE);
                for (unsigned int k = 0; k < geom.size(); k++)
                {

                    geom[k].SetLock();
                    geom[k].FastGetSolutionStepValue(DENSITY) += N[k] * density_particle;
                    geom[k].GetValue(POISSON_RATIO) += N[k];
                    geom[k].UnSetLock();


                    if (geom[k].IsFixed(VELOCITY_X) == false)
                    {
                        geom[k].SetLock();
                        geom[k].FastGetSolutionStepValue(VELOCITY) += N[k] * vel_particle;
                        geom[k].GetValue(YOUNG_MODULUS) += N[k];
                        geom[k].UnSetLock();
                    }

                }
            }
        }


        for (ModelPart::NodesContainerType::iterator node_it = rEulerianModelPart.NodesBegin();
                node_it != rEulerianModelPart.NodesEnd(); node_it++)
        {

            const double NN1 = (node_it)->GetValue(POISSON_RATIO);
            if (NN1 != 0.0)
            {
                //double pepe;
                (node_it)->FastGetSolutionStepValue(DENSITY) /= NN1;
                //(node_it)->FastGetSolutionStepValue(FLAG_VARIABLE) /= NN1;
            }
            else
            {
                std::cout << node_it->Id() << " coeff = " << NN1 << std::endl;
            }



            if (node_it->IsFixed(VELOCITY_X) == false)
            {
                const double NN = (node_it)->GetValue(YOUNG_MODULUS);
                if (NN != 0.0)
                {
                    (node_it)->FastGetSolutionStepValue(VELOCITY) /= NN;
                }
                else
                {
                    std::cout << node_it->Id() << " coeff = " << NN << std::endl;
                }
            }

        }


        Timer::Stop("Interpolacion");
        //KRATOS_WATCH(time)

        KRATOS_CATCH("")
    }

    //restarting the step from the beginning

    void RestartStep(ModelPart & rModelPart)
    {
        KRATOS_TRY;

        //setting the variables to their value at the beginning of the time step
        rModelPart.OverwriteSolutionStepData(1, 0);

        //setting the coordinates to their value at the beginning of the step
        for (ModelPart::NodesContainerType::iterator node_it = rModelPart.NodesBegin();
                node_it != rModelPart.NodesEnd(); node_it++)
        {
            array_1d<double, 3 > & coords = node_it->Coordinates();
            const array_1d<double, 3 > & old_disp = node_it->FastGetSolutionStepValue(DISPLACEMENT, 1);

            coords[0] = node_it->X0() + old_disp[0];
            coords[1] = node_it->Y0() + old_disp[1];
            coords[2] = node_it->Z0() + old_disp[2];
        }


        KRATOS_CATCH("");
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

    inline void CalculateCenterAndSearchRadius(Geometry<Node < 3 > >&geom,
            double& xc, double& yc, double& zc, double& R, array_1d<double, 3 > & N
                                              )
    {
        double x0 = geom[0].X();
        double y0 = geom[0].Y();
        double x1 = geom[1].X();
        double y1 = geom[1].Y();
        double x2 = geom[2].X();
        double y2 = geom[2].Y();


        xc = 0.3333333333333333333 * (x0 + x1 + x2);
        yc = 0.3333333333333333333 * (y0 + y1 + y2);
        zc = 0.0;

        double R1 = (xc - x0)*(xc - x0) + (yc - y0)*(yc - y0);
        double R2 = (xc - x1)*(xc - x1) + (yc - y1)*(yc - y1);
        double R3 = (xc - x2)*(xc - x2) + (yc - y2)*(yc - y2);

        R = R1;
        if (R2 > R) R = R2;
        if (R3 > R) R = R3;

        R = 1.01 * sqrt(R);
    }
    //***************************************
    //***************************************

    inline void CalculateCenterAndSearchRadius(Geometry<Node < 3 > >&geom,
            double& xc, double& yc, double& zc, double& R, array_1d<double, 4 > & N

                                              )
    {
        double x0 = geom[0].X();
        double y0 = geom[0].Y();
        double z0 = geom[0].Z();
        double x1 = geom[1].X();
        double y1 = geom[1].Y();
        double z1 = geom[1].Z();
        double x2 = geom[2].X();
        double y2 = geom[2].Y();
        double z2 = geom[2].Z();
        double x3 = geom[3].X();
        double y3 = geom[3].Y();
        double z3 = geom[3].Z();


        xc = 0.25 * (x0 + x1 + x2 + x3);
        yc = 0.25 * (y0 + y1 + y2 + y3);
        zc = 0.25 * (z0 + z1 + z2 + z3);

        double R1 = (xc - x0)*(xc - x0) + (yc - y0)*(yc - y0) + (zc - z0)*(zc - z0);
        double R2 = (xc - x1)*(xc - x1) + (yc - y1)*(yc - y1) + (zc - z1)*(zc - z1);
        double R3 = (xc - x2)*(xc - x2) + (yc - y2)*(yc - y2) + (zc - z2)*(zc - z2);
        double R4 = (xc - x3)*(xc - x3) + (yc - y3)*(yc - y3) + (zc - z3)*(zc - z3);

        R = R1;
        if (R2 > R) R = R2;
        if (R3 > R) R = R3;
        if (R4 > R) R = R4;

        R = sqrt(R);
    }

    //***************************************
    //***************************************

    inline bool CalculatePosition(Geometry<Node < 3 > >&geom,
                                  const double xc, const double yc, const double zc,
                                  array_1d<double, 3 > & N
                                 )
    {
        double x0 = geom[0].X();
        double y0 = geom[0].Y();
        double x1 = geom[1].X();
        double y1 = geom[1].Y();
        double x2 = geom[2].X();
        double y2 = geom[2].Y();

        double area = CalculateVol(x0, y0, x1, y1, x2, y2);
        double inv_area = 0.0;
        if (area == 0.0)
        {
            KRATOS_THROW_ERROR(std::logic_error, "element with zero area found", "");
        }
        else
        {
            inv_area = 1.0 / area;
        }


        N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
        N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
        N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;


        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
            return true;

        return false;
    }

    //***************************************
    //***************************************

    inline bool CalculatePosition(Geometry<Node < 3 > >&geom,
                                  const double xc, const double yc, const double zc,
                                  array_1d<double, 4 > & N
                                 )
    {

        double x0 = geom[0].X();
        double y0 = geom[0].Y();
        double z0 = geom[0].Z();
        double x1 = geom[1].X();
        double y1 = geom[1].Y();
        double z1 = geom[1].Z();
        double x2 = geom[2].X();
        double y2 = geom[2].Y();
        double z2 = geom[2].Z();
        double x3 = geom[3].X();
        double y3 = geom[3].Y();
        double z3 = geom[3].Z();

        double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

        double inv_vol = 0.0;
        if (vol < 0.0000000000001)
        {
            KRATOS_THROW_ERROR(std::logic_error, "element with zero vol found", "");
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

    inline double CalculateVol(const double x0, const double y0,
                               const double x1, const double y1,
                               const double x2, const double y2
                              )
    {
        return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
    }
    //***************************************
    //***************************************

    inline double CalculateVol(const double x0, const double y0, const double z0,
                               const double x1, const double y1, const double z1,
                               const double x2, const double y2, const double z2,
                               const double x3, const double y3, const double z3
                              )
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

    void ComputeGaussPointPositions(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 4, 3 > & pos, boost::numeric::ublas::bounded_matrix<double, 4, 3 > & N)
    {
        double one_third = 1.0 / 3.0;
        double one_sixt = 1.0 / 6.0;
        double two_third = 2.0 * one_third;

        N(0, 0) = one_sixt;
        N(0, 1) = one_sixt;
        N(0, 2) = two_third;
        N(1, 0) = two_third;
        N(1, 1) = one_sixt;
        N(1, 2) = one_sixt;
        N(2, 0) = one_sixt;
        N(2, 1) = two_third;
        N(2, 2) = one_sixt;
        N(3, 0) = one_third;
        N(3, 1) = one_third;
        N(3, 2) = one_third;


        //first
        pos(0, 0) = one_sixt * geom[0].X() + one_sixt * geom[1].X() + two_third * geom[2].X();
        pos(0, 1) = one_sixt * geom[0].Y() + one_sixt * geom[1].Y() + two_third * geom[2].Y();
        pos(0, 2) = one_sixt * geom[0].Z() + one_sixt * geom[1].Z() + two_third * geom[2].Z();

        //second
        pos(1, 0) = two_third * geom[0].X() + one_sixt * geom[1].X() + one_sixt * geom[2].X();
        pos(1, 1) = two_third * geom[0].Y() + one_sixt * geom[1].Y() + one_sixt * geom[2].Y();
        pos(1, 2) = two_third * geom[0].Z() + one_sixt * geom[1].Z() + one_sixt * geom[2].Z();

        //third
        pos(2, 0) = one_sixt * geom[0].X() + two_third * geom[1].X() + one_sixt * geom[2].X();
        pos(2, 1) = one_sixt * geom[0].Y() + two_third * geom[1].Y() + one_sixt * geom[2].Y();
        pos(2, 2) = one_sixt * geom[0].Z() + two_third * geom[1].Z() + one_sixt * geom[2].Z();

        //fourth
        pos(3, 0) = one_third * geom[0].X() + one_third * geom[1].X() + one_third * geom[2].X();
        pos(3, 1) = one_third * geom[0].Y() + one_third * geom[1].Y() + one_third * geom[2].Y();
        pos(3, 2) = one_third * geom[0].Z() + one_third * geom[1].Z() + one_third * geom[2].Z();

    }

    void ComputeGaussPointPositions(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 16, 3 > & pos, boost::numeric::ublas::bounded_matrix<double, 16, 3 > & N)
    {
        //lower diagonal terms
        double ypos = 1.0 / 12.0;
        int pos_counter = 0;
        for (unsigned int i = 0; i < 4; i++)
        {
            double xpos = 1.0 / 12.0;
            for (unsigned int j = 0; j < 4 - i; j++)
            {
                double N1 = xpos;
                double N2 = ypos;
                double N3 = 1.0 - xpos - ypos;

                pos(pos_counter, 0) = N1 * geom[0].X() + N2 * geom[1].X() + N3 * geom[2].X();
                pos(pos_counter, 1) = N1 * geom[0].Y() + N2 * geom[1].Y() + N3 * geom[2].Y();
                pos(pos_counter, 2) = N1 * geom[0].Z() + N2 * geom[1].Z() + N3 * geom[2].Z();

                N(pos_counter, 0) = N1;
                N(pos_counter, 1) = N2;
                N(pos_counter, 2) = N3;

                xpos += 1.0 / 4.0;
                pos_counter += 1;

            }
            ypos += 1.0 / 4.0;
        }

        //lower diagonal terms
        ypos = 2.0 / 12.0;
        // pos_counter = 8;
        for (unsigned int i = 0; i < 3; i++)
        {
            double xpos = 2.0 / 12.0;
            for (unsigned int j = 0; j < 4 - i; j++)
            {
                double N1 = xpos;
                double N2 = ypos;
                double N3 = 1.0 - xpos - ypos;

                pos(pos_counter, 0) = N1 * geom[0].X() + N2 * geom[1].X() + N3 * geom[2].X();
                pos(pos_counter, 1) = N1 * geom[0].Y() + N2 * geom[1].Y() + N3 * geom[2].Y();
                pos(pos_counter, 2) = N1 * geom[0].Z() + N2 * geom[1].Z() + N3 * geom[2].Z();

                N(pos_counter, 0) = N1;
                N(pos_counter, 1) = N2;
                N(pos_counter, 2) = N3;

                xpos += 1.0 / 4.0;
                pos_counter += 1;

            }
            ypos += 1.0 / 4.0;
        }
    }

    void ConsistentMassMatrix(const double A, boost::numeric::ublas::bounded_matrix<double, 3, 3 > & M)
    {
        double c1 = A / 12.0;
        double c2 = 2.0 * c1;
        M(0, 0) = c2;
        M(0, 1) = c1;
        M(0, 2) = c1;
        M(1, 0) = c1;
        M(1, 1) = c2;
        M(1, 2) = c1;
        M(2, 0) = c1;
        M(2, 1) = c1;
        M(2, 2) = c2;
    }
};

} // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED  defined


