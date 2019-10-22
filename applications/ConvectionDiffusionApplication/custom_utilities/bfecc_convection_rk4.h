// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

#if !defined(KRATOS_BFECC_CONVECTION_RK4_INCLUDED )
#define  KRATOS_BFECC_CONVECTION_RK4_INCLUDED

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
#include "utilities/binbased_fast_point_locator.h"


#include <boost/timer.hpp>
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"


namespace Kratos
{

template<std::size_t TDim> class BFECCConvectionRK4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BFECCConvectionRK4<TDim>);

    BFECCConvectionRK4(typename BinBasedFastPointLocator<TDim>::Pointer pSearchStructure)
        : mpSearchStructure(pSearchStructure)
    {
    }

    ~BFECCConvectionRK4()
    {
    }

    //**********************************************************************************************
    //**********************************************************************************************
    void BFECCconvectRK4(ModelPart& rModelPart, const Variable< double >& rVar, const Variable<array_1d<double,3> >& conv_var)
    {
        KRATOS_TRY
        const double dt = rModelPart.GetProcessInfo()[DELTA_TIME];

        //do movement
        Vector N(TDim + 1);
        Vector N_valid(TDim + 1);
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        const int nparticles = rModelPart.Nodes().size();

        PointerVector< Element > elem_backward( rModelPart.Nodes().size());
        std::vector< Vector > Ns( rModelPart.Nodes().size());
        std::vector< bool > found( rModelPart.Nodes().size());
        int counter=1;

        //FIRST LOOP: estimate rVar(n+1)
        #pragma omp parallel for firstprivate(results,N,N_valid)
        for (int i = 0; i < nparticles; i++)
        {
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

            ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;

            Element::Pointer pelement;
            Element::Pointer pelement_valid;

            array_1d<double,3> bckPos = iparticle->Coordinates();
            const array_1d<double,3>& vel = iparticle->FastGetSolutionStepValue(conv_var);
//             if (counter==4)
//             {
//             KRATOS_WATCH(*iparticle);
//             KRATOS_WATCH(vel);
//             KRATOS_WATCH(rVar);
//             }
            bool has_valid_elem_pointer = false;
            bool is_found = ConvectByRK4(dt,bckPos,vel, N,N_valid, pelement,pelement_valid, result_begin, max_results, -1.0, conv_var, has_valid_elem_pointer,counter);

            found[i] = is_found;
            
//            if (counter==4)
//             KRATOS_WATCH(is_found);

            if(is_found) {
                //save position backwards
                elem_backward(i) = pelement;
                Ns[i] = N;

                Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                double phi1 = N[0] * ( geom[0].FastGetSolutionStepValue(rVar,1));
//                 if (counter==4)
//                 {
//                  KRATOS_WATCH(phi1);
//                 }
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi1 += N[k] * ( geom[k].FastGetSolutionStepValue(rVar,1) );
                }

                iparticle->FastGetSolutionStepValue(rVar) = phi1;
            }
            else if(has_valid_elem_pointer)
            {
                //save position backwards
                elem_backward(i) = pelement_valid;
                Ns[i] = N_valid;

                Geometry< Node < 3 > >& geom = pelement_valid->GetGeometry();
                double phi1 = N[0] * ( geom[0].FastGetSolutionStepValue(rVar,1));
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi1 += N_valid[k] * ( geom[k].FastGetSolutionStepValue(rVar,1) );
                }

                iparticle->FastGetSolutionStepValue(rVar) = phi1;
            }
                    
                    
         counter++;
        }
        


        //now obtain the value AT TIME STEP N by taking it from N+1
        #pragma omp parallel for firstprivate(results,N,N_valid)
        for (int i = 0; i < nparticles; i++)
        {
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

            ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;

            Element::Pointer pelement;
            Element::Pointer pelement_valid;

            array_1d<double,3> fwdPos = iparticle->Coordinates();
            const array_1d<double,3>& vel = iparticle->FastGetSolutionStepValue(conv_var,1);
            bool has_valid_elem_pointer = false;
            bool is_found = ConvectByRK4(dt,fwdPos,vel, N, N_valid, pelement, pelement_valid, result_begin, max_results, 1.0, conv_var,has_valid_elem_pointer,counter);

            if(is_found) {
                Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                double phi_old = N[0] * ( geom[0].FastGetSolutionStepValue(rVar));

                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi_old  += N[k] * ( geom[k].FastGetSolutionStepValue(rVar) );
                }

                //store correction
                iparticle->GetValue(rVar) = 1.5*iparticle->FastGetSolutionStepValue(rVar,1) - 0.5 * phi_old;
//                 iparticle->FastGetSolutionStepValue(rVar) = iparticle->GetValue(rVar) - 0.5 * (phi2 - iparticle->FastGetSolutionStepValue(rVar,1));
            }
            else
            {
                iparticle->GetValue(rVar) = iparticle->FastGetSolutionStepValue(rVar,1);
            }
        }

         #pragma omp parallel for
        for (int i = 0; i < nparticles; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = rModelPart.NodesBegin() + i;
            bool is_found = found[i];
            if(is_found) {
                Vector N = Ns[i];
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

        KRATOS_CATCH("")
    }

    bool ConvectByRK4(
        const double dt,
        array_1d<double,3>& position, //IT WILL BE MODIFIED
        const array_1d<double,3>& initial_velocity,
        Vector& N,
        Vector& N_valid,
        Element::Pointer& pelement,
        Element::Pointer& pelement_valid,
        typename BinBasedFastPointLocator<TDim>::ResultIteratorType& result_begin,
        const unsigned int max_results,
        const double velocity_sign,
        const Variable<array_1d<double,3> >& conv_var,
        bool& has_valid_elem_pointer,
        int& counter
                             )
    {
        bool is_found = false;
        bool is_found2 = false;
        bool is_found3 = false;
        bool is_found4 = false;
        array_1d<double,3> veulerian;
        Element::Pointer pelement_aux;
        array_1d<double,3> position_aux;
        array_1d<double,3> k1=ZeroVector(3);
        array_1d<double,3> k2=ZeroVector(3);
        array_1d<double,3> k3=ZeroVector(3);
        array_1d<double,3> k4=ZeroVector(3);
        array_1d<double,3> k2_aux=ZeroVector(3);
        array_1d<double,3> k3_aux=ZeroVector(3);
        array_1d<double,3> k4_aux=ZeroVector(3);
        
//                    if (counter==4)
//             KRATOS_WATCH("insideconvectbysubstepping");

        if(velocity_sign > 0.0) //going from the past to the future
        {
            //RK4 first step
            k1 = initial_velocity*dt;
            k2_aux=position+k1/2.0;
            //RK4 second step
            is_found2 = mpSearchStructure->FindPointOnMesh(k2_aux, N, pelement_aux, result_begin, max_results);
            if (is_found2 == true)
            {
             Geometry< Node < 3 > >& geom = pelement_aux->GetGeometry(); 
             noalias(veulerian) = N[0] * (0.5*geom[0].FastGetSolutionStepValue(conv_var)+0.5*geom[0].FastGetSolutionStepValue(conv_var,1));
             for (unsigned int k = 1; k < geom.size(); k++)
              noalias(veulerian) += N[k]*(0.5*geom[k].FastGetSolutionStepValue(conv_var)+0.5*geom[k].FastGetSolutionStepValue(conv_var,1)); 
             k2 = veulerian*dt;
             k3_aux=position+k2/2.0;
            }
            //RK4 third step
            is_found3 = mpSearchStructure->FindPointOnMesh(k3_aux, N, pelement_aux, result_begin, max_results);
            if (is_found3 == true)
            {
             Geometry< Node < 3 > >& geom = pelement_aux->GetGeometry(); 
             noalias(veulerian) = N[0] * (0.5*geom[0].FastGetSolutionStepValue(conv_var)+0.5*geom[0].FastGetSolutionStepValue(conv_var,1));
             for (unsigned int k = 1; k < geom.size(); k++)
              noalias(veulerian) += N[k]*(0.5*geom[k].FastGetSolutionStepValue(conv_var)+0.5*geom[k].FastGetSolutionStepValue(conv_var,1)); 
             k3 = veulerian*dt;
             k4_aux=position+k3;
            }
            //RK4 fourth step
            is_found4 = mpSearchStructure->FindPointOnMesh(k4_aux, N, pelement_aux, result_begin, max_results);
            if (is_found4 == true)
            {
             Geometry< Node < 3 > >& geom = pelement_aux->GetGeometry(); 
             noalias(veulerian) = N[0] * (geom[0].FastGetSolutionStepValue(conv_var));
             for (unsigned int k = 1; k < geom.size(); k++)
              noalias(veulerian) += N[k]*(geom[k].FastGetSolutionStepValue(conv_var)); 
             k4 = veulerian*dt;
            }
            if (is_found2==1 & is_found3==1 & is_found4==1)
            {
              position_aux=position+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
              is_found = mpSearchStructure->FindPointOnMesh(position_aux, N, pelement, result_begin, max_results);
              if (is_found == true)
              {
               noalias(position) += (1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
               N_valid  = N;
               pelement_valid = pelement;
               has_valid_elem_pointer = true; 
              }
            }
        }
        else //going from the future to the past
        {
            //RK4 first step
            k1 = initial_velocity*dt;
            k2_aux=position-k1/2.0;
            //RK4 second step
            is_found2 = mpSearchStructure->FindPointOnMesh(k2_aux, N, pelement_aux, result_begin, max_results);
            if (is_found2 == true)
            {
             Geometry< Node < 3 > >& geom = pelement_aux->GetGeometry(); 
             noalias(veulerian) = N[0] * (0.5*geom[0].FastGetSolutionStepValue(conv_var)+0.5*geom[0].FastGetSolutionStepValue(conv_var,1));
             for (unsigned int k = 1; k < geom.size(); k++)
              noalias(veulerian) += N[k]*(0.5*geom[k].FastGetSolutionStepValue(conv_var)+0.5*geom[k].FastGetSolutionStepValue(conv_var,1)); 
             k2 = veulerian*dt;
             k3_aux=position-k2/2.0;
            }
            //RK4 third step
            is_found3 = mpSearchStructure->FindPointOnMesh(k3_aux, N, pelement_aux, result_begin, max_results);
            if (is_found3 == true)
            {
             Geometry< Node < 3 > >& geom = pelement_aux->GetGeometry(); 
             noalias(veulerian) = N[0] * (0.5*geom[0].FastGetSolutionStepValue(conv_var)+0.5*geom[0].FastGetSolutionStepValue(conv_var,1));
             for (unsigned int k = 1; k < geom.size(); k++)
              noalias(veulerian) += N[k]*(0.5*geom[k].FastGetSolutionStepValue(conv_var)+0.5*geom[k].FastGetSolutionStepValue(conv_var,1)); 
             k3 = veulerian*dt;
             k4_aux=position-k3;
            }
            //RK4 fourth step
            is_found4 = mpSearchStructure->FindPointOnMesh(k4_aux, N, pelement_aux, result_begin, max_results);
            if (is_found4 == true)
            {
             Geometry< Node < 3 > >& geom = pelement_aux->GetGeometry(); 
             noalias(veulerian) = N[0] * (geom[0].FastGetSolutionStepValue(conv_var,1));
             for (unsigned int k = 1; k < geom.size(); k++)
              noalias(veulerian) += N[k]*(geom[k].FastGetSolutionStepValue(conv_var,1)); 
             k4 = veulerian*dt;
            }
            if (is_found2==1 & is_found3==1 & is_found4==1)
            {
              position_aux=position-(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
              is_found = mpSearchStructure->FindPointOnMesh(position_aux, N, pelement, result_begin, max_results);
              if (is_found == true)
              {
               noalias(position) -= (1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
               N_valid  = N;
               pelement_valid = pelement;
               has_valid_elem_pointer = true; 
              }
            }
        }

                return is_found;

    }


        void ResetBoundaryConditions(ModelPart& rModelPart, const Variable< double >& rVar)
        {
                KRATOS_TRY

                ModelPart::NodesContainerType::iterator inodebegin = rModelPart.NodesBegin();
                vector<unsigned int> node_partition;
                #ifdef _OPENMP
                    int number_of_threads = omp_get_max_threads();
                #else
                    int number_of_threads = 1;
                #endif
                OpenMPUtils::CreatePartition(number_of_threads, rModelPart.Nodes().size(), node_partition);

                #pragma omp parallel for
                for(int kkk=0; kkk<number_of_threads; kkk++)
                {
                    for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
                    {
                            ModelPart::NodesContainerType::iterator inode = inodebegin+ii;

                            if (inode->IsFixed(rVar))
                            {
                                inode->FastGetSolutionStepValue(rVar)=inode->GetSolutionStepValue(rVar,1);
                            }
                    }
                }

                KRATOS_CATCH("")
        }

        void CopyScalarVarToPreviousTimeStep(ModelPart& rModelPart, const Variable< double >& rVar)
        {
            KRATOS_TRY
            ModelPart::NodesContainerType::iterator inodebegin = rModelPart.NodesBegin();
            vector<unsigned int> node_partition;
            #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
            #else
                int number_of_threads = 1;
            #endif
            OpenMPUtils::CreatePartition(number_of_threads, rModelPart.Nodes().size(), node_partition);

            #pragma omp parallel for
            for(int kkk=0; kkk<number_of_threads; kkk++)
            {
                for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
                {
                    ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
                    inode->GetSolutionStepValue(rVar,1) = inode->FastGetSolutionStepValue(rVar);
                }
            }
            KRATOS_CATCH("")
        }
private:
    typename BinBasedFastPointLocator<TDim>::Pointer mpSearchStructure;



};

} // namespace Kratos.

#endif // KRATOS_BFECC_CONVECTION_RK4_INCLUDED  defined


