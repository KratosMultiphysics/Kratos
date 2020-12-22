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

#if !defined(KRATOS_BFECC_CONVECTION_INCLUDED )
#define  KRATOS_BFECC_CONVECTION_INCLUDED

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
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/variables.h"
#include "utilities/timer.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

template<std::size_t TDim> class BFECCConvection
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BFECCConvection<TDim>);

    BFECCConvection(typename BinBasedFastPointLocator<TDim>::Pointer pSearchStructure)
        : mpSearchStructure(pSearchStructure)
    {
    }

    ~BFECCConvection()
    {
    }

    //**********************************************************************************************
    //**********************************************************************************************
    void BFECCconvectPartially(ModelPart& rModelPart, const Variable< double >& rVar, const Variable<array_1d<double,3> >& conv_var, const double substeps, const double dt_factor)
    {
        KRATOS_TRY
        const double dt = dt_factor*rModelPart.GetProcessInfo()[DELTA_TIME];

        //do movement
        Vector N(TDim + 1);
        Vector N_valid(TDim + 1);
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        const int nparticles = rModelPart.Nodes().size();

        mSigmaPlus.resize(nparticles);
        mSigmaMinus.resize(nparticles);

        PointerVector< Element > elem_backward( rModelPart.Nodes().size());
        std::vector< Vector > Ns( rModelPart.Nodes().size());
        std::vector< bool > found( rModelPart.Nodes().size());

        // Allocate non-historical variables and update old velocity as per dt_factor
        #pragma omp parallel for
        for (int i = 0; i < nparticles; i++) {
            auto iparticle = rModelPart.NodesBegin() + i;
            iparticle->SetValue(rVar, 0.0);
            const auto old_velocity = iparticle->FastGetSolutionStepValue(conv_var, 1);
            const auto current_velocity = iparticle->FastGetSolutionStepValue(conv_var);
            iparticle->SetValue(conv_var, old_velocity);
            iparticle->FastGetSolutionStepValue(conv_var, 1) = dt_factor*old_velocity + (1.0 - dt_factor)*current_velocity;
        }

        // ****************************************************************************************
        // ****************************************************************************************
        // Calculating nodal limiter using \beta_ij = 1 (works fine on symmetric structural meshes)
        // D. Kuzmin et al. / Comput. Methods Appl. Mech. Engrg. 322 (2017) 23–41
        const double epsilon = 1.0e-15;
        const double power = 4;

        #pragma omp parallel for
        for (unsigned int i_node = 0; i_node < static_cast<int>(rModelPart.NumberOfNodes()); ++i_node){
            auto it_node = rModelPart.NodesBegin() + i_node;
            const auto X_i = it_node->Coordinates();
            const auto grad_i = it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT);

            double S_plus = 0.0;
            double S_minus = 0.0;

            for( GlobalPointersVector< Node<3> >::iterator j_node = it_node->GetValue(NEIGHBOUR_NODES).begin();
                j_node != it_node->GetValue(NEIGHBOUR_NODES).end(); ++j_node){

                if (it_node->Id() == j_node->Id())
                    continue;

                const auto X_j = j_node->Coordinates();

                S_plus += std::max(0.0, inner_prod(grad_i, X_i-X_j));
                S_minus += std::min(0.0, inner_prod(grad_i, X_i-X_j));
            }

            mSigmaPlus[i_node] = std::min(1.0, (std::abs(S_minus)+epsilon)/(S_plus+epsilon));
            mSigmaMinus[i_node] = std::min(1.0, (S_plus+epsilon)/(std::abs(S_minus)+epsilon));
        }

        #pragma omp parallel for
        for (unsigned int i_node = 0; i_node < static_cast<int>(rModelPart.NumberOfNodes()); ++i_node){
            auto it_node = rModelPart.NodesBegin() + i_node;
            const double distance_i = it_node->FastGetSolutionStepValue(rVar);
            const auto X_i = it_node->Coordinates();
            const auto grad_i = it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT);

            double numerator = 0.0;
            double denominator = 0.0;

            for( GlobalPointersVector< Node<3> >::iterator j_node = it_node->GetValue(NEIGHBOUR_NODES).begin();
                j_node != it_node->GetValue(NEIGHBOUR_NODES).end(); ++j_node){

                if (it_node->Id() == j_node->Id())
                    continue;

                const double distance_j = j_node->FastGetSolutionStepValue(rVar);
                const auto X_j = j_node->Coordinates();

                double beta_ij = 1.0;
                if (inner_prod(grad_i, X_i-X_j) > 0)
                    beta_ij = mSigmaPlus[i_node];
                else if (inner_prod(grad_i, X_i-X_j) < 0)
                    beta_ij = mSigmaMinus[i_node];

                numerator += beta_ij*(distance_i - distance_j);
                denominator += beta_ij*std::abs(distance_i - distance_j);
            }

            const double fraction = (std::abs(numerator)+epsilon) / (denominator + epsilon);
            const double limiter_i = 1.0 - std::pow(fraction, power);
            it_node->SetValue(LIMITER_COEFFICIENT, limiter_i);
        }

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
            bool has_valid_elem_pointer = false;
            bool is_found = ConvectBySubstepping(dt,bckPos,vel, N,N_valid, pelement,pelement_valid, result_begin, max_results, -1.0, substeps, conv_var, has_valid_elem_pointer);
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
            bool is_found = ConvectBySubstepping(dt,fwdPos,vel, N, N_valid, pelement, pelement_valid, result_begin, max_results, 1.0, substeps, conv_var,has_valid_elem_pointer);

            if(is_found) {
                Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                double phi_old = N[0] * ( geom[0].FastGetSolutionStepValue(rVar));

                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi_old  += N[k] * ( geom[k].FastGetSolutionStepValue(rVar) );
                }

                //store correction
                const double compensated_error = iparticle->GetValue(LIMITER_COEFFICIENT)*0.5*(
                    phi_old - iparticle->FastGetSolutionStepValue(rVar,1));

                iparticle->GetValue(rVar) = iparticle->FastGetSolutionStepValue(rVar,1) - compensated_error; //1.5*iparticle->FastGetSolutionStepValue(rVar,1) - 0.5*phi_old;
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

            iparticle->FastGetSolutionStepValue(conv_var, 1) = iparticle->GetValue(conv_var); // Restoring the old velocity
        }

        KRATOS_CATCH("")
    }

    void BFECCconvect(ModelPart& rModelPart, const Variable< double >& rVar, const Variable<array_1d<double,3> >& conv_var, const double substeps)
    {
        BFECCconvectPartially(rModelPart, rVar, conv_var, substeps, 1.0);
    }

    bool ConvectBySubstepping(
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
        const double subdivisions,
        const Variable<array_1d<double,3> >& conv_var,
        bool& has_valid_elem_pointer)
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

                    noalias(veulerian) = N[0] * ( new_step_factor*geom[0].FastGetSolutionStepValue(conv_var) + old_step_factor*geom[0].FastGetSolutionStepValue(conv_var,1));
                    for (unsigned int k = 1; k < geom.size(); k++)
                        noalias(veulerian) += N[k] * ( new_step_factor*geom[k].FastGetSolutionStepValue(conv_var) + old_step_factor*geom[k].FastGetSolutionStepValue(conv_var,1) );

                    noalias(position) += small_dt*veulerian;

                    N_valid  = N;
                    pelement_valid = pelement;
                    has_valid_elem_pointer = true;

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

                    noalias(veulerian) = N[0] * ( new_step_factor*geom[0].FastGetSolutionStepValue(conv_var) + old_step_factor*geom[0].FastGetSolutionStepValue(conv_var,1));
                    for (unsigned int k = 1; k < geom.size(); k++)
                        noalias(veulerian) += N[k] * ( new_step_factor*geom[k].FastGetSolutionStepValue(conv_var) + old_step_factor*geom[k].FastGetSolutionStepValue(conv_var,1) );

                    noalias(position) -= small_dt*veulerian;

                    N_valid  = N;
                    pelement_valid = pelement;
                    has_valid_elem_pointer = true;


                }
             else
                 break;
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

protected:
    std::vector< double > mSigmaPlus, mSigmaMinus;

private:
    typename BinBasedFastPointLocator<TDim>::Pointer mpSearchStructure;



};

} // namespace Kratos.

#endif // KRATOS_BFECC_CONVECTION_INCLUDED  defined
