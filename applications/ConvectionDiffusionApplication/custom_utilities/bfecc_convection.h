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
#include "processes/compute_nodal_gradient_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/pointer_map_communicator.h"

namespace Kratos
{

template<std::size_t TDim>
class BFECCConvection
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BFECCConvection<TDim>);

    BFECCConvection(
        typename BinBasedFastPointLocator<TDim>::Pointer pSearchStructure,
        const bool PartialDt = false,
        const bool ActivateLimiter = false)
        : mpSearchStructure(pSearchStructure), mActivateLimiter(ActivateLimiter)
    {
    }

    ~BFECCConvection()
    {
    }

    //**********************************************************************************************
    //**********************************************************************************************
    void BFECCconvect(
        ModelPart& rModelPart,
        const Variable< double >& rVar,
        const Variable<array_1d<double,3> >& conv_var,
        const double substeps)
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

        // Allocate non-historical variables
        block_for_each(rModelPart.Nodes(), [&](Node& rNode){
            rNode.SetValue(rVar, 0.0);
        });

        mLimiter.resize(nparticles);
        if (mActivateLimiter){
            CalculateLimiter(rModelPart, rVar);
        } else{
            for (int i = 0; i < nparticles; i++){
                mLimiter[i] = 1.0;
            }
        }

        //FIRST LOOP: estimate rVar(n+1)
        #pragma omp parallel for firstprivate(results,N,N_valid)
        for (int i = 0; i < nparticles; i++)
        {
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

            ModelPart::NodesContainerType::iterator it_particle = rModelPart.NodesBegin() + i;

            Element::Pointer pelement;
            Element::Pointer pelement_valid;

            array_1d<double,3> bckPos = it_particle->Coordinates();
            const array_1d<double,3>& vel = it_particle->FastGetSolutionStepValue(conv_var);
            bool has_valid_elem_pointer = false;
            bool is_found = ConvectBySubstepping(dt,bckPos,vel, N,N_valid, pelement,pelement_valid, result_begin, max_results, -1.0, substeps, conv_var, has_valid_elem_pointer);
            found[i] = is_found;

            if(is_found) {
                //save position backwards
                elem_backward(i) = pelement;
                Ns[i] = N;

                Geometry< Node >& geom = pelement->GetGeometry();
                double phi1 = N[0] * ( geom[0].FastGetSolutionStepValue(rVar,1));
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi1 += N[k] * ( geom[k].FastGetSolutionStepValue(rVar,1) );
                }

                it_particle->FastGetSolutionStepValue(rVar) = phi1;
            }
            else if(has_valid_elem_pointer)
            {
                //save position backwards
                elem_backward(i) = pelement_valid;
                Ns[i] = N_valid;

                Geometry< Node >& geom = pelement_valid->GetGeometry();
                double phi1 = N[0] * ( geom[0].FastGetSolutionStepValue(rVar,1));
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi1 += N_valid[k] * ( geom[k].FastGetSolutionStepValue(rVar,1) );
                }

                it_particle->FastGetSolutionStepValue(rVar) = phi1;
            }
        }

        //now obtain the value AT TIME STEP N by taking it from N+1
        #pragma omp parallel for firstprivate(results,N,N_valid)
        for (int i = 0; i < nparticles; i++)
        {
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

            ModelPart::NodesContainerType::iterator it_particle = rModelPart.NodesBegin() + i;

            Element::Pointer pelement;
            Element::Pointer pelement_valid;

            array_1d<double,3> fwdPos = it_particle->Coordinates();
            const array_1d<double,3>& vel = it_particle->FastGetSolutionStepValue(conv_var,1);
            bool has_valid_elem_pointer = false;
            bool is_found = ConvectBySubstepping(dt,fwdPos,vel, N, N_valid, pelement, pelement_valid, result_begin, max_results, 1.0, substeps, conv_var,has_valid_elem_pointer);

            if(is_found) {
                Geometry< Node >& geom = pelement->GetGeometry();
                double phi_old = N[0] * ( geom[0].FastGetSolutionStepValue(rVar));

                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi_old  += N[k] * ( geom[k].FastGetSolutionStepValue(rVar) );
                }

                //store correction
                const auto limiter_factor = 0.5*mLimiter[i];
                it_particle->GetValue(rVar) = (1.0 + limiter_factor)*it_particle->FastGetSolutionStepValue(rVar,1) - limiter_factor*phi_old;
//                 iparticle->FastGetSolutionStepValue(rVar) = iparticle->GetValue(rVar) - 0.5 * (phi2 - iparticle->FastGetSolutionStepValue(rVar,1));
            }
            else
            {
                it_particle->GetValue(rVar) = it_particle->FastGetSolutionStepValue(rVar,1);
            }
        }

        #pragma omp parallel for
        for (int i = 0; i < nparticles; i++)
        {
            ModelPart::NodesContainerType::iterator it_particle = rModelPart.NodesBegin() + i;
            bool is_found = found[i];
            if(is_found) {
                Vector N = Ns[i];
                Geometry< Node >& geom = elem_backward[i].GetGeometry();
                double phi1 = N[0] * ( geom[0].GetValue(rVar));
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi1 += N[k] * ( geom[k].GetValue(rVar) );
                }

                it_particle->FastGetSolutionStepValue(rVar) = phi1;
            }
//             else
//                 std::cout << "it should find it" << std::endl;
        }

        KRATOS_CATCH("")
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
                    Geometry< Node >& geom = pelement->GetGeometry();

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
                    Geometry< Node >& geom = pelement->GetGeometry();

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

    // ************************************************************************************************************
    // See [Kuzmin et al., Comput. Methods Appl. Mech. Engrg., 322 (2017) 23–41] for more info about this limiter
    // Befor calling make sure that non-historical variable "DISTANCE_GRADIENT" contains the nodal gradient of rVar
    void CalculateLimiter(
        ModelPart& rModelPart,
        const Variable< double >& rVar)
    {
        const double epsilon = 1.0e-15;
        const double power = 2.0;

        const int nparticles = rModelPart.Nodes().size();

        if(static_cast<int>(mSigmaPlus.size()) != nparticles){
            mSigmaPlus.resize(nparticles);
            mSigmaMinus.resize(nparticles);
        }

        auto& r_default_comm = rModelPart.GetCommunicator().GetDataCommunicator();
        GlobalPointersVector< Node > gp_list;

        for (int i_node = 0; i_node < static_cast<int>(rModelPart.NumberOfNodes()); ++i_node){
            auto it_node = rModelPart.NodesBegin() + i_node;
            GlobalPointersVector< Node >& global_pointer_list = it_node->GetValue(NEIGHBOUR_NODES);

            for (unsigned int j = 0; j< global_pointer_list.size(); ++j)
            {
                auto& global_pointer = global_pointer_list(j);
                gp_list.push_back(global_pointer);
            }
        }

        GlobalPointerCommunicator< Node > pointer_comm(r_default_comm, gp_list);

        auto coordinate_proxy = pointer_comm.Apply(
            [](GlobalPointer<Node >& global_pointer) -> Point::CoordinatesArrayType
            {
                return global_pointer->Coordinates();
            }
        );

        auto distance_proxy = pointer_comm.Apply(
            [&](GlobalPointer<Node >& global_pointer) -> double
            {
                return global_pointer->FastGetSolutionStepValue(rVar);
            }
        );

        IndexPartition<int>(nparticles).for_each(
        [&](int i_node){
            auto it_node = rModelPart.NodesBegin() + i_node;
            const auto& X_i = it_node->Coordinates();
            const auto& grad_i = it_node->GetValue(DISTANCE_GRADIENT);

            double S_plus = 0.0;
            double S_minus = 0.0;

            GlobalPointersVector< Node >& global_pointer_list = it_node->GetValue(NEIGHBOUR_NODES);

            for (unsigned int j = 0; j< global_pointer_list.size(); ++j)
            {

                /* if (it_node->Id() == j_node->Id())
                    continue; */

                auto& global_pointer = global_pointer_list(j);
                auto X_j = coordinate_proxy.Get(global_pointer);

                S_plus += std::max(0.0, inner_prod(grad_i, X_i-X_j));
                S_minus += std::min(0.0, inner_prod(grad_i, X_i-X_j));
            }

            mSigmaPlus[i_node] = std::min(1.0, (std::abs(S_minus)+epsilon)/(S_plus+epsilon));
            mSigmaMinus[i_node] = std::min(1.0, (S_plus+epsilon)/(std::abs(S_minus)+epsilon));
        }
        );

        IndexPartition<int>(nparticles).for_each(
        [&](int i_node){
            auto it_node = rModelPart.NodesBegin() + i_node;
            const double distance_i = it_node->FastGetSolutionStepValue(rVar);
            const auto& X_i = it_node->Coordinates();
            const auto& grad_i = it_node->GetValue(DISTANCE_GRADIENT);

            double numerator = 0.0;
            double denominator = 0.0;

            GlobalPointersVector< Node >& global_pointer_list = it_node->GetValue(NEIGHBOUR_NODES);

            for (unsigned int j = 0; j< global_pointer_list.size(); ++j)
            {

                /* if (it_node->Id() == j_node->Id())
                    continue; */

                auto& global_pointer = global_pointer_list(j);
                auto X_j = coordinate_proxy.Get(global_pointer);
                const double distance_j = distance_proxy.Get(global_pointer);

                double beta_ij = 1.0;
                if (inner_prod(grad_i, X_i-X_j) > 0)
                    beta_ij = mSigmaPlus[i_node];
                else if (inner_prod(grad_i, X_i-X_j) < 0)
                    beta_ij = mSigmaMinus[i_node];

                numerator += beta_ij*(distance_i - distance_j);
                denominator += beta_ij*std::abs(distance_i - distance_j);
            }

            const double fraction = (std::abs(numerator)/* +epsilon */) / (denominator + epsilon);
            mLimiter[i_node] = 1.0 - std::pow(fraction, power);
        }
        );
    }


    void ResetBoundaryConditions(ModelPart& rModelPart, const Variable< double >& rVar)
    {
            KRATOS_TRY

            ModelPart::NodesContainerType::iterator inodebegin = rModelPart.NodesBegin();
            vector<int> node_partition;
            #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
            #else
                int number_of_threads = 1;
            #endif
            OpenMPUtils::CreatePartition(number_of_threads, rModelPart.Nodes().size(), node_partition);

            #pragma omp parallel for
            for(int kkk=0; kkk<number_of_threads; kkk++)
            {
                for(int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
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
        vector<int> node_partition;
        #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
        #else
            int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, rModelPart.Nodes().size(), node_partition);

        #pragma omp parallel for
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
                inode->GetSolutionStepValue(rVar,1) = inode->FastGetSolutionStepValue(rVar);
            }
        }
        KRATOS_CATCH("")
    }

protected:
    Kratos::Vector mSigmaPlus, mSigmaMinus, mLimiter;

private:
    typename BinBasedFastPointLocator<TDim>::Pointer mpSearchStructure;
    //const bool mPartialDt;
    const bool mActivateLimiter;



};

} // namespace Kratos.

#endif // KRATOS_BFECC_CONVECTION_INCLUDED  defined
