//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    ME
//
//

// System includes

// External includes

// Project includes
// See header file

// Application includes
#include "hyperbolic_distance_reinitialization.h"

namespace Kratos
{

/* Public functions *******************************************************/

// Constructor
template<std::size_t TDim>
HyperbolicDistanceReinitialization<TDim>::HyperbolicDistanceReinitialization(
    ModelPart& rModelPart,
    typename BinBasedFastPointLocator<TDim>::Pointer pSearchStructure,
    ComputeNodalGradientProcess<true>::Pointer pComputeGradient,
    const int maxNumIterations,
    const double tolerance,
    const double pseudoTimeStep)
    : Process(), mrModelPart(rModelPart), mpSearchStructure(pSearchStructure), mpComputeGradient(pComputeGradient)
    {
        mMaxNumIterations = maxNumIterations;
        mTolerance = tolerance;
        mPseudoTimeStep = pseudoTimeStep;
    }

template<std::size_t TDim>
void HyperbolicDistanceReinitialization<TDim>::Execute(){

    KRATOS_TRY;

    const unsigned int num_nodes = mrModelPart.NumberOfNodes();
    double min_h = 1.0/mesh_tolerance;

    #pragma omp parallel for
    for (unsigned int j = 0; j < num_nodes; ++j) {
        auto it_node = mrModelPart.NodesBegin() + j;
        const auto distance1 = it_node->FastGetSolutionStepValue(DISTANCE, 1);
        it_node->SetValue(AUX_DISTANCE, distance1);
        const auto distance = it_node->FastGetSolutionStepValue(DISTANCE);
        it_node->SetValue(DISTANCE, distance);
        const auto distance_gradient = it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT);
        it_node->SetValue(DISTANCE_GRADIENT, distance_gradient);
        const double nodal_h = it_node->GetValue(NODAL_H);
        if (nodal_h < min_h)
        {
            #pragma omp critical
            min_h = nodal_h;
        }

        it_node -> Free(AUX_DISTANCE);  // Do not work with DISTANCE in order not to make any confusion in LS
    }

    const unsigned int num_elements = mrModelPart.NumberOfElements();
    //array_1d<double,TDim+1> distance_array;

    #pragma omp parallel for //private(distance_array)
    for (unsigned int k = 0; k < num_elements; ++k) {
        auto it_elem = mrModelPart.ElementsBegin() + k;
        unsigned int n_neg = 0;
        unsigned int n_pos = 0;
        auto& elem_geometry = it_elem->GetGeometry();

        array_1d<double,TDim+1> distance_array;

        for (unsigned int j = 0; j < elem_geometry.size(); j++) {
            distance_array[j] = elem_geometry[j].GetValue(DISTANCE);//FastGetSolutionStepValue(DISTANCE);
            if (distance_array[j] > 0.0)
                n_pos++;
            else
                n_neg++;
        }

        if (n_pos > 0 && n_neg > 0){
            /* GeometryUtils::CalculateTetrahedraDistances */
            GeometryUtils::CalculateExactDistancesToPlane(elem_geometry, distance_array);

            for (unsigned int j = 0; j < elem_geometry.size(); j++){
                elem_geometry[j].FastGetSolutionStepValue(DISTANCE) = distance_array[j];

                #pragma omp critical
                elem_geometry[j].Fix(AUX_DISTANCE); // Do not work with DISTANCE in order not to make any confusion in LS
            }
        }
    }

    mPseudoTimeStep = std::min( mPseudoTimeStep, min_h/20.0);
    KRATOS_INFO("Hyperbolic Redistance") << "min_h = " << min_h << std::endl;
    KRATOS_INFO("Hyperbolic Redistance") << "mPseudoTimeStep = " << mPseudoTimeStep << std::endl;

    // Updating distance gradient at nodes
    mpComputeGradient->Execute();

    double dist_error = mTolerance + 1.0;
    unsigned int outer_iteration = 0;
    while (dist_error > mTolerance && outer_iteration < mMaxNumIterations){
        // Resetting the error
        dist_error = 0.0;
        outer_iteration++;

        // Updating distance gradient at nodes
        mpComputeGradient->Execute();

        Vector N(TDim + 1);
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        #pragma omp parallel for firstprivate(results,N)
        for (unsigned int i = 0; i < num_nodes; i++){
            auto i_node = mrModelPart.NodesBegin() + i;
            const auto dist = i_node->FastGetSolutionStepValue(DISTANCE);
            i_node->FastGetSolutionStepValue(DISTANCE, 1) = dist;
        }

        #pragma omp parallel for firstprivate(results,N)
        for (unsigned int i = 0; i < num_nodes; i++)
        {
            auto i_node = mrModelPart.NodesBegin() + i;

            if (i_node->IsFixed(AUX_DISTANCE))
                continue; // Do not work with DISTANCE in order not to make any confusion in LS

            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
            Element::Pointer pelement;

            const array_1d<double,3> node_position = i_node->Coordinates();
            const auto dist0 = i_node->GetValue(DISTANCE);
            const auto dist = i_node->FastGetSolutionStepValue(DISTANCE);
            const auto dist_grad = i_node->FastGetSolutionStepValue(DISTANCE_GRADIENT);
            const auto dist_grad_norm = norm_2(dist_grad);
            if (dist_grad_norm < mesh_tolerance)
                KRATOS_INFO("Hyperbolic Redistance") << "dist_grad_norm = " << dist_grad_norm << std::endl;

            double sgn0 = 0.0;
            if (dist0 > mesh_tolerance){sgn0 = 1.0;}
            else if (dist0 < -mesh_tolerance){sgn0 = -1.0;}

            const array_1d<double,3> vel = sgn0/dist_grad_norm*dist_grad;
            //KRATOS_INFO("Hyperbolic Redistance") << "vel = " << vel << std::endl;

            bool is_found = false;

            const auto position = node_position - mPseudoTimeStep*vel;
            is_found = mpSearchStructure->FindPointOnMesh(position, N, pelement, result_begin, max_results, mesh_tolerance);

            if (is_found) {
                Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                double phi1 = N[0] * ( geom[0].FastGetSolutionStepValue(DISTANCE,1) );
                for (unsigned int k = 1; k < geom.size(); k++) {
                    phi1 += N[k] * ( geom[k].FastGetSolutionStepValue(DISTANCE,1) );
                }

                i_node->FastGetSolutionStepValue(DISTANCE) = phi1 + sgn0*mPseudoTimeStep;
                //KRATOS_INFO("Hyperbolic Redistance") << "Found a regular node" << std::endl;
            }
            else if (i_node->Is(SLIP) && sgn0 != 0.0)
            {
                auto normal = i_node->FastGetSolutionStepValue(NORMAL);
                normal /= norm_2(normal);
                const array_1d<double,3> vel_blind = /* sgn0* */normal/dist_grad_norm;
                //KRATOS_INFO("Hyperbolic Redistance") << "vel_blind = " << vel_blind << std::endl;

                const auto position_blind = node_position - mPseudoTimeStep*vel_blind;
                bool is_found_blind = false;
                is_found_blind = mpSearchStructure->FindPointOnMesh(position_blind, N, pelement, result_begin, max_results, mesh_tolerance);

                if (is_found_blind) {
                    Geometry< Node < 3 > >& geom = pelement->GetGeometry();
                    double phi1 = N[0] * ( geom[0].FastGetSolutionStepValue(DISTANCE,1) );
                    for (unsigned int k = 1; k < geom.size(); k++) {
                        phi1 += N[k] * ( geom[k].FastGetSolutionStepValue(DISTANCE,1) );
                    }

                    const auto dist_grad0 = i_node->GetValue(DISTANCE_GRADIENT);
                    const auto normal_dist_grad0 = inner_prod(dist_grad0, normal);

                    i_node->FastGetSolutionStepValue(DISTANCE) = phi1 +
                        /* sgn0* */normal_dist_grad0/norm_2(dist_grad0)*mPseudoTimeStep;
                } else {
                    KRATOS_INFO("Hyperbolic Redistance") << "not found at all!, SLIP" << position_blind << std::endl;
                    KRATOS_INFO("Hyperbolic Redistance") << "vel_blind = " << vel_blind << std::endl;
                    KRATOS_INFO("Hyperbolic Redistance") << "normal = " << normal << std::endl;
                    KRATOS_INFO("Hyperbolic Redistance") << "dist_grad_norm = " << dist_grad_norm << std::endl;
                }
            }
            else
            {
                if (sgn0 != 0.0)
                    KRATOS_INFO("Hyperbolic Redistance") << "not found at all:" << position << std::endl;
                    KRATOS_INFO("Hyperbolic Redistance") << "vel = " << vel << std::endl;
                    KRATOS_INFO("Hyperbolic Redistance") << "dist_grad_norm = " << dist_grad_norm << std::endl;
                    KRATOS_INFO("Hyperbolic Redistance") << "dist_grad = " << dist_grad << std::endl;
            }

            const auto d_dist = std::abs( dist - i_node->FastGetSolutionStepValue(DISTANCE) );
            if (d_dist > dist_error)
            {
                #pragma omp critical
                dist_error = d_dist;
            }
        }
    }

    #pragma omp parallel for
    for (unsigned int k = 0; k < num_nodes; ++k) {
        auto it_node = mrModelPart.NodesBegin() + k;
        it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT)=
            it_node->GetValue(DISTANCE_GRADIENT);
        it_node->FastGetSolutionStepValue(DISTANCE, 1) =
            it_node->GetValue(AUX_DISTANCE);
        it_node -> Free(AUX_DISTANCE);
    }

    KRATOS_CATCH("")
}

template class Kratos::HyperbolicDistanceReinitialization<2>;
template class Kratos::HyperbolicDistanceReinitialization<3>;

};  // namespace Kratos.