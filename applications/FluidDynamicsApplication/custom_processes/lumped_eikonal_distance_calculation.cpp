//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Simon Wenczowski, ME
//
//

// System includes

// External includes

// Project includes
// See header file

// Application includes
#include "lumped_eikonal_distance_calculation.h"


namespace Kratos
{

/* Public functions *******************************************************/

/// constructor
LumpedEikonalDistanceCalculation::LumpedEikonalDistanceCalculation(
        ModelPart& rModelPart,
        const int maxNumIterations,
        const double tolerance,
        const double pseudoTimeStep)
    : Process(), mrModelPart(rModelPart) {

    mMaxNumIterations = maxNumIterations;
    mTolerance = tolerance;
    mPseudoTimeStep = pseudoTimeStep;
}


/// Initialization function to find the initial volumes and print first lines in the log-file
void LumpedEikonalDistanceCalculation::Initialize(){

    const auto it_node_begin = mrModelPart.NodesBegin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = it_node_begin + i;
        it_node->SetValue(DISTANCE, it_node->FastGetSolutionStepValue(DISTANCE));
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LumpedEikonalDistanceCalculation::Execute(){

    KRATOS_TRY;

    const auto it_node_begin = mrModelPart.NodesBegin();

    // First element iterator
    const auto it_element_begin = mrModelPart.ElementsBegin();

    // Current domain size
    const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = it_node_begin + i;

        double distance = it_node->FastGetSolutionStepValue(DISTANCE);

        if ( distance == 0.0 ){
            it_node->FastGetSolutionStepValue(DISTANCE) = 1.0e-8;
            distance = 1.0e-8;
        }

        it_node->SetValue(DISTANCE, distance);
    }

    // Auxiliar containers
    Matrix DN_DX, J0;
    Vector N;

    int iter = 0;
    double tol = 1.0;

    while (iter < mMaxNumIterations && tol > mTolerance) {
        tol = mTolerance;
        iter++;

        // Set to zero
        ClearVariables();

        // Iterate over the elements
        #pragma omp parallel for firstprivate(DN_DX, N, J0)
        for(int i_elem=0; i_elem<static_cast<int>(mrModelPart.Elements().size()); ++i_elem) {
            auto it_elem = it_element_begin + i_elem;
            auto& r_geometry = it_elem->GetGeometry();

            // Current geometry information
            const std::size_t local_space_dimension = r_geometry.LocalSpaceDimension();
            const std::size_t number_of_nodes = r_geometry.PointsNumber();
            const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();

            Vector distances( number_of_nodes );
            Vector distances0( number_of_nodes );
            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node){
                distances(i_node) = r_geometry[i_node].GetValue(DISTANCE);
                distances0(i_node) = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE);
            }

            unsigned int nneg=0, npos=0;
            for(unsigned int i = 0; i<4; ++i)
                if(distances0(i) >= 1.0e-8) npos += 1;
                else if(distances0(i) <= -1.0e-8) nneg += 1;
        
            if(nneg==0 || npos==0)
            {
                // Resize if needed
                if (DN_DX.size1() != number_of_nodes || DN_DX.size2() != dimension)
                    DN_DX.resize(number_of_nodes, dimension);
                if (N.size() != number_of_nodes)
                    N.resize(number_of_nodes);
                if (J0.size1() != dimension || J0.size2() != local_space_dimension)
                    J0.resize(dimension, local_space_dimension);

                // The integration points
                const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
                const std::size_t number_of_integration_points = r_integration_points.size();

                // The containers of the shape functions and its local gradient
                const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(r_integration_method);
                const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(r_integration_method);

                for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                    // Getting the shape functions
                    noalias(N) = row(rNcontainer, point_number);

                    // Getting the jacobians and local gradients
                    GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[point_number], J0);
                    double detJ0;
                    Matrix InvJ0;
                    MathUtils<double>::GeneralizedInvertMatrix(J0, InvJ0, detJ0);
                    const Matrix& rDN_De = rDN_DeContainer[point_number];
                    GeometryUtils::ShapeFunctionsGradients(rDN_De, InvJ0, DN_DX);

                    const Vector dist_grad = prod(trans(DN_DX), distances);

                    const double gauss_point_volume = r_integration_points[point_number].Weight() * detJ0;

                    double rhs = gauss_point_volume * (1 - norm_2(dist_grad));
                    if (nneg > 0)
                        rhs = -rhs;

                    for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {

                        double& r_dist_diff = GetDDistance(r_geometry, i_node);

                        #pragma omp atomic
                        r_dist_diff += mPseudoTimeStep * N[i_node] * rhs;

                        double& vol = r_geometry[i_node].GetValue(AREA_VARIABLE_AUX);

                        #pragma omp atomic
                        vol += N[i_node] * gauss_point_volume;
                    }
                }
            }else{
                Matrix neg_shape_functions, pos_shape_functions;
                GeometryType::ShapeFunctionsGradientsType neg_shape_derivatives, pos_shape_derivatives;
                Vector neg_gp_w, pos_gp_w;

                Kratos::unique_ptr<ModifiedShapeFunctions> p_modified_sh_func = nullptr;
                p_modified_sh_func = Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(it_elem->pGetGeometry(), distances0);
                
                p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
                        pos_shape_functions,        // N
                        pos_shape_derivatives,      // DN_DX
                        pos_gp_w,                   // weight * detJ
                        r_integration_method);      // Default

                const std::size_t number_of_pos_integration_points = pos_gp_w.size();

                for ( IndexType point_number = 0; point_number < number_of_pos_integration_points; ++point_number ) {
                    const Matrix& DN_DX_gp = pos_shape_derivatives[point_number];
                    const Vector dist_grad = prod(trans(DN_DX_gp), distances);
                    double rhs = pos_gp_w(point_number) * (1 - norm_2(dist_grad));

                    for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {

                        double& r_dist_diff = GetDDistance(r_geometry, i_node);

                        #pragma omp atomic
                        r_dist_diff += mPseudoTimeStep * pos_shape_functions(point_number, i_node) * rhs;

                        double& vol = r_geometry[i_node].GetValue(AREA_VARIABLE_AUX);

                        #pragma omp atomic
                        vol += pos_shape_functions(point_number, i_node) * pos_gp_w(point_number);
                    }
                }

                p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
                        neg_shape_functions,        // N
                        neg_shape_derivatives,      // DN_DX
                        neg_gp_w,                   // weight * detJ
                        r_integration_method);      // Default

                const std::size_t number_of_neg_integration_points = neg_gp_w.size();

                for ( IndexType point_number = 0; point_number < number_of_neg_integration_points; ++point_number ) {
                    const Matrix& DN_DX_gp = neg_shape_derivatives[point_number];
                    const Vector dist_grad = prod(trans(DN_DX_gp), distances);
                    double rhs = - neg_gp_w(point_number) * (1 - norm_2(dist_grad));

                    for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {

                        double& r_dist_diff = GetDDistance(r_geometry, i_node);

                        #pragma omp atomic
                        r_dist_diff += mPseudoTimeStep * neg_shape_functions(point_number, i_node) * rhs;

                        double& vol = r_geometry[i_node].GetValue(AREA_VARIABLE_AUX);

                        #pragma omp atomic
                        vol += neg_shape_functions(point_number, i_node) * neg_gp_w(point_number);
                    }
                }
            }   
        }

        PonderateDDistance();

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
            auto it_node = it_node_begin + i;

            const double ddistance = it_node->FastGetSolutionStepValue(DISTANCE_AUX);
            if (ddistance > tol)
                tol = ddistance;

            double distance = it_node->GetValue(DISTANCE) + ddistance;

            /* if ( distance == 0.0 ){
                distance = 1.0e-8;
            } */

            it_node->SetValue(DISTANCE, distance);
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = it_node_begin + i;

        const double distance = it_node->GetValue(DISTANCE);
        it_node->FastGetSolutionStepValue(DISTANCE) = distance;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void LumpedEikonalDistanceCalculation::ClearVariables()
{
    const auto it_node_begin = mrModelPart.NodesBegin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=it_node_begin + i;
        it_node->SetValue(AREA_VARIABLE_AUX, 0.0);
        it_node->FastGetSolutionStepValue(DISTANCE_AUX) = 0.0;
    }
}

/***********************************************************************************/
/***********************************************************************************/

double& LumpedEikonalDistanceCalculation::GetDDistance(
    Element::GeometryType& rThisGeometry,
    unsigned int i
    )
{
    return rThisGeometry[i].FastGetSolutionStepValue(DISTANCE_AUX);
}

/***********************************************************************************/
/***********************************************************************************/

void LumpedEikonalDistanceCalculation::PonderateDDistance()
{
    const auto it_node_begin = mrModelPart.NodesBegin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = it_node_begin + i;
        it_node->FastGetSolutionStepValue(DISTANCE_AUX) /= it_node->GetValue(AREA_VARIABLE_AUX);
    }
}

};  // namespace Kratos.