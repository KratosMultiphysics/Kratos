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
#include "lumped_interface_curvature_calculation.h"


namespace Kratos
{

/* Public functions *******************************************************/

/// constructor
LumpedInterfaceCurvatureCalculation::LumpedInterfaceCurvatureCalculation(
        ModelPart& rModelPart)
    : Process(), mrModelPart(rModelPart) {}

/***********************************************************************************/
/***********************************************************************************/

void LumpedInterfaceCurvatureCalculation::Execute(){

    KRATOS_TRY;

    // First node iterator
    const auto it_node_begin = mrModelPart.NodesBegin();

    // First element iterator
    const auto it_element_begin = mrModelPart.ElementsBegin();

    // Current domain size
    const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = it_node_begin + i;

        const double distance = it_node->FastGetSolutionStepValue(DISTANCE);

        if ( distance == 0.0 ){
            it_node->FastGetSolutionStepValue(DISTANCE) = 1.0e-15;
        }
    }

    // Auxiliar containers
    /* Matrix DN_DX, J0;
    Vector N; */

    std::vector<unsigned int> cut_element_id;
    std::vector<Matrix> cut_element_int_N;
    std::vector<GeometryType::ShapeFunctionsGradientsType> cut_element_int_DN_DX;
    std::vector<Vector> cut_element_int_weights;

    ClearVariables();

    // Iterate over the elements, calculate gradient
    #pragma omp parallel for //firstprivate(DN_DX, N, J0)
    for(int i_elem=0; i_elem<static_cast<int>(mrModelPart.Elements().size()); ++i_elem) {
        auto it_elem = it_element_begin + i_elem;
        auto& r_geometry = it_elem->GetGeometry();

        // Current geometry information
        const std::size_t local_space_dimension = r_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_geometry.PointsNumber();
        const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();

        Vector distances( number_of_nodes );
        for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node){
            distances(i_node) = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE);
        }

        unsigned int nneg=0, npos=0;
        for(unsigned int i = 0; i < number_of_nodes; ++i)
        {
            if(distances(i) >= 1.0e-14) npos += 1;
            else if(distances(i) <= -1.0e-14) nneg += 1;
        }
        
        if(nneg != 0 && npos != 0)
        {
            /* for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node){
                r_geometry[i_node].SetValue(IS_NEAR_CUT, 1.0);
            } */

            /* #pragma omp atomic
            cut_element_id.push_back(i_elem); */

            Kratos::unique_ptr<ModifiedShapeFunctions> p_modified_sh_func = nullptr;
            p_modified_sh_func = Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(it_elem->pGetGeometry(), distances);

            Matrix int_N;
            GeometryType::ShapeFunctionsGradientsType int_DN_DX;
            Vector int_weights;

            p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
            int_N,
            int_DN_DX,
            int_weights,
            r_integration_method);

            /* #pragma omp atomic
            cut_element_int_N.push_back(int_N);

            #pragma omp atomic
            cut_element_int_DN_DX.push_back(int_DN_DX);

            #pragma omp atomic
            cut_element_int_weights.push_back(int_weights); */

            const unsigned int number_of_integration_points = int_weights.size();

            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                const Matrix DN_DX = int_DN_DX[point_number];
                const Vector grad_phi = prod(trans(DN_DX), distances);

                for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
                    array_1d<double, 3>& r_gradient = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT/* _AUX */);
                    for(std::size_t k=0; k<dimension; ++k) {
                        #pragma omp atomic
                        r_gradient[k] += int_N(point_number,i_node) * int_weights(point_number) * grad_phi[k];
                    }

                    double& area = r_geometry[i_node].GetValue(AREA_VARIABLE_AUX);

                    #pragma omp atomic
                    area += int_N(point_number,i_node) * int_weights(point_number);
                }
            }

            //KRATOS_INFO("LumpedInterfaceCurvatureCalculation, gradient") << "end a cut element" << std::endl;
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = it_node_begin + i;
        if (it_node->GetValue(AREA_VARIABLE_AUX) > 1.0e-14)
            it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT/* _AUX */) /= it_node->GetValue(AREA_VARIABLE_AUX);
    }

    //KRATOS_INFO("LumpedInterfaceCurvatureCalculation") << "ponderate gradient" << std::endl;

    // #pragma omp parallel for
    // for(int i=0; i < cut_element_id.size(); ++i) {

    //     auto it_elem = it_element_begin + cut_element_id[i];
    //     auto& r_geometry = it_elem->GetGeometry();
    //     KRATOS_INFO("LumpedInterfaceCurvatureCalculation, area variable") << r_geometry[0].GetValue(AREA_VARIABLE_AUX) << std::endl;

    //     const std::size_t number_of_nodes = r_geometry.PointsNumber();
    //     const unsigned int number_of_integration_points = (cut_element_int_weights[i]).size();

    //     std::vector <Vector> values(number_of_nodes);
    //     for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node){
    //         values[i_node] = ZeroVector(dimension);
    //         double norm = 0.0;
    //         const array_1d<double, 3> i_value = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT/* _AUX */);

    //         for(std::size_t i_dim=0; i_dim<dimension; ++i_dim){
    //             (values[i_node])(i_dim) = i_value(i_dim);
    //             norm += i_value(i_dim)*i_value(i_dim);
    //         }
    //         norm = std::sqrt(norm);
            
    //         for(std::size_t i_dim=0; i_dim<dimension; ++i_dim){
    //             (values[i_node])(i_dim) /= norm;
    //         }
    //     }

    //     for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
    //         double divergence = 0.0;
    //         for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
    //             for(std::size_t i_dim=0; i_dim<dimension; ++i_dim) {
    //                 divergence += ((cut_element_int_DN_DX[i])[point_number])(i_node,i_dim)*(values[i_node])(i_dim);
    //             }
    //         }

    //         for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {

    //             double& r_divergence = r_geometry[i_node].FastGetSolutionStepValue(CURVATURE);

    //             #pragma omp atomic
    //             r_divergence += (cut_element_int_N[i])(point_number,i_node) * 
    //                 (cut_element_int_weights[i])(point_number) * divergence;
    //         }
    //     }
    // }

    // Iterate over the elements, calculate divergence
    #pragma omp parallel for //firstprivate(DN_DX, N, J0)
    for(int i_elem=0; i_elem<static_cast<int>(mrModelPart.Elements().size()); ++i_elem) {
        auto it_elem = it_element_begin + i_elem;
        auto& r_geometry = it_elem->GetGeometry();

        // Current geometry information
        const std::size_t local_space_dimension = r_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_geometry.PointsNumber();
        const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();

        Vector distances( number_of_nodes );
        for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node){
            distances(i_node) = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE);
        }

        unsigned int nneg=0, npos=0;
        for(unsigned int i = 0; i < number_of_nodes; ++i)
        {
            if(distances(i) >= 1.0e-14) npos += 1;
            else if(distances(i) <= -1.0e-14) nneg += 1;
        }
        
        if(nneg != 0 && npos != 0)
        {
            Kratos::unique_ptr<ModifiedShapeFunctions> p_modified_sh_func = nullptr;
            p_modified_sh_func = Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(it_elem->pGetGeometry(), distances);

            Matrix int_N;
            GeometryType::ShapeFunctionsGradientsType int_DN_DX;
            Vector int_weights;

            p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
            int_N,
            int_DN_DX,
            int_weights,
            r_integration_method);

            std::vector <Vector> values(number_of_nodes);
            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node){
                values[i_node] = ZeroVector(dimension);
                double norm = 0.0;
                const array_1d<double, 3> i_value = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT/* _AUX */);

                for(std::size_t i_dim=0; i_dim<dimension; ++i_dim){
                    (values[i_node])(i_dim) = i_value(i_dim);
                    norm += i_value(i_dim)*i_value(i_dim);
                }
                norm = std::sqrt(norm);
                
                for(std::size_t i_dim=0; i_dim<dimension; ++i_dim){
                    (values[i_node])(i_dim) /= norm;
                }
            }

            const unsigned int number_of_integration_points = int_weights.size();

            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                const Matrix DN_DX = int_DN_DX[point_number];
                const Vector grad_phi = prod(trans(DN_DX), distances);

                double divergence = 0.0;
                for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
                    for(std::size_t i_dim=0; i_dim<dimension; ++i_dim) {
                        divergence += DN_DX(i_node,i_dim)*(values[i_node])(i_dim);
                    }
                }

                for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {

                    double& r_divergence = r_geometry[i_node].FastGetSolutionStepValue(CURVATURE);

                    #pragma omp atomic
                    r_divergence += int_N(point_number,i_node) * int_weights(point_number) * divergence;

                    //double& area = r_geometry[i_node].GetValue(AREA_VARIABLE_AUX);

                    //#pragma omp atomic
                    //area += int_N(point_number,i_node) * int_weights(point_number);
                }
            }

            //KRATOS_INFO("LumpedInterfaceCurvatureCalculation, divergence") << "end a cut element" << std::endl;
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = it_node_begin + i;
        if (it_node->GetValue(AREA_VARIABLE_AUX) > 1.0e-14)
            it_node->FastGetSolutionStepValue(CURVATURE) /= it_node->GetValue(AREA_VARIABLE_AUX);
    }

    //KRATOS_INFO("LumpedInterfaceCurvatureCalculation") << "ponderate divergence" << std::endl;       

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void LumpedInterfaceCurvatureCalculation::ClearVariables()
{
    const auto it_node_begin = mrModelPart.NodesBegin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=it_node_begin + i;
        it_node->SetValue(AREA_VARIABLE_AUX, 0.0);
        it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT/* _AUX */).clear();
        it_node->FastGetSolutionStepValue(CURVATURE, 0.0);
        /* it_node->SetValue(IS_NEAR_CUT, 0.0); */
    }
}

};  // namespace Kratos.
