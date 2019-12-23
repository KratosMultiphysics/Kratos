//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Me
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "processes/compute_nodal_normal_divergence_process.h"
// See the header file

namespace Kratos
{
template<bool THistorical>
void ComputeNodalNormalDivergenceProcess<THistorical>::Execute()
{
    KRATOS_TRY;

    // Set to zero
    ClearDivergence();

    // Auxiliar containers
    Matrix DN_DX, J0;
    Vector N;

    // First element iterator
    const auto it_element_begin = mrModelPart.ElementsBegin();

    // Current domain size
    const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // Iterate over the elements
    #pragma omp parallel for firstprivate(DN_DX,  N, J0)
    for(int i_elem=0; i_elem<static_cast<int>(mrModelPart.Elements().size()); ++i_elem) {
        auto it_elem = it_element_begin + i_elem;
        auto& r_geometry = it_elem->GetGeometry();

        // Current geometry information
        const std::size_t local_space_dimension = r_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_geometry.PointsNumber();

        // Resize if needed
        if (DN_DX.size1() != number_of_nodes || DN_DX.size2() != dimension)
            DN_DX.resize(number_of_nodes, dimension);
        if (N.size() != number_of_nodes)
            N.resize(number_of_nodes);
        if (J0.size1() != dimension || J0.size2() != local_space_dimension)
            J0.resize(dimension, local_space_dimension);

        // The integration points
        const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
        const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
        const std::size_t number_of_integration_points = r_integration_points.size();

        std::vector <Vector> values(number_of_nodes);
        if (!mNonHistoricalVariable) {
            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node){

                values[i_node] = ZeroVector(dimension);
                double norm = 0.0;
                //KRATOS_INFO("Reading OriginVariable") << *mpOriginVariable << std::endl;
                const array_1d<double, 3> i_value = r_geometry[i_node].FastGetSolutionStepValue(*mpOriginVariable);
                //KRATOS_INFO("Read OriginVariable") << i_value << std::endl;

                for(std::size_t i_dim=0; i_dim<dimension; ++i_dim){
                    //KRATOS_INFO("Values") << i_node << ", " << i_dim << std::endl;
                    //KRATOS_INFO("Value =") << i_value(i_dim) << std::endl;
                    (values[i_node])(i_dim) = i_value(i_dim);
                    norm += i_value(i_dim)*i_value(i_dim);
                }
                norm = std::sqrt(norm);
                //KRATOS_INFO("norm") << norm << std::endl;
                for(std::size_t i_dim=0; i_dim<dimension; ++i_dim){
                    (values[i_node])(i_dim) /= norm;
                }
            }
        } else {
            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node){
                values.push_back(ZeroVector(dimension));
                double norm = 0.0;
                const array_1d<double, 3> i_value = r_geometry[i_node].GetValue(*mpOriginVariable);

                for(std::size_t i_dim=0; i_dim<dimension; ++i_dim){
                    (values[i_node])(i_dim) = i_value(i_dim);
                    norm += i_value(i_dim)*i_value(i_dim);
                }
                norm = std::sqrt(norm);
                for(std::size_t i_dim=0; i_dim<dimension; ++i_dim){
                    (values[i_node])(i_dim) /= norm;
                }
            }
        }

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

            double divergence = 0.0;
            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
                for(std::size_t i_dim=0; i_dim<dimension; ++i_dim) {
                    divergence += DN_DX(i_node,i_dim)*(values[i_node])(i_dim);
                }
            }

            const double gauss_point_volume = r_integration_points[point_number].Weight() * detJ0;

            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {

                double& r_divergence = GetDivergence(r_geometry, i_node);

                #pragma omp atomic
                r_divergence += N[i_node] * gauss_point_volume * divergence;

                double& vol = r_geometry[i_node].GetValue(*mpAreaVariable);

                #pragma omp atomic
                vol += N[i_node] * gauss_point_volume;
            }
        }
    }

    PonderateDivergence();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>::ComputeNodalNormalDivergenceProcess(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rOriginVariable,
    const Variable<double>& rDivergenceVariable,
    const Variable<double>& rAreaVariable,
    const bool NonHistoricalVariable
    ) : mrModelPart(rModelPart),
        mpOriginVariable(&rOriginVariable),
        mpDivergenceVariable(&rDivergenceVariable),
        mpAreaVariable(&rAreaVariable),
        mNonHistoricalVariable(NonHistoricalVariable)
{
    KRATOS_TRY

    // Doing several checks
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rDivergenceVariable, mrModelPart.Nodes());

    // In case the area variable is not initialized we initialize it
    auto& r_nodes = rModelPart.Nodes();
    if (!r_nodes.begin()->Has( rAreaVariable )) {
        VariableUtils().SetNonHistoricalVariable(rAreaVariable, 0.0, r_nodes);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>::ComputeNodalNormalDivergenceProcess(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rOriginVariable,
    const Variable<double>& rDivergenceVariable,
    const Variable<double>& rAreaVariable,
    const bool NonHistoricalVariable
    ) : mrModelPart(rModelPart),
        mpOriginVariable(&rOriginVariable),
        mpDivergenceVariable(&rDivergenceVariable),
        mpAreaVariable(&rAreaVariable),
        mNonHistoricalVariable(NonHistoricalVariable)
{
    KRATOS_TRY

    // Doing several checks
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());

    // In case the area or gradient variable is not initialized we initialize it
    auto& r_nodes = rModelPart.Nodes();
    if (!r_nodes.begin()->Has( rDivergenceVariable )) {
        VariableUtils().SetNonHistoricalVariable(rDivergenceVariable, 0.0, r_nodes);
    }
    if (!r_nodes.begin()->Has( rAreaVariable )) {
        VariableUtils().SetNonHistoricalVariable(rAreaVariable, 0.0, r_nodes);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>::ClearDivergence()
{
    const auto it_node_begin = mrModelPart.NodesBegin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=it_node_begin + i;
        it_node->SetValue(*mpAreaVariable, 0.0);
        it_node->FastGetSolutionStepValue(*mpDivergenceVariable) = 0.0;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>::ClearDivergence()
{
    const auto it_node_begin = mrModelPart.NodesBegin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node= it_node_begin + i;
        it_node->SetValue(*mpAreaVariable, 0.0);
        it_node->SetValue(*mpDivergenceVariable, 0.0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
double& ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>::GetDivergence(
    Element::GeometryType& rThisGeometry,
    unsigned int i
    )
{
    return rThisGeometry[i].FastGetSolutionStepValue(*mpDivergenceVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <>
double& ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>::GetDivergence(
    Element::GeometryType& rThisGeometry,
    unsigned int i
    )
{
    return rThisGeometry[i].GetValue(*mpDivergenceVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>::PonderateDivergence()
{
    const auto it_node_begin = mrModelPart.NodesBegin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = it_node_begin + i;
        it_node->FastGetSolutionStepValue(*mpDivergenceVariable) /= it_node->GetValue(*mpAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>::PonderateDivergence()
{
    const auto it_node_begin = mrModelPart.NodesBegin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = it_node_begin + i;
        it_node->GetValue(*mpDivergenceVariable) /= it_node->GetValue(*mpAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
Parameters ComputeNodalNormalDivergenceProcess<THistorical>::GetDefaultParameters() const
{
    Parameters default_parameters = Parameters(R"(
    {
        "origin_variable"                : "PLEASE_DEFINE_A_VARIABLE",
        "non_historical_origin_variable" :  false
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>;
template class ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>;

} /* namespace Kratos.*/
