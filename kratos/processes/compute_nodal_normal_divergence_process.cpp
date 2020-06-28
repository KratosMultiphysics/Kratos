//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi (based on the work by Riccardo Rossi and Vicente Mataix Ferrandiz)
//

/* Project includes */
#include "processes/compute_nodal_normal_divergence_process.h"

namespace Kratos
{
template<bool THistorical>
void ComputeNodalNormalDivergenceProcess<THistorical>::Execute()
{
    KRATOS_TRY;

    // Set to zero
    ClearDivergence();

    // Auxiliary containers
    Matrix J0, InvJ0, DN_DX;

    // First element iterator
    const auto it_element_begin = mrModelPart.ElementsBegin();

    if (!mNonHistoricalOriginVariable) {
        // Iterate over the elements
        #pragma omp parallel for firstprivate(J0, InvJ0, DN_DX)
        for(int i_elem=0; i_elem<static_cast<int>(mrModelPart.Elements().size()); ++i_elem) {
            auto it_elem = it_element_begin + i_elem;
            auto& r_geometry = it_elem->GetGeometry();

            // Current geometry information
            const std::size_t number_of_nodes = r_geometry.PointsNumber();

            // The integration points
            const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
            const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
            const std::size_t number_of_integration_points = r_integration_points.size();

            // The containers of the shape functions and their local gradient
            const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(r_integration_method);
            const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(r_integration_method);

            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Getting the shape functions
                const auto& N = row(rNcontainer, point_number);

                // Getting the jacobians and local shape functions gradient
                double detJ0;
                GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[point_number], J0);
                MathUtils<double>::GeneralizedInvertMatrix(J0, InvJ0, detJ0);
                const Matrix& rDN_De = rDN_DeContainer[point_number];
                GeometryUtils::ShapeFunctionsGradients(rDN_De, InvJ0, DN_DX);

                double divergence = 0.0;
                for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
                    const auto& vector_field = r_geometry[i_node].FastGetSolutionStepValue(*mpOriginVariable);

                    const double norm = norm_2(vector_field);

                    divergence += inner_prod( row(DN_DX, i_node), vector_field ) / norm;
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
    } else{
        // Iterate over the elements
        #pragma omp parallel for firstprivate(J0, InvJ0, DN_DX)
        for(int i_elem=0; i_elem<static_cast<int>(mrModelPart.Elements().size()); ++i_elem) {
            auto it_elem = it_element_begin + i_elem;
            auto& r_geometry = it_elem->GetGeometry();

            // Current geometry information
            const std::size_t number_of_nodes = r_geometry.PointsNumber();

            // The integration points
            const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
            const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
            const std::size_t number_of_integration_points = r_integration_points.size();

            // The containers of the shape functions and their local gradient
            const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(r_integration_method);
            const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(r_integration_method);

            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Getting the shape functions
                const auto& N = row(rNcontainer, point_number);

                // Getting the jacobians and local shape functions gradient
                double detJ0;
                GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[point_number], J0);
                MathUtils<double>::GeneralizedInvertMatrix(J0, InvJ0, detJ0);
                const Matrix& rDN_De = rDN_DeContainer[point_number];
                GeometryUtils::ShapeFunctionsGradients(rDN_De, InvJ0, DN_DX);

                double divergence = 0.0;
                for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
                    const auto& vector_field = r_geometry[i_node].GetValue(*mpOriginVariable);

                    const double norm = norm_2(vector_field);

                    divergence += inner_prod( row(DN_DX, i_node), vector_field ) / norm;
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
    const bool NonHistoricalOriginVariable
    ) : mrModelPart(rModelPart),
        mpOriginVariable(&rOriginVariable),
        mpDivergenceVariable(&rDivergenceVariable),
        mpAreaVariable(&rAreaVariable),
        mNonHistoricalOriginVariable(NonHistoricalOriginVariable)
{
    KRATOS_TRY

    // Doing several checks
    if (!mNonHistoricalOriginVariable) {
        VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    } else{
        KRATOS_ERROR_IF_NOT(mrModelPart.Nodes().begin()->Has(rOriginVariable)) << "Variable " << rOriginVariable.Name() << " not defined on non-historial database" << std::endl;
    }
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
    const bool NonHistoricalOriginVariable
    ) : mrModelPart(rModelPart),
        mpOriginVariable(&rOriginVariable),
        mpDivergenceVariable(&rDivergenceVariable),
        mpAreaVariable(&rAreaVariable),
        mNonHistoricalOriginVariable(NonHistoricalOriginVariable)
{
    KRATOS_TRY

    // Doing several checks
    if (!mNonHistoricalOriginVariable) {
        VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    } else{
        KRATOS_ERROR_IF_NOT(mrModelPart.Nodes().begin()->Has(rOriginVariable)) << "Variable " << rOriginVariable.Name() << " not defined on non-historial database" << std::endl;
    }

    // In case the area or divergence variable is not initialized we initialize it
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
    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
            rNode.SetValue(*mpAreaVariable, 0.0);
            rNode.FastGetSolutionStepValue(*mpDivergenceVariable) = 0.0;
        });
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>::ClearDivergence()
{
    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
            rNode.SetValue(*mpAreaVariable, 0.0);
            rNode.SetValue(*mpDivergenceVariable, 0.0);
        });
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
    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
            rNode.FastGetSolutionStepValue(*mpDivergenceVariable) /=
                rNode.GetValue(*mpAreaVariable);
        });
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>::PonderateDivergence()
{
    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
            rNode.GetValue(*mpDivergenceVariable) /=
                rNode.GetValue(*mpAreaVariable);
        });
}

/***********************************************************************************/
/***********************************************************************************/

template class ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>;
template class ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>;

} /* namespace Kratos.*/