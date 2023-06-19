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

/* Sysytem includes */
#include <functional>

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
    struct TLSType
    {
        Matrix J0, InvJ0, DN_DX;
    };
    TLSType tls;

    // First element iterator
    const auto it_element_begin = mrModelPart.ElementsBegin();

    // Current domain size
    const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // Initial resize
    const auto& r_first_element_geometry = it_element_begin->GetGeometry();
    const std::size_t number_of_nodes_first_element = r_first_element_geometry.PointsNumber();
    const std::size_t local_space_dimension_first_element = r_first_element_geometry.LocalSpaceDimension();
    if (tls.DN_DX.size1() != number_of_nodes_first_element || tls.DN_DX.size2() != dimension)
        tls.DN_DX.resize(number_of_nodes_first_element, dimension);
    if (tls.J0.size1() != dimension || tls.J0.size2() != local_space_dimension_first_element)
        tls.J0.resize(dimension, local_space_dimension_first_element);

    std::function<array_1d<double,3>(const Node&, const Variable<array_1d<double,3>>&)> get_vector_field;
    if (mNonHistoricalOriginVariable) {
        if (mNormalizeDivergence) {
            get_vector_field = GetNonHistoricalNormalVectorField;
        } else {
            get_vector_field = GetNonHistoricalVectorField;
        }
    } else {
        if (mNormalizeDivergence) {
            get_vector_field = GetHistoricalNormalVectorField;
        } else {
            get_vector_field = GetHistoricalVectorField;
        }
    }

    // Iterate over the elements
    block_for_each(mrModelPart.Elements(), tls, [&](Element& rElem, TLSType& rTls){
        auto& r_geometry = rElem.GetGeometry();

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
            GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[point_number], rTls.J0);
            MathUtils<double>::GeneralizedInvertMatrix(rTls.J0, rTls.InvJ0, detJ0);
            const Matrix& rDN_De = rDN_DeContainer[point_number];
            GeometryUtils::ShapeFunctionsGradients(rDN_De, rTls.InvJ0, rTls.DN_DX);

            double divergence = 0.0;
            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
                const auto& vector_field = get_vector_field(r_geometry[i_node], *mpOriginVariable);

                for(std::size_t k=0; k<dimension; ++k) {
                    divergence += rTls.DN_DX(i_node, k)*vector_field[k];
                }
            }

            const double gauss_point_volume = r_integration_points[point_number].Weight() * detJ0;

            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
                double& r_divergence = GetDivergence(r_geometry, i_node);

                AtomicAdd(r_divergence, N[i_node] * gauss_point_volume * divergence );

                double& vol = r_geometry[i_node].GetValue(*mpAreaVariable);

                AtomicAdd(vol, N[i_node] * gauss_point_volume );
            }
        }
    });

    SynchronizeDivergenceAndVolume();

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
    const bool NormalizeDivergence,
    const bool NonHistoricalOriginVariable
    ) : mrModelPart(rModelPart),
        mpOriginVariable(&rOriginVariable),
        mpDivergenceVariable(&rDivergenceVariable),
        mpAreaVariable(&rAreaVariable),
        mNormalizeDivergence(NormalizeDivergence),
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
    const bool NormalizeDivergence,
    const bool NonHistoricalOriginVariable
    ) : mrModelPart(rModelPart),
        mpOriginVariable(&rOriginVariable),
        mpDivergenceVariable(&rDivergenceVariable),
        mpAreaVariable(&rAreaVariable),
        mNormalizeDivergence(NormalizeDivergence),
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
    block_for_each(mrModelPart.Nodes(), [&](Node& rNode){
            rNode.SetValue(*mpAreaVariable, 0.0);
            rNode.FastGetSolutionStepValue(*mpDivergenceVariable) = 0.0;
        });
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>::ClearDivergence()
{
    block_for_each(mrModelPart.Nodes(), [&](Node& rNode){
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
void ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>::SynchronizeDivergenceAndVolume()
{
    mrModelPart.GetCommunicator().AssembleCurrentData(*mpDivergenceVariable);
    mrModelPart.GetCommunicator().AssembleNonHistoricalData(*mpAreaVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>::SynchronizeDivergenceAndVolume()
{
    mrModelPart.GetCommunicator().AssembleNonHistoricalData(*mpDivergenceVariable);
    mrModelPart.GetCommunicator().AssembleNonHistoricalData(*mpAreaVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>::PonderateDivergence()
{
    block_for_each(mrModelPart.Nodes(), [&](Node& rNode){
            rNode.FastGetSolutionStepValue(*mpDivergenceVariable) /=
                rNode.GetValue(*mpAreaVariable);
        });
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>::PonderateDivergence()
{
    block_for_each(mrModelPart.Nodes(), [&](Node& rNode){
            rNode.GetValue(*mpDivergenceVariable) /=
                rNode.GetValue(*mpAreaVariable);
        });
}

/***********************************************************************************/
/***********************************************************************************/

template <bool THistorical>
array_1d<double,3> ComputeNodalNormalDivergenceProcess<THistorical>::GetHistoricalNormalVectorField(
    const Node& rNode,
    const Variable<array_1d<double,3>>& rVariable)
{
    const auto& vector_field = rNode.FastGetSolutionStepValue(rVariable);
    const double norm = norm_2(vector_field);
    #ifdef KRATOS_DEBUG
    KRATOS_WARNING_IF("NodalNormalDivergenceProcess", norm < 1.0e-12) << "Unexpected zero norm" <<std::endl;
    #endif
    return vector_field / norm;
}

template <bool THistorical>
array_1d<double,3> ComputeNodalNormalDivergenceProcess<THistorical>::GetNonHistoricalNormalVectorField(
    const Node& rNode,
    const Variable<array_1d<double,3>>& rVariable)
{
    const auto& vector_field = rNode.GetValue(rVariable);
    const double norm = norm_2(vector_field);
    #ifdef KRATOS_DEBUG
    KRATOS_WARNING_IF("NodalNormalDivergenceProcess", norm < 1.0e-12) << "Unexpected zero norm" <<std::endl;
    #endif
    return vector_field / norm;
}

template <bool THistorical>
array_1d<double,3> ComputeNodalNormalDivergenceProcess<THistorical>::GetHistoricalVectorField(
    const Node& rNode,
    const Variable<array_1d<double,3>>& rVariable)
{
    return rNode.FastGetSolutionStepValue(rVariable);
}

template <bool THistorical>
array_1d<double,3> ComputeNodalNormalDivergenceProcess<THistorical>::GetNonHistoricalVectorField(
    const Node& rNode,
    const Variable<array_1d<double,3>>& rVariable)
{
    return rNode.GetValue(rVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template class ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>;
template class ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>;

} /* namespace Kratos.*/
