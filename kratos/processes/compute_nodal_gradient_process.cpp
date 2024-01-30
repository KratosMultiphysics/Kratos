//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "utilities/variable_utils.h"
#include "utilities/geometry_utilities.h"
#include "processes/compute_nodal_gradient_process.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
template<bool THistorical>
void ComputeNodalGradientProcess<THistorical>::Execute()
{
    KRATOS_TRY;

    // Set to zero
    ClearGradient();

    if (mrModelPart.NumberOfElements() != 0) {
        ComputeElementalContributionsAndVolume();
    }

    SynchronizeGradientAndVolume();

    PonderateGradient();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void ComputeNodalGradientProcess<THistorical>::ComputeElementalContributionsAndVolume() {

    // Auxiliar containers
    struct TLSType
    {
        Matrix DN_DX, J0, InvJ0;
        Vector N, values;
        double detJ0 = 0.0;
    };
    TLSType tls;

    // First element iterator
    const auto it_element_begin = mrModelPart.ElementsBegin();

    // Current domain size
    KRATOS_ERROR_IF_NOT(mrModelPart.GetProcessInfo().Has(DOMAIN_SIZE)) << "'DOMAIN_SIZE' is not found in '" << mrModelPart.FullName() << "' ProcessInfo container." << std::endl;
    const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // Initial resize
    const auto& r_first_element_geometry = it_element_begin->GetGeometry();
    const std::size_t number_of_nodes_first_element = r_first_element_geometry.PointsNumber();
    const std::size_t local_space_dimension_first_element = r_first_element_geometry.LocalSpaceDimension();
    if (tls.DN_DX.size1() != number_of_nodes_first_element || tls.DN_DX.size2() != dimension)
        tls.DN_DX.resize(number_of_nodes_first_element, dimension);
    if (tls.N.size() != number_of_nodes_first_element)
        tls.N.resize(number_of_nodes_first_element);
    if (tls.values.size() != number_of_nodes_first_element)
        tls.values.resize(number_of_nodes_first_element);
    if (tls.J0.size1() != dimension || tls.J0.size2() != local_space_dimension_first_element)
        tls.J0.resize(dimension, local_space_dimension_first_element);

    // Variable retriever
    Kratos::unique_ptr<AuxiliarVariableVectorRetriever> p_variable_retriever;
    if (mNonHistoricalVariable) {
        p_variable_retriever = Kratos::make_unique<VariableVectorRetriever<ComputeNodalGradientProcessSettings::GetAsNonHistoricalVariable>>();
    } else {
        p_variable_retriever = Kratos::make_unique<VariableVectorRetriever<ComputeNodalGradientProcessSettings::GetAsHistoricalVariable>>();
    }

    // Iterate over the elements
    block_for_each(mrModelPart.Elements(), tls, [&](Element& rElem, TLSType& rTls){
        auto& r_geometry = rElem.GetGeometry();

        // Current geometry information
        const std::size_t number_of_nodes = r_geometry.PointsNumber();

        // Resize if needed
        if (rTls.N.size() != number_of_nodes)
            rTls.N.resize(number_of_nodes);
        if (rTls.values.size() != number_of_nodes)
            rTls.values.resize(number_of_nodes);

        // The integration points
        const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
        const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
        const std::size_t number_of_integration_points = r_integration_points.size();

        // Fill vector
        p_variable_retriever->GetVariableVector(r_geometry, *mpOriginVariable, rTls.values);

        // The containers of the shape functions and the local gradients
        const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(r_integration_method);
        const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(r_integration_method);

        for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
            // Getting the shape functions
            noalias(rTls.N) = row(rNcontainer, point_number);

            // Getting the jacobians and local gradients
            GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[point_number], rTls.J0);
            MathUtils<double>::GeneralizedInvertMatrix(rTls.J0, rTls.InvJ0, rTls.detJ0);
            const Matrix& rDN_De = rDN_DeContainer[point_number];
            GeometryUtils::ShapeFunctionsGradients(rDN_De, rTls.InvJ0, rTls.DN_DX);

            const Vector grad = prod(trans(rTls.DN_DX), rTls.values);
            const double gauss_point_volume = r_integration_points[point_number].Weight() * (rTls.detJ0);

            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
                array_1d<double, 3>& r_gradient = GetGradient(r_geometry, i_node);
                for(std::size_t k=0; k<dimension; ++k) {
                    AtomicAdd(r_gradient[k], (rTls.N)[i_node] * gauss_point_volume*grad[k] );
                }

                double& vol = r_geometry[i_node].GetValue(*mpAreaVariable);

                AtomicAdd(vol, (rTls.N)[i_node] * gauss_point_volume );
            }
        }
    });

}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
ComputeNodalGradientProcess<THistorical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart,
    const Variable<double>& rOriginVariable,
    const Variable<array_1d<double,3> >& rGradientVariable,
    const Variable<double>& rAreaVariable,
    const bool NonHistoricalVariable
    ) : mrModelPart(rModelPart),
        mpOriginVariable(&rOriginVariable),
        mpGradientVariable(&rGradientVariable),
        mpAreaVariable(&rAreaVariable),
        mNonHistoricalVariable(NonHistoricalVariable)
{
    KRATOS_TRY

    CheckOriginAndAreaVariables();

    if constexpr (THistorical) {
        // Checking historical gradient variable
        VariableUtils().CheckVariableExists(rGradientVariable, mrModelPart.Nodes());
    } else {
        // In case the area or gradient variable is not initialized we initialize it
        if (mrModelPart.NumberOfNodes() != 0) {
            if (!mrModelPart.Nodes().begin()->Has( rGradientVariable )) {
                const array_1d<double,3> zero_vector = ZeroVector(3);
                VariableUtils().SetNonHistoricalVariable(rGradientVariable, zero_vector, mrModelPart.Nodes());
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
ComputeNodalGradientProcess<THistorical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : mrModelPart(rModelPart)
{
    KRATOS_TRY

    // We check the parameters
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // We get the gradient variable
    const std::string& r_origin_variable_name = ThisParameters["origin_variable"].GetString();
    const std::string& r_gradient_variable_name = ThisParameters["gradient_variable"].GetString();
    const std::string& r_area_variable_name = ThisParameters["area_variable"].GetString();

    // We push the list of double variables
    if (KratosComponents<Variable<double>>::Has(r_origin_variable_name) && KratosComponents<Variable<double>>::Has(r_area_variable_name)) {
        mpOriginVariable = &KratosComponents<Variable<double>>::Get(r_origin_variable_name);
        mpAreaVariable = &KratosComponents<Variable<double>>::Get(r_area_variable_name);
    } else {
        KRATOS_ERROR << "Only doubles are allowed as variables, given variables: " <<
            r_origin_variable_name << " " << r_area_variable_name << std::endl;
    }

    // We push the list of double variables
    if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_gradient_variable_name)) {
        mpGradientVariable = &KratosComponents<Variable<array_1d<double, 3>>>::Get(r_gradient_variable_name);
    } else {
        KRATOS_ERROR << "Only vectors are allowed as variables, given variable: " <<
            r_gradient_variable_name << std::endl;
    }

    // Setting the non-historical flag
    mNonHistoricalVariable = ThisParameters["non_historical_origin_variable"].GetBool();

    ComputeNodalGradientProcess<THistorical>(mrModelPart,
        *mpOriginVariable,
        *mpGradientVariable,
        *mpAreaVariable,
        mNonHistoricalVariable);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void ComputeNodalGradientProcess<THistorical>::CheckOriginAndAreaVariables()
{
    KRATOS_TRY

    auto& r_nodes = mrModelPart.Nodes();
    // Doing several checks
    if (!mNonHistoricalVariable) {
        VariableUtils().CheckVariableExists(*mpOriginVariable, r_nodes);
    } else {
        bool hasVariable = 0;
        if(mrModelPart.NumberOfNodes() != 0)
            hasVariable = r_nodes.begin()->Has(*mpOriginVariable);
        hasVariable = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(hasVariable);

        KRATOS_ERROR_IF_NOT(hasVariable) << "Variable " << mpOriginVariable->Name() << " not defined on non-historial database" << std::endl;
    }

    // In case the area or gradient variable is not initialized we initialize it
    if (mrModelPart.NumberOfNodes() != 0) {
        if (!r_nodes.begin()->Has(*mpAreaVariable)) {
            VariableUtils().SetNonHistoricalVariable(*mpAreaVariable, 0.0, r_nodes);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void ComputeNodalGradientProcess<THistorical>::ClearGradient()
{
    if constexpr (THistorical) {
        block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
            rNode.SetValue(*mpAreaVariable, 0.0);
            rNode.FastGetSolutionStepValue(*mpGradientVariable).clear();
        });
    } else {
        const array_1d<double, 3> aux_zero_vector = ZeroVector(3);
        block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
            rNode.SetValue(*mpAreaVariable, 0.0);
            rNode.SetValue(*mpGradientVariable, aux_zero_vector);
        });
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
array_1d<double, 3>& ComputeNodalGradientProcess<THistorical>::GetGradient(
    Element::GeometryType& rThisGeometry,
    unsigned int i
    )
{
    if constexpr (THistorical) {
        return rThisGeometry[i].FastGetSolutionStepValue(*mpGradientVariable);
    } else {
        return rThisGeometry[i].GetValue(*mpGradientVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void ComputeNodalGradientProcess<THistorical>::PonderateGradient()
{
    if constexpr (THistorical) {
        block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
            rNode.FastGetSolutionStepValue(*mpGradientVariable) /=
            rNode.GetValue(*mpAreaVariable);
        });
    } else {
        block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
            rNode.GetValue(*mpGradientVariable) /=
            rNode.GetValue(*mpAreaVariable);
        });
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void ComputeNodalGradientProcess<THistorical>::SynchronizeGradientAndVolume()
{
    if constexpr (THistorical) {
        mrModelPart.GetCommunicator().AssembleCurrentData(*mpGradientVariable);
        mrModelPart.GetCommunicator().AssembleNonHistoricalData(*mpAreaVariable);
    } else {
        mrModelPart.GetCommunicator().AssembleNonHistoricalData(*mpGradientVariable);
        mrModelPart.GetCommunicator().AssembleNonHistoricalData(*mpAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void VariableVectorRetriever<THistorical>::GetVariableVector(
    const Geometry<Node>& rGeometry,
    const Variable<double>& rVariable,
    Vector& rVector
    )
{
    if constexpr (THistorical) {
        for(std::size_t i_node=0; i_node < rGeometry.size(); ++i_node) {
            rVector[i_node] = rGeometry[i_node].FastGetSolutionStepValue(rVariable);
        }
    } else {
        for(std::size_t i_node=0; i_node < rGeometry.size(); ++i_node) {
            rVector[i_node] = rGeometry[i_node].GetValue(rVariable);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
const Parameters ComputeNodalGradientProcess<THistorical>::GetDefaultParameters() const
{
    Parameters default_parameters = Parameters(R"(
    {
        "origin_variable"                : "PLEASE_DEFINE_A_VARIABLE",
        "gradient_variable"              : "PLEASE_DEFINE_A_VARIABLE",
        "area_variable"                  : "NODAL_AREA",
        "non_historical_origin_variable" :  false
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template struct VariableVectorRetriever<ComputeNodalGradientProcessSettings::GetAsHistoricalVariable>;
template struct VariableVectorRetriever<ComputeNodalGradientProcessSettings::GetAsNonHistoricalVariable>;

/***********************************************************************************/
/***********************************************************************************/

template class ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>;
template class ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>;

} /* namespace Kratos.*/
