//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_ADJOINT_TEST_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_ADJOINT_TEST_UTILITIES_H_INCLUDED

// System includes
#include <functional>
#include <tuple>

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/geometrical_sensitivity_utility.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/test_utilities.h"

namespace Kratos
{
namespace RansApplicationTestUtilities
{
template <int I, class... TArgs>
typename std::enable_if<(I == sizeof...(TArgs)), void>::type
fill_value(
    std::tuple<TArgs...>& rSensitivity,
    const std::tuple<TArgs...>& rValue,
    const std::tuple<TArgs...>& rRefValue,
    const double Delta)
{
}

template <int I, class... TArgs>
typename std::enable_if<(I < sizeof...(TArgs)), void>::type
fill_value(
    std::tuple<TArgs...>& rSensitivity,
    const std::tuple<TArgs...>& rValue,
    const std::tuple<TArgs...>& rRefValue,
    const double Delta)
{
    std::get<I>(rSensitivity) = (std::get<I>(rValue) - std::get<I>(rRefValue)) / Delta;
    fill_value<I + 1, TArgs...>(rSensitivity, rValue, rRefValue, Delta);
}

template<class... TArgs>
std::tuple<TArgs...> ComputeFiniteDifferenceSensitivities(
    const std::tuple<TArgs...>& rValues,
    const std::tuple<TArgs...>& rRefValues,
    const double Delta)
{
    std::tuple<TArgs...> result;
    fill_value<0, TArgs...>(result, rValues, rRefValues, Delta);
    return result;
}

template<class TDataType>
void CompareValues(const TDataType& rA, const TDataType& rB, const double Tolerance);

template <int I, std::size_t TDerivativesSize, class... TArgs>
typename std::enable_if<(I == sizeof...(TArgs)), void>::type
compare_tuples(
    const std::tuple<TArgs...>& rFiniteDifferenceSensitivities,
    const std::tuple<BoundedVector<TArgs, TDerivativesSize>...>& rAdjointSensitivities,
    const unsigned int CheckIndex,
    const double RelativeTolerance)
{
}

template <int I, std::size_t TDerivativesSize, class... TArgs>
typename std::enable_if<(I < sizeof...(TArgs)), void>::type
compare_tuples(
    const std::tuple<TArgs...>& rFiniteDifferenceSensitivities,
    const std::tuple<BoundedVector<TArgs, TDerivativesSize>...>& rAdjointSensitivities,
    const unsigned int CheckIndex,
    const double RelativeTolerance)
{
    CompareValues(std::get<I>(rFiniteDifferenceSensitivities), std::get<I>(rAdjointSensitivities)[CheckIndex], RelativeTolerance);
    compare_tuples<I + 1, TDerivativesSize, TArgs...>(rFiniteDifferenceSensitivities, rAdjointSensitivities, CheckIndex, RelativeTolerance);
}


template<std::size_t TDerivativesSize, class... TArgs>
void CompareSensitivities(
    const std::tuple<TArgs...>& rFiniteDifferenceSensitivities,
    const std::tuple<BoundedVector<TArgs, TDerivativesSize>...>& rAdjointSensitivities,
    const unsigned int CheckIndex,
    const double RelativeTolerance)
{
    compare_tuples<0, TDerivativesSize, TArgs...>(rFiniteDifferenceSensitivities, rAdjointSensitivities, CheckIndex, RelativeTolerance);
}

template<class TDataType>

std::function<double&(ModelPart::NodeType&, const IndexType)> GetPerturbationMethod(
    const Variable<TDataType>& rPerturbationVariable);

template <class TPrimalElementDataType, class TAdjointElementDataType, class TAdjointDerivativeType, class... TArgs>
void RunAdjointDataTest(
    Model& rModel,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFunction,
    const std::function<void(ModelPart& rModelPart)>& rSetVariableDataFunction,
    const std::function<void(ModelPart& rModelPart)>& rUpdateFunction,
    const std::function<std::tuple<TArgs...>(const TPrimalElementDataType&, const Vector&, const Matrix&)>& rPrimalTupleRetrievalMethod,
    const std::function<std::tuple<BoundedVector<TArgs, TAdjointDerivativeType::TDerivativesSize>...>(const TAdjointDerivativeType&, const Vector&, const Matrix&, const IndexType)>& rAdjointTupleRetrievalMethod,
    const int BufferSize,
    const double Delta,
    const double RelativeTolerance)
{
    constexpr unsigned int derivatives_size = TAdjointDerivativeType::TDerivativesSize;
    constexpr unsigned int derivative_dimension = derivatives_size / 3;

    const auto& perturbation_method = GetPerturbationMethod(TAdjointDerivativeType::GetDerivativeVariable());

    // setup primal model part
    ModelPart& r_primal_model_part = CreateScalarVariableTestModelPart(
        rModel, "Element2D3N", "LineCondition2D2N", rAddNodalSolutionStepVariablesFunction,
        TPrimalElementDataType::GetScalarVariable(), BufferSize, false, false);
    rSetVariableDataFunction(r_primal_model_part);

    auto& r_primal_geometry = r_primal_model_part.Elements().front().GetGeometry();

    // setup primal element data
    TPrimalElementDataType primal_element_data(r_primal_geometry);

    // setup adjoint model part
    ModelPart& r_adjoint_model_part = CreateScalarVariableTestModelPart(
        rModel, "Element2D3N", "LineCondition2D2N", rAddNodalSolutionStepVariablesFunction,
        TAdjointElementDataType::GetAdjointScalarVariable(), BufferSize, false, false);
    rSetVariableDataFunction(r_adjoint_model_part);
    rUpdateFunction(r_adjoint_model_part);

    auto& r_adjoint_geometry = r_adjoint_model_part.Elements().front().GetGeometry();

    // setup adjoint element data
    TAdjointElementDataType adjoint_element_data(r_adjoint_geometry);
    TAdjointDerivativeType derivatives(adjoint_element_data);

    TPrimalElementDataType::Check(r_primal_geometry, r_primal_model_part.GetProcessInfo());
    TAdjointElementDataType::Check(r_adjoint_geometry, r_adjoint_model_part.GetProcessInfo());

    adjoint_element_data.CalculateConstants(r_adjoint_model_part.GetProcessInfo());
    primal_element_data.CalculateConstants(r_primal_model_part.GetProcessInfo());

    // calculate shape function values
    Vector adjoint_gauss_weights;
    Matrix adjoint_shape_functions;
    ModelPart::GeometryType::ShapeFunctionsGradientsType adjoint_shape_function_derivatives;
    RansCalculationUtilities::CalculateGeometryData(
        r_adjoint_geometry, GeometryData::GI_GAUSS_2,
        adjoint_gauss_weights, adjoint_shape_functions, adjoint_shape_function_derivatives);

    Vector primal_gauss_weights;
    Matrix primal_shape_functions;
    ModelPart::GeometryType::ShapeFunctionsGradientsType primal_shape_function_derivatives;

    BoundedMatrix<double, derivatives_size, 3> effective_velocity_derivatives;
    BoundedVector<double, derivatives_size> effective_kinematic_viscosity_derivatives, reaction_term_derivatives, source_term_derivatives;

    // calculate derivative values
    for (int g = 0; g < static_cast<int>(adjoint_gauss_weights.size()); ++g) {
        const Vector& adjoint_gauss_shape_functions = row(adjoint_shape_functions, g);
        const Matrix& adjoint_gauss_shape_function_derivatives =
            adjoint_shape_function_derivatives[g];

        // calculating adjoint sensitivities
        adjoint_element_data.CalculateGaussPointData(adjoint_gauss_shape_functions, adjoint_gauss_shape_function_derivatives);
        const auto& analytical_values =  rAdjointTupleRetrievalMethod(derivatives, adjoint_gauss_shape_functions, adjoint_gauss_shape_function_derivatives, g);

        // calculating primal finite difference sensitivities
        RansCalculationUtilities::CalculateGeometryData(
            r_primal_geometry, GeometryData::GI_GAUSS_2,
            primal_gauss_weights, primal_shape_functions, primal_shape_function_derivatives);

        rUpdateFunction(r_primal_model_part);

        const Vector& primal_gauss_shape_functions = row(primal_shape_functions, g);
        const Matrix& primal_gauss_shape_function_derivatives = primal_shape_function_derivatives[g];

        primal_element_data.CalculateGaussPointData(primal_gauss_shape_functions, primal_gauss_shape_function_derivatives);
        const auto& ref_values = rPrimalTupleRetrievalMethod(primal_element_data, primal_gauss_shape_functions, primal_gauss_shape_function_derivatives);

        for (IndexType i_node = 0; i_node < 3; ++i_node) {
            auto& r_node = r_primal_geometry[i_node];
            const auto block_row = i_node * derivative_dimension;
            for (IndexType k = 0; k < derivative_dimension; ++k) {
                perturbation_method(r_node, k) += Delta;

                RansCalculationUtilities::CalculateGeometryData(
                    r_primal_geometry, GeometryData::GI_GAUSS_2,
                    primal_gauss_weights, primal_shape_functions,
                    primal_shape_function_derivatives);

                const Vector& fd_gauss_shape_functions = row(primal_shape_functions, g);
                const Matrix& fd_gauss_shape_function_derivatives = primal_shape_function_derivatives[g];

                rUpdateFunction(r_primal_model_part);

                primal_element_data.CalculateGaussPointData(fd_gauss_shape_functions, fd_gauss_shape_function_derivatives);
                const auto& values = rPrimalTupleRetrievalMethod(primal_element_data, fd_gauss_shape_functions, fd_gauss_shape_function_derivatives);
                const auto& fd_values = ComputeFiniteDifferenceSensitivities(values, ref_values, Delta);
                CompareSensitivities(fd_values, analytical_values, block_row + k, RelativeTolerance);

                perturbation_method(r_node, k) -= Delta;
            }
        }
    }
}

template<class TPrimalElementDataType>
std::tuple<array_1d<double, 3>, double, double, double> RetrievePrimalValues(
    const TPrimalElementDataType& rElementData,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives)
{
    return std::make_tuple<array_1d<double, 3>, double, double, double>(
        rElementData.CalculateEffectiveVelocity(rShapeFunctions, rShapeFunctionDerivatives),
        rElementData.CalculateEffectiveKinematicViscosity(rShapeFunctions, rShapeFunctionDerivatives),
        rElementData.CalculateReactionTerm(rShapeFunctions, rShapeFunctionDerivatives),
        rElementData.CalculateSourceTerm(rShapeFunctions, rShapeFunctionDerivatives)
    );
}

template <class TPrimalElementDataType, class TAdjointElementDataType, class TAdjointDerivativeType>
void RunAdjointElementDataTest(
    Model& rModel,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFunction,
    const std::function<void(ModelPart& rModelPart)>& rSetVariableDataFunction,
    const std::function<void(ModelPart& rModelPart)>& rUpdateFunction,
    const int BufferSize,
    const double Delta,
    const double RelativeTolerance)
{
    using BArrN = BoundedVector<array_1d<double, 3>, TAdjointDerivativeType::TDerivativesSize>;
    using BDN = BoundedVector<double, TAdjointDerivativeType::TDerivativesSize>;

    const std::function<std::tuple<BArrN, BDN, BDN, BDN>(const TAdjointDerivativeType&, const Vector&, const Matrix&, const IndexType)>& adjoint_values = [](
        const TAdjointDerivativeType& rDerivatives,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const IndexType GaussPointIndex) -> std::tuple<BArrN, BDN, BDN, BDN> {

        BoundedMatrix<double, TAdjointDerivativeType::TDerivativesSize, 3> effective_velocity_derivatives;
        BoundedVector<double, TAdjointDerivativeType::TDerivativesSize> effective_kinematic_viscosity_derivatives, reaction_term_derivatives, source_term_derivatives;
        BoundedVector<array_1d<double, 3>, TAdjointDerivativeType::TDerivativesSize> effective_velocity_derivatives_arr;

        rDerivatives.CalculateEffectiveVelocityDerivatives(effective_velocity_derivatives, rShapeFunctions, rShapeFunctionDerivatives);
        rDerivatives.CalculateEffectiveKinematicViscosityDerivatives(effective_kinematic_viscosity_derivatives, rShapeFunctions, rShapeFunctionDerivatives);
        rDerivatives.CalculateReactionTermDerivatives(reaction_term_derivatives, rShapeFunctions, rShapeFunctionDerivatives);
        rDerivatives.CalculateSourceTermDerivatives(source_term_derivatives, rShapeFunctions, rShapeFunctionDerivatives);

        for (IndexType i = 0; i < TAdjointDerivativeType::TDerivativesSize; ++i) {
            noalias(effective_velocity_derivatives_arr[i]) =  row(effective_velocity_derivatives, i);
        }

        return std::make_tuple(
                effective_velocity_derivatives_arr,
                effective_kinematic_viscosity_derivatives,
                reaction_term_derivatives,
                source_term_derivatives
            );
    };

    RunAdjointDataTest<TPrimalElementDataType, TAdjointElementDataType, TAdjointDerivativeType, array_1d<double, 3>, double, double, double>(
        rModel,
        rAddNodalSolutionStepVariablesFunction,
        rSetVariableDataFunction,
        rUpdateFunction,
        RetrievePrimalValues<TPrimalElementDataType>,
        adjoint_values,
        BufferSize,
        Delta,
        RelativeTolerance
    );
}

template <class TPrimalElementDataType, class TAdjointElementDataType, class TAdjointDerivativeType>
void RunAdjointSensitivityDataTest(
    Model& rModel,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFunction,
    const std::function<void(ModelPart& rModelPart)>& rSetVariableDataFunction,
    const std::function<void(ModelPart& rModelPart)>& rUpdateFunction,
    const int BufferSize,
    const double Delta,
    const double RelativeTolerance)
{
    using BArrN = BoundedVector<array_1d<double, 3>, TAdjointDerivativeType::TDerivativesSize>;
    using BDN = BoundedVector<double, TAdjointDerivativeType::TDerivativesSize>;

    const std::function<std::tuple<BArrN, BDN, BDN, BDN>(
        const TAdjointDerivativeType&, const Vector&, const Matrix&, const IndexType)>& adjoint_values =
        [](const TAdjointDerivativeType& rDerivatives,
           const Vector& rShapeFunctions, const Matrix& rShapeFunctionDerivatives,
           const IndexType GaussPointIndex) -> std::tuple<BArrN, BDN, BDN, BDN> {

        BoundedVector<array_1d<double, 3>, TAdjointDerivativeType::TDerivativesSize> effective_velocity_derivatives;
        BoundedVector<double, TAdjointDerivativeType::TDerivativesSize> effective_kinematic_viscosity_derivatives,
            reaction_term_derivatives, source_term_derivatives;

        Geometry<Point>::JacobiansType J;
        rDerivatives.GetElementData().GetGeometry().Jacobian(J, GeometryData::GI_GAUSS_2);
        const auto& DN_De = rDerivatives.GetElementData().GetGeometry().ShapeFunctionsLocalGradients(GeometryData::GI_GAUSS_2);

        GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_deriv;
        const Matrix& rJ = J[GaussPointIndex];
        const Matrix& rDN_De = DN_De[GaussPointIndex];
        GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

        ShapeParameter deriv;
        for (int c = 0; c < 3; ++c) {
            const auto block_row = c * 2;
            for (int k = 0; k < 2; ++k) {
                deriv.NodeIndex = c;
                deriv.Direction = k;

                const auto value_row = block_row + k;

                double detJ_deriv;
                geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, DN_DX_deriv);
                effective_velocity_derivatives[value_row] = rDerivatives.CalculateEffectiveVelocityDerivatives(deriv, rShapeFunctions, rShapeFunctionDerivatives, detJ_deriv, DN_DX_deriv);
                effective_kinematic_viscosity_derivatives[value_row] = rDerivatives.CalculateEffectiveKinematicViscosityDerivatives(deriv, rShapeFunctions, rShapeFunctionDerivatives, detJ_deriv, DN_DX_deriv);
                reaction_term_derivatives[value_row] = rDerivatives.CalculateReactionTermDerivatives(deriv, rShapeFunctions, rShapeFunctionDerivatives, detJ_deriv, DN_DX_deriv);
                source_term_derivatives[value_row] = rDerivatives.CalculateSourceTermDerivatives(deriv, rShapeFunctions, rShapeFunctionDerivatives, detJ_deriv, DN_DX_deriv);
            }
        }

        return std::make_tuple(
                effective_velocity_derivatives,
                effective_kinematic_viscosity_derivatives,
                reaction_term_derivatives,
                source_term_derivatives
            );
    };

    RunAdjointDataTest<TPrimalElementDataType, TAdjointElementDataType, TAdjointDerivativeType, array_1d<double, 3>, double, double, double>(
        rModel,
        rAddNodalSolutionStepVariablesFunction,
        rSetVariableDataFunction,
        rUpdateFunction,
        RetrievePrimalValues<TPrimalElementDataType>,
        adjoint_values,
        BufferSize,
        Delta,
        RelativeTolerance
    );
}

} // namespace RansApplicationTestUtilities
} // namespace Kratos

#endif // KRATOS_RANS_ADJOINT_TEST_UTILITIES_H_INCLUDED