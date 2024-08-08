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

#if !defined(KRATOS_ADJOINT_TEST_UTILITIES_H_INCLUDED)
#define KRATOS_ADJOINT_TEST_UTILITIES_H_INCLUDED

// System includes
#include <cmath>
#include <sstream>
#include <iomanip>
#include <functional>
#include <tuple>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_utilities/fluid_test_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"
#include "custom_utilities/test_utilities.h"


/**
 * @brief Based on python version of is_close
 *
 * The core KRATOS_CHECK_RELATIVE_NEAR fails when your first argument is zero and
 * second argument is close to zero (but more than the numeric_limits<double>::epsilon()).
 * So then the check fails. There are situations where I cannot use NEAR_CHECKS because this
 * testing framework is used in all element data derivative testing, and some of them has to check
 * large numbers so KRATOS_CHECK_NEAR will fail at those test. So this hybrid version is implemented
 * by following the python is_close.
 *
 * Based on the following formula:
 *      abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
 *
 **/
#define KRATOS_CHECK_IS_CLOSE(a, b, relative_tolerance, absolute_tolerance)   \
    {                                                                         \
        if (std::abs(a - b) >                                                 \
            std::max(relative_tolerance * std::max(std::abs(a), std::abs(b)), \
                     absolute_tolerance)) {                                   \
            std::stringstream msg;                                            \
            msg << std::scientific << std::setprecision(16);                  \
            msg << "Check failed because " << #a << " = " << a                \
                << " is not close to " << #b << " = " << b                    \
                << " within the relative tolerance " << relative_tolerance    \
                << " or absolute tolerance " << absolute_tolerance;           \
            KRATOS_ERROR << msg.str();                                        \
        }                                                                     \
    }

/**
 * @brief Based on python version of is_close
 *
 * The core KRATOS_CHECK_VECTOR_RELATIVE_NEAR fails when your first argument is zero and
 * second argument is close to zero (but more than the numeric_limits<double>::epsilon()).
 * So then the check fails. There are situations where I cannot use NEAR_CHECKS because this
 * testing framework is used in all element data derivative testing, and some of them has to check
 * large numbers so KRATOS_CHECK_VECTOR_NEAR will fail at those test. So this hybrid version is implemented
 * by following the python is_close.
 *
 * Based on the following formula:
 *      abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
 *
 **/
#define KRATOS_CHECK_VECTOR_IS_CLOSE(a, b, relative_tolerance, absolute_tolerance)      \
    {                                                                                   \
        KRATOS_ERROR_IF_NOT(a.size() == b.size())                                       \
            << "Check failed because vector arguments do not have the same "            \
               "size:"                                                                  \
            << std::endl                                                                \
            << "First argument has size " << a.size() << ", "                           \
            << "second argument has size " << b.size() << "." << std::endl;             \
        for (std::size_t i = 0; i < a.size(); i++) {                                    \
            if (std::abs(a[i] - b[i]) >                                                 \
                std::max(relative_tolerance * std::max(std::abs(a[i]), std::abs(b[i])), \
                         absolute_tolerance)) {                                         \
                std::stringstream msg;                                                  \
                msg << std::scientific << std::setprecision(16);                        \
                msg << "Check failed because " << #a << "[" << i                        \
                    << "]  = " << a[i] << " is not close to " << #b << "["              \
                    << i << "] = " << b[i] << " within the relative tolerance "         \
                    << relative_tolerance << " or absolute tolerance "                  \
                    << absolute_tolerance << "." << std::endl;                          \
                msg << std::endl << "Vector " << #a << " values:" << std::endl;         \
                for (std::size_t i = 0; i < a.size(); ++i) {                            \
                    msg << #a << "[" << i << "] = " << a[i] << std::endl;               \
                }                                                                       \
                msg << std::endl << "Vector " << #b << " values:" << std::endl;         \
                for (std::size_t i = 0; i < b.size(); ++i) {                            \
                    msg << #b << "[" << i << "] = " << b[i] << std::endl;               \
                }                                                                       \
                KRATOS_ERROR << msg.str();                                              \
            }                                                                           \
        }                                                                               \
    }

namespace Kratos
{
class AdjointTestUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    using GeometryType = ModelPart::GeometryType;

    using PropertiesType = ModelPart::PropertiesType;

    using ShapeFunctionsGradientsType = GeometryType::ShapeFunctionsGradientsType;

    using IndexType = std::size_t;

    using ArrayD = array_1d<double, 2>;

    ///@}
    ///@name Static Operations
    ///@{

    template<class TPrimalElementDataType, class TAdjointElementDerivativeDataType>
    static void RunElementDataDerivativeTest(
        Model& rModel,
        const std::function<void(ModelPart&)>& rAddNodalSolutionStepVariablesFunction,
        const std::function<void(ModelPart&)>& rSetVariableDataFunction,
        const std::function<void(Properties&)>& rSetProperties,
        const std::function<void(ModelPart&)>& rUpdateFunction,
        const IndexType BufferSize,
        const double Perturbation,
        const double RelativeTolerance,
        const double AbsoluteTolerance = -1.0)
    {
        const double absolute_tolerance = (AbsoluteTolerance == -1.0) ? RelativeTolerance * 1e-3 : AbsoluteTolerance;
        using TAdjointElementDataType = typename TAdjointElementDerivativeDataType::DataType;
        constexpr IndexType derivative_dimension = TAdjointElementDerivativeDataType::TDerivativeDimension;

        const auto& add_dofs = [&](NodeType& rNode) {
            rNode.AddDof(TPrimalElementDataType::GetScalarVariable());
            rNode.AddDof(TAdjointElementDataType::GetAdjointScalarVariable());
        };

        // setup primal model part
        auto& r_model_part = FluidTestUtilities::CreateTestModelPart(
            rModel, "test", "Element2D3N", "LineCondition2D2N", rSetProperties,
            [](PropertiesType&) {}, rAddNodalSolutionStepVariablesFunction,
            add_dofs, BufferSize);
        rSetVariableDataFunction(r_model_part);

        RansVariableUtilities::SetElementConstitutiveLaws(r_model_part.Elements());

        auto& r_process_info = r_model_part.GetProcessInfo();
        auto& r_element = r_model_part.Elements().front();
        auto& r_geometry = r_element.GetGeometry();
        auto& r_properties = r_element.GetProperties();

        r_element.Initialize(r_process_info);

        const auto& r_integration_method = r_element.GetIntegrationMethod();

        // setup primal element data
        TPrimalElementDataType primal_element_data(r_geometry, r_properties, r_process_info);

        // setup adjoint element data
        TAdjointElementDataType adjoint_element_data(r_geometry, r_properties, r_process_info);

        TPrimalElementDataType::Check(r_element, r_process_info);
        TAdjointElementDataType::Check(r_element, r_process_info);

        adjoint_element_data.CalculateConstants(r_process_info);
        primal_element_data.CalculateConstants(r_process_info);

        // calculate shape function values
        Vector adjoint_W_values;
        Matrix adjoint_N_vectors;
        ShapeFunctionsGradientsType adjoint_dNdX_matrices;
        RansCalculationUtilities::CalculateGeometryData(
            r_geometry, r_integration_method,
            adjoint_W_values, adjoint_N_vectors, adjoint_dNdX_matrices);

        Vector primal_W_values;
        Matrix primal_N_vectors;
        ShapeFunctionsGradientsType primal_dNdX_matrices;

        ShapeParameter deriv;

        // calculate derivative values
        for (IndexType g = 0; g < adjoint_W_values.size(); ++g) {
            // calculate adjoint shape function data
            const double W = adjoint_W_values[g];
            const Vector& adjoint_N = row(adjoint_N_vectors, g);
            const Matrix& adjoint_dNdX = adjoint_dNdX_matrices[g];

            // calculate reference values for finite difference derivative computation
            rUpdateFunction(r_model_part);

            ArrayD ref_velocity, velocity, adjoint_velocity;
            double ref_viscosity, ref_reaction_term, ref_source_term, viscosity,
                reaction_term, source_term, adjoint_viscosity,
                adjoint_reaction_term, adjoint_source_term;
            std::tie(ref_velocity, ref_viscosity, ref_reaction_term, ref_source_term) =
                ComputeElementDataValues(r_geometry, g, primal_W_values,
                                         primal_N_vectors, primal_dNdX_matrices,
                                         r_process_info, r_integration_method, primal_element_data);

            std::tie(adjoint_velocity, adjoint_viscosity, adjoint_reaction_term,
                     adjoint_source_term) =
                ComputeElementDataValues(r_geometry, g, adjoint_W_values,
                                         adjoint_N_vectors, adjoint_dNdX_matrices,
                                         r_process_info, r_integration_method, adjoint_element_data);

            // check whether primal and derivatives have the same coefficients
            KRATOS_CHECK_VECTOR_IS_CLOSE(ref_velocity, adjoint_velocity, RelativeTolerance, absolute_tolerance);
            KRATOS_CHECK_IS_CLOSE(ref_viscosity, adjoint_viscosity, RelativeTolerance, absolute_tolerance);
            KRATOS_CHECK_IS_CLOSE(ref_reaction_term, adjoint_reaction_term, RelativeTolerance, absolute_tolerance);
            KRATOS_CHECK_IS_CLOSE(ref_source_term, adjoint_source_term, RelativeTolerance, absolute_tolerance);

            GeometryType::JacobiansType J;
            r_geometry.Jacobian(J, r_integration_method);
            const auto& DN_De = r_geometry.ShapeFunctionsLocalGradients(r_integration_method);

            GeometricalSensitivityUtility::ShapeFunctionsGradientType adjoint_dNdX_derivative;
            const Matrix& rJ = J[g];
            const Matrix& rDN_De = DN_De[g];
            const double inv_detJ = 1.0 / MathUtils<double>::Det(rJ);
            GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

            for (IndexType c = 0; c < 3; ++c) {
                for (IndexType k = 0; k < derivative_dimension; ++k) {
                    deriv.NodeIndex = c;
                    deriv.Direction = k;

                    double adjoint_det_j_derivative;
                    geom_sensitivity.CalculateSensitivity(deriv, adjoint_det_j_derivative, adjoint_dNdX_derivative);
                    const double adjoint_W_derivative = adjoint_det_j_derivative * inv_detJ * W;

                    TAdjointElementDerivativeDataType derivative_data(
                        adjoint_element_data, c, k, W, adjoint_N, adjoint_dNdX,
                        adjoint_W_derivative, adjoint_det_j_derivative,
                        adjoint_dNdX_derivative);

                    // compute effective velocity analytical derivative
                    const auto& analytical_velocity_derivative = derivative_data.CalculateEffectiveVelocityDerivative();

                    // compute effective kinematic analytical derivative
                    const auto& analytical_viscosity_derivative = derivative_data.CalculateEffectiveKinematicViscosityDerivative();

                    // compute reaction term analytical derivative
                    const auto& analytical_reaction_term_derivative = derivative_data.CalculateReactionTermDerivative();

                    // compute source term analytical derivative
                    const auto& analytical_source_term_derivative = derivative_data.CalculateSourceTermDerivative();

                    const auto& derivative_variable = derivative_data.GetDerivativeVariable();

                    // finite difference derivative computations
                    PerturbNodalVariable(r_geometry[c], derivative_variable, Perturbation);

                    rUpdateFunction(r_model_part);

                    TPrimalElementDataType perturbed_primal_element_data(r_geometry, r_properties, r_process_info);
                    perturbed_primal_element_data.CalculateConstants(r_process_info);

                    std::tie(velocity, viscosity, reaction_term, source_term) =
                        ComputeElementDataValues(
                            r_geometry, g, primal_W_values, primal_N_vectors,
                            primal_dNdX_matrices, r_process_info, r_integration_method, perturbed_primal_element_data);

                    const ArrayD fd_velocity_derivative = (velocity - ref_velocity) / Perturbation;
                    const double fd_viscosity_derivative = (viscosity - ref_viscosity) / Perturbation;
                    const double fd_reaction_term_derivative = (reaction_term - ref_reaction_term) / Perturbation;
                    const double fd_source_term_derivative = (source_term - ref_source_term) / Perturbation;

                    PerturbNodalVariable(r_geometry[c], derivative_variable, -Perturbation);

                    // compare finite difference values with analytical derivatives
                    KRATOS_CHECK_VECTOR_IS_CLOSE(fd_velocity_derivative, analytical_velocity_derivative, RelativeTolerance, absolute_tolerance);
                    KRATOS_CHECK_IS_CLOSE(fd_viscosity_derivative, analytical_viscosity_derivative, RelativeTolerance, absolute_tolerance);
                    KRATOS_CHECK_IS_CLOSE(fd_reaction_term_derivative, analytical_reaction_term_derivative, RelativeTolerance, absolute_tolerance);
                    KRATOS_CHECK_IS_CLOSE(fd_source_term_derivative, analytical_source_term_derivative, RelativeTolerance, absolute_tolerance);
                }
            }
        }
    }

    ///@}

private:
    ///@name Private Static Operations
    ///@{

    static void PerturbNodalVariable(
        NodeType& rNode,
        const Variable<double>& rVariable,
        const double Perturbation);

    template <class TElementDataType>
    static std::tuple<ArrayD, double, double, double> ComputeElementDataValues(
        const GeometryType& rGeometry,
        const IndexType GaussPointIndex,
        Vector& rWs,
        Matrix& rNs,
        ShapeFunctionsGradientsType& rdNdXs,
        const ProcessInfo& rProcessInfo,
        const GeometryData::IntegrationMethod& rIntegrationMethod,
        TElementDataType& rElementData)
    {
        RansCalculationUtilities::CalculateGeometryData(
            rGeometry, rIntegrationMethod, rWs, rNs, rdNdXs);
        const Vector& N = row(rNs, GaussPointIndex);
        const Matrix& dNdX = rdNdXs[GaussPointIndex];

        rElementData.CalculateConstants(rProcessInfo);
        rElementData.CalculateGaussPointData(N, dNdX);

        return std::make_tuple(rElementData.GetEffectiveVelocity(),
                               rElementData.GetEffectiveKinematicViscosity(),
                               rElementData.GetReactionTerm(),
                               rElementData.GetSourceTerm());
    }

    ///@}

};
} // namespace Kratos

#endif // KRATOS_ADJOINT_TEST_UTILITIES_H_INCLUDED