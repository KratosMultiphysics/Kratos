//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

#if !defined(KRATOS_RANS_ADJOINT_TEST_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_ADJOINT_TEST_UTILITIES_H_INCLUDED

// System includes
#include <functional>

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "test_utilities.h"

namespace Kratos
{
namespace RansApplicationTestUtilities
{
using NodeType = Node<3>;
using ElementType = ModelPart::ElementType;

/**
 * @brief based on python's math.IsNear()
 */
bool IsNear(const double ValueA, const double ValueB, const double RelTol, const double AbsTol);

void CheckNear(const double ValueA, const double ValueB, const double RelTol, const double AbsTol);

void CheckNear(const Matrix& rA, const Matrix& rB, const double RelTol, const double AbsTol);

void CheckNear(const Vector& rA, const Vector& rB, const double RelTol, const double AbsTol);

template <class TPrimalElementDataType, class TAdjointElementDataType>
void RunAdjointElementDataTest(
    Model& rModel,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFunction,
    const std::function<void(ModelPart& rModelPart)>& rSetVariableDataFunction,
    const std::function<void(ModelPart& rModelPart)>& rUpdateFunction,
    const std::function<double(const TPrimalElementDataType&, const Vector&, const Matrix&)>& rPrimalValueFunction,
    const std::function<double&(NodeType&)>& rPerturbVariableFunction,
    const std::function<void(BoundedVector<double, 3>&, const TAdjointElementDataType&, const Vector&, const Matrix&, const int)>& rAdjointDerivativesFunction,
    const int BufferSize,
    const double Delta,
    const double RelativeTolerance,
    const double AbsoluteTolerance)
{
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

    // check for integration method
    KRATOS_CHECK_EQUAL(TPrimalElementDataType::GetIntegrationMethod(),
                       TAdjointElementDataType::GetIntegrationMethod())
        << "Integration type mismatch between " << primal_element_data.GetName()
        << " and " << adjoint_element_data.GetName() << ".";

    TPrimalElementDataType::Check(r_primal_geometry, r_primal_model_part.GetProcessInfo());
    TAdjointElementDataType::Check(r_adjoint_geometry,
                                   r_adjoint_model_part.GetProcessInfo());

    adjoint_element_data.CalculateConstants(r_adjoint_model_part.GetProcessInfo());
    primal_element_data.CalculateConstants(r_primal_model_part.GetProcessInfo());

    // calculate shape function values
    Vector adjoint_gauss_weights;
    Matrix adjoint_shape_functions;
    ModelPart::GeometryType::ShapeFunctionsGradientsType adjoint_shape_function_derivatives;
    RansCalculationUtilities::CalculateGeometryData(
        r_adjoint_geometry, TPrimalElementDataType::GetIntegrationMethod(),
        adjoint_gauss_weights, adjoint_shape_functions, adjoint_shape_function_derivatives);

    Vector primal_gauss_weights;
    Matrix primal_shape_functions;
    ModelPart::GeometryType::ShapeFunctionsGradientsType primal_shape_function_derivatives;

    // calculate derivative values
    for (int g = 0; g < static_cast<int>(adjoint_gauss_weights.size()); ++g)
    {
        const Vector& adjoint_gauss_shape_functions = row(adjoint_shape_functions, g);
        const Matrix& adjoint_gauss_shape_function_derivatives =
            adjoint_shape_function_derivatives[g];

        // calculating adjoint sensitivities
        adjoint_element_data.CalculateGaussPointData(
            adjoint_gauss_shape_functions, adjoint_gauss_shape_function_derivatives);
        BoundedVector<double, 3> adjoint_derivatives;
        rAdjointDerivativesFunction(adjoint_derivatives, adjoint_element_data,
                                    adjoint_gauss_shape_functions,
                                    adjoint_gauss_shape_function_derivatives, g);

        // calculating primal finite difference sensitivities
        RansCalculationUtilities::CalculateGeometryData(
            r_primal_geometry, TPrimalElementDataType::GetIntegrationMethod(),
            primal_gauss_weights, primal_shape_functions, primal_shape_function_derivatives);

        rUpdateFunction(r_primal_model_part);

        const Vector& primal_gauss_shape_functions = row(primal_shape_functions, g);
        const Matrix& primal_gauss_shape_function_derivatives =
            primal_shape_function_derivatives[g];

        primal_element_data.CalculateGaussPointData(
            primal_gauss_shape_functions, primal_gauss_shape_function_derivatives);

        const double ref_value =
            rPrimalValueFunction(primal_element_data, primal_gauss_shape_functions,
                                 primal_gauss_shape_function_derivatives);

        BoundedVector<double, 3> finite_difference_derivatives;
        for (int i_node = 0; i_node < 3; ++i_node)
        {
            NodeType& r_node = r_primal_geometry[i_node];
            rPerturbVariableFunction(r_node) += Delta;

            RansCalculationUtilities::CalculateGeometryData(
                r_primal_geometry, TPrimalElementDataType::GetIntegrationMethod(),
                primal_gauss_weights, primal_shape_functions,
                primal_shape_function_derivatives);

            const Vector& fd_gauss_shape_functions = row(primal_shape_functions, g);
            const Matrix& fd_gauss_shape_function_derivatives =
                primal_shape_function_derivatives[g];

            rUpdateFunction(r_primal_model_part);

            primal_element_data.CalculateGaussPointData(
                fd_gauss_shape_functions, fd_gauss_shape_function_derivatives);

            const double perturbed_value =
                rPrimalValueFunction(primal_element_data, fd_gauss_shape_functions,
                                     fd_gauss_shape_function_derivatives);
            finite_difference_derivatives[i_node] = (perturbed_value - ref_value) / Delta;

            rPerturbVariableFunction(r_node) -= Delta;
        }

        CheckNear(finite_difference_derivatives, adjoint_derivatives,
                  RelativeTolerance, AbsoluteTolerance);
    }
}

template <class TPrimalElementDataType, class TAdjointElementDataType>
void RunAdjointElementDataTest(
    Model& rModel,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFunction,
    const std::function<void(ModelPart& rModelPart)>& rSetVariableDataFunction,
    const std::function<void(ModelPart& rModelPart)>& rUpdateFunction,
    const std::function<double(const TPrimalElementDataType&, const Vector&, const Matrix&)>& rPrimalValueFunction,
    const std::function<double&(NodeType&, const int)>& rPerturbVariableFunction,
    const std::function<void(BoundedMatrix<double, 3, 2>&, const TAdjointElementDataType&, const Vector&, const Matrix&, const int)>& rAdjointDerivativesFunction,
    const int BufferSize,
    const double Delta,
    const double RelativeTolerance,
    const double AbsoluteTolerance)
{
    const int domain_size = 2;

    for (int i_dim = 0; i_dim < domain_size; ++i_dim)
    {
        const auto& adjoint_derivatives_func =
            [&rAdjointDerivativesFunction, i_dim](
                BoundedVector<double, 3>& rOutput,
                const TAdjointElementDataType& rElementData, const Vector& rShapeFunctions,
                const Matrix& rShapeFunctionDerivatives, const int GaussPointIndex) {
                BoundedMatrix<double, 3, 2> output_matrix;
                rAdjointDerivativesFunction(output_matrix, rElementData, rShapeFunctions,
                                            rShapeFunctionDerivatives, GaussPointIndex);
                noalias(rOutput) = column(output_matrix, i_dim);
            };

        const auto& perturb_variable = [&rPerturbVariableFunction,
                                        i_dim](NodeType& rNode) -> double& {
            return rPerturbVariableFunction(rNode, i_dim);
        };

        RunAdjointElementDataTest<TPrimalElementDataType, TAdjointElementDataType>(
            rModel, rAddNodalSolutionStepVariablesFunction, rSetVariableDataFunction,
            rUpdateFunction, rPrimalValueFunction, perturb_variable, adjoint_derivatives_func,
            BufferSize, Delta, RelativeTolerance, AbsoluteTolerance);
    }
}

std::function<double&(NodeType&)> GetPerturbationMethod(const Variable<double>& rPerturbationVariable);

std::function<double&(NodeType&, const int)> GetPerturbationMethod(
    const Variable<array_1d<double, 3>>& rPerturbationVariable);

template <class TPrimalElementDataType, class TAdjointElementDataType>
void RunAdjointElementDataTest(
    Model& rModel,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFunction,
    const std::function<void(ModelPart& rModelPart)>& rSetVariableDataFunction,
    const std::function<void(ModelPart& rModelPart)>& rUpdateFunction,
    const Variable<double>& rDerivativeVariable,
    double (TPrimalElementDataType::*pPrimalValueFunction)(const Vector&, const Matrix&) const,
    void (TAdjointElementDataType::AdjointBaseType::ScalarDerivative::*pAdjointDerivativesFunction)(
        BoundedVector<double, 3>&, const Vector&, const Matrix&) const,
    const int BufferSize,
    const double Delta,
    const double RelativeTolerance,
    const double AbsoluteTolerance)
{
    const auto& r_primal_value_func =
        [pPrimalValueFunction](const TPrimalElementDataType& rElementData,
                               const Vector& rShapeFunctions,
                               const Matrix& rShapeFunctionDerivatives) -> double {
        return (rElementData.*pPrimalValueFunction)(rShapeFunctions, rShapeFunctionDerivatives);
    };

    const auto& r_adjoint_derivative_func =
        [pAdjointDerivativesFunction, &rDerivativeVariable](
            BoundedVector<double, 3>& rOutput,
            const TAdjointElementDataType& rElementData, const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives, const int) {
            (rElementData.GetScalarDerivativeData(rDerivativeVariable).*
             pAdjointDerivativesFunction)(rOutput, rShapeFunctions, rShapeFunctionDerivatives);
        };

    RunAdjointElementDataTest<TPrimalElementDataType, TAdjointElementDataType>(
        rModel, rAddNodalSolutionStepVariablesFunction,
        rSetVariableDataFunction, rUpdateFunction, r_primal_value_func,
        GetPerturbationMethod(rDerivativeVariable), r_adjoint_derivative_func,
        BufferSize, Delta, RelativeTolerance, AbsoluteTolerance);
}

template <class TPrimalElementDataType, class TAdjointElementDataType>
void RunAdjointElementDataTest(
    Model& rModel,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFunction,
    const std::function<void(ModelPart& rModelPart)>& rSetVariableDataFunction,
    const std::function<void(ModelPart& rModelPart)>& rUpdateFunction,
    const Variable<array_1d<double, 3>>& rDerivativeVariable,
    double (TPrimalElementDataType::*pPrimalValueFunction)(const Vector&, const Matrix&) const,
    void (TAdjointElementDataType::AdjointBaseType::VectorDerivative::*pAdjointDerivativesFunction)(
        BoundedMatrix<double, 3, 2>&, const Vector&, const Matrix&) const,
    const int BufferSize,
    const double Delta,
    const double RelativeTolerance,
    const double AbsoluteTolerance)
{
    const auto& r_primal_value_func =
        [pPrimalValueFunction](const TPrimalElementDataType& rElementData,
                               const Vector& rShapeFunctions,
                               const Matrix& rShapeFunctionDerivatives) -> double {
        return (rElementData.*pPrimalValueFunction)(rShapeFunctions, rShapeFunctionDerivatives);
    };

    const auto& r_adjoint_derivative_func =
        [pAdjointDerivativesFunction, &rDerivativeVariable](
            BoundedMatrix<double, 3, 2>& rOutput,
            const TAdjointElementDataType& rElementData, const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives, const int) {
            (rElementData.GetVectorDerivativeData(rDerivativeVariable).*
             pAdjointDerivativesFunction)(rOutput, rShapeFunctions, rShapeFunctionDerivatives);
        };

    RunAdjointElementDataTest<TPrimalElementDataType, TAdjointElementDataType>(
        rModel, rAddNodalSolutionStepVariablesFunction,
        rSetVariableDataFunction, rUpdateFunction, r_primal_value_func,
        GetPerturbationMethod(rDerivativeVariable), r_adjoint_derivative_func,
        BufferSize, Delta, RelativeTolerance, AbsoluteTolerance);
}

template <class TPrimalElementDataType, class TAdjointElementDataType>
void RunAdjointElementDataTest(
    Model& rModel,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFunction,
    const std::function<void(ModelPart& rModelPart)>& rSetVariableDataFunction,
    const std::function<void(ModelPart& rModelPart)>& rUpdateFunction,
    double (TPrimalElementDataType::*pPrimalValueFunction)(const Vector&, const Matrix&) const,
    double (TAdjointElementDataType::AdjointBaseType::ShapeDerivative::*pAdjointDerivativesFunction)(
        const ShapeParameter&,
        const Vector&,
        const Matrix&,
        const double,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType&) const,
    const int BufferSize,
    const double Delta,
    const double RelativeTolerance,
    const double AbsoluteTolerance)
{
    const auto& r_primal_value_func =
        [pPrimalValueFunction](const TPrimalElementDataType& rElementData,
                               const Vector& rShapeFunctions,
                               const Matrix& rShapeFunctionDerivatives) -> double {
        return (rElementData.*pPrimalValueFunction)(rShapeFunctions, rShapeFunctionDerivatives);
    };

    const auto& r_adjoint_derivative_func =
        [pAdjointDerivativesFunction](
            BoundedMatrix<double, 3, 2>& rOutput,
            const TAdjointElementDataType& rElementData, const Vector& rShapeFunctions,
            const Matrix& rShapeFunctionDerivatives, const int GaussPointIndex) {
            Geometry<Point>::JacobiansType J;
            rElementData.GetGeometry().Jacobian(
                J, TAdjointElementDataType::GetIntegrationMethod());
            const auto& DN_De = rElementData.GetGeometry().ShapeFunctionsLocalGradients(
                TAdjointElementDataType::GetIntegrationMethod());

            GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_deriv;
            const Matrix& rJ = J[GaussPointIndex];
            const Matrix& rDN_De = DN_De[GaussPointIndex];
            GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

            ShapeParameter deriv;
            for (int c = 0; c < 3; ++c)
            {
                for (int k = 0; k < 2; ++k)
                {
                    deriv.NodeIndex = c;
                    deriv.Direction = k;

                    double detJ_deriv;
                    geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, DN_DX_deriv);

                    rOutput(c, k) =
                        (rElementData.GetShapeDerivativeData(SHAPE_SENSITIVITY).*
                         pAdjointDerivativesFunction)(deriv, rShapeFunctions,
                                                      rShapeFunctionDerivatives,
                                                      detJ_deriv, DN_DX_deriv);
                }
            }
        };

    RunAdjointElementDataTest<TPrimalElementDataType, TAdjointElementDataType>(
        rModel, rAddNodalSolutionStepVariablesFunction,
        rSetVariableDataFunction, rUpdateFunction, r_primal_value_func,
        GetPerturbationMethod(SHAPE_SENSITIVITY), r_adjoint_derivative_func,
        BufferSize, Delta, RelativeTolerance, AbsoluteTolerance);
}

template <class TPrimalElementDataType, class TAdjointElementDataType>
void RunAdjointElementDataTest(
    Model& rModel,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFunction,
    const std::function<void(ModelPart& rModelPart)>& rSetVariableDataFunction,
    const std::function<void(ModelPart& rModelPart)>& rUpdateFunction,
    const Variable<double>& rDerivativeVariable,
    array_1d<double, 3> (TPrimalElementDataType::*pPrimalValueFunction)(const Vector&, const Matrix&)
        const,
    void (TAdjointElementDataType::AdjointBaseType::ScalarDerivative::*pAdjointDerivativesFunction)(
        BoundedMatrix<double, 3, 2>&, const Vector&, const Matrix&) const,
    const int BufferSize,
    const double Delta,
    const double RelativeTolerance,
    const double AbsoluteTolerance)
{
    for (int i_dim = 0; i_dim < 2; ++i_dim)
    {
        const auto& r_primal_value_func =
            [pPrimalValueFunction, i_dim](
                const TPrimalElementDataType& rElementData, const Vector& rShapeFunctions,
                const Matrix& rShapeFunctionDerivatives) -> double {
            return (rElementData.*pPrimalValueFunction)(
                rShapeFunctions, rShapeFunctionDerivatives)[i_dim];
        };

        const auto& r_adjoint_derivative_func =
            [pAdjointDerivativesFunction, i_dim, &rDerivativeVariable](
                BoundedVector<double, 3>& rOutput,
                const TAdjointElementDataType& rElementData, const Vector& rShapeFunctions,
                const Matrix& rShapeFunctionDerivatives, const int) {
                BoundedMatrix<double, 3, 2> output;
                (rElementData.GetScalarDerivativeData(rDerivativeVariable).*
                 pAdjointDerivativesFunction)(output, rShapeFunctions, rShapeFunctionDerivatives);
                noalias(rOutput) = column(output, i_dim);
            };

        RunAdjointElementDataTest<TPrimalElementDataType, TAdjointElementDataType>(
            rModel, rAddNodalSolutionStepVariablesFunction,
            rSetVariableDataFunction, rUpdateFunction, r_primal_value_func,
            GetPerturbationMethod(rDerivativeVariable), r_adjoint_derivative_func,
            BufferSize, Delta, RelativeTolerance, AbsoluteTolerance);
    }
}

template <class TPrimalElementDataType, class TAdjointElementDataType>
void RunAdjointElementDataTest(
    Model& rModel,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFunction,
    const std::function<void(ModelPart& rModelPart)>& rSetVariableDataFunction,
    const std::function<void(ModelPart& rModelPart)>& rUpdateFunction,
    const Variable<array_1d<double, 3>>& rDerivativeVariable,
    array_1d<double, 3> (TPrimalElementDataType::*pPrimalValueFunction)(const Vector&, const Matrix&)
        const,
    void (TAdjointElementDataType::AdjointBaseType::VectorDerivative::*pAdjointDerivativesFunction)(
        BoundedMatrix<double, 6, 2>&, const Vector&, const Matrix&) const,
    const int BufferSize,
    const double Delta,
    const double RelativeTolerance,
    const double AbsoluteTolerance)
{
    for (int i_dim = 0; i_dim < 2; ++i_dim)
    {
        const auto& r_primal_value_func =
            [pPrimalValueFunction, i_dim](
                const TPrimalElementDataType& rElementData, const Vector& rShapeFunctions,
                const Matrix& rShapeFunctionDerivatives) -> double {
            return (rElementData.*pPrimalValueFunction)(
                rShapeFunctions, rShapeFunctionDerivatives)[i_dim];
        };

        const auto& r_adjoint_derivative_func =
            [pAdjointDerivativesFunction, i_dim, &rDerivativeVariable](
                BoundedMatrix<double, 3, 2>& rOutput,
                const TAdjointElementDataType& rElementData, const Vector& rShapeFunctions,
                const Matrix& rShapeFunctionDerivatives, const int) {
                BoundedMatrix<double, 6, 2> output;
                (rElementData.GetVectorDerivativeData(rDerivativeVariable).*
                 pAdjointDerivativesFunction)(output, rShapeFunctions, rShapeFunctionDerivatives);

                for (int c = 0; c < 3; ++c)
                {
                    for (int k = 0; k < 2; ++k)
                    {
                        rOutput(c, k) = output(c * 2 + k, i_dim);
                    }
                }
            };

        RunAdjointElementDataTest<TPrimalElementDataType, TAdjointElementDataType>(
            rModel, rAddNodalSolutionStepVariablesFunction,
            rSetVariableDataFunction, rUpdateFunction, r_primal_value_func,
            GetPerturbationMethod(rDerivativeVariable), r_adjoint_derivative_func,
            BufferSize, Delta, RelativeTolerance, AbsoluteTolerance);
    }
}

template <class TPrimalElementDataType, class TAdjointElementDataType>
void RunAdjointElementDataTest(
    Model& rModel,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFunction,
    const std::function<void(ModelPart& rModelPart)>& rSetVariableDataFunction,
    const std::function<void(ModelPart& rModelPart)>& rUpdateFunction,
    array_1d<double, 3> (TPrimalElementDataType::*pPrimalValueFunction)(const Vector&, const Matrix&)
        const,
    array_1d<double, 3> (TAdjointElementDataType::AdjointBaseType::ShapeDerivative::*pAdjointDerivativesFunction)(
        const ShapeParameter&,
        const Vector&,
        const Matrix&,
        const double,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType&) const,
    const int BufferSize,
    const double Delta,
    const double RelativeTolerance,
    const double AbsoluteTolerance)
{
    for (int i_dim = 0; i_dim < 2; ++i_dim)
    {
        const auto& r_primal_value_func =
            [pPrimalValueFunction, i_dim](
                const TPrimalElementDataType& rElementData, const Vector& rShapeFunctions,
                const Matrix& rShapeFunctionDerivatives) -> double {
            return (rElementData.*pPrimalValueFunction)(
                rShapeFunctions, rShapeFunctionDerivatives)[i_dim];
        };

        const auto& r_adjoint_derivative_func =
            [pAdjointDerivativesFunction, i_dim](
                BoundedMatrix<double, 3, 2>& rOutput,
                const TAdjointElementDataType& rElementData, const Vector& rShapeFunctions,
                const Matrix& rShapeFunctionDerivatives, const int GaussPointIndex) {
                Geometry<Point>::JacobiansType J;
                rElementData.GetGeometry().Jacobian(
                    J, TAdjointElementDataType::GetIntegrationMethod());
                const auto& DN_De = rElementData.GetGeometry().ShapeFunctionsLocalGradients(
                    TAdjointElementDataType::GetIntegrationMethod());

                GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_deriv;
                const Matrix& rJ = J[GaussPointIndex];
                const Matrix& rDN_De = DN_De[GaussPointIndex];
                GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

                ShapeParameter deriv;
                for (int c = 0; c < 3; ++c)
                {
                    for (int k = 0; k < 2; ++k)
                    {
                        deriv.NodeIndex = c;
                        deriv.Direction = k;

                        double detJ_deriv;
                        geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, DN_DX_deriv);

                        const array_1d<double, 3>& r_value =
                            (rElementData.GetShapeDerivativeData(SHAPE_SENSITIVITY).*
                             pAdjointDerivativesFunction)(deriv, rShapeFunctions,
                                                          rShapeFunctionDerivatives,
                                                          detJ_deriv, DN_DX_deriv);
                        rOutput(c, k) = r_value[i_dim];
                    }
                }
            };

        RunAdjointElementDataTest<TPrimalElementDataType, TAdjointElementDataType>(
            rModel, rAddNodalSolutionStepVariablesFunction,
            rSetVariableDataFunction, rUpdateFunction, r_primal_value_func,
            GetPerturbationMethod(SHAPE_SENSITIVITY), r_adjoint_derivative_func,
            BufferSize, Delta, RelativeTolerance, AbsoluteTolerance);
    }
}

template <class TClassType>
void CalculateResidual(Vector& residual, TClassType& rClassTypeObject, const ProcessInfo& rProcessInfo);

template <typename TContainer>
TContainer& GetContainerItems(ModelPart& rModelPart);

template <typename TContainer>
void RunResidualSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    const std::vector<Process*>& rPrimalProcesses,
    const std::vector<Process*>& rAdjointProcesses,
    const std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    const std::function<void(Matrix&, typename TContainer::data_type&, const ProcessInfo&)> CalculateElementResidualScalarSensitivity,
    const std::function<double&(NodeType&)> PerturbVariable,
    const double Delta,
    const double Tolerance,
    const int DerivativesOffset = 0,
    const int EquationOffset = 0);

template <typename TContainer>
void RunResidualSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    const std::vector<Process*>& rPrimalProcesses,
    const std::vector<Process*>& rAdjointProcesses,
    const std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    const std::function<void(Matrix&, typename TContainer::data_type&, const ProcessInfo&)> CalculateElementResidualVectorSensitivity,
    const std::function<double&(NodeType&, const int)> PerturbVariable,
    const double Delta,
    const double Tolerance,
    const int DerivativesOffset = 0,
    const int EquationOffset = 0);

} // namespace RansApplicationTestUtilities
} // namespace Kratos

#endif // KRATOS_RANS_ADJOINT_TEST_UTILITIES_H_INCLUDED