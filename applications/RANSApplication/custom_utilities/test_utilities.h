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

#if !defined(KRATOS_RANS_TEST_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_TEST_UTILITIES_H_INCLUDED

// System includes
#include <functional>

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/geometrical_sensitivity_utility.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

namespace Kratos
{
namespace RansModellingApplicationTestUtilities
{
typedef ModelPart::NodeType NodeType;
typedef ModelPart::ElementType ElementType;
typedef ModelPart::ConditionType ConditionType;
typedef Geometry<NodeType> GeometryType;
typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;
typedef std::size_t IndexType;

bool IsNear(double ValueA, double ValueB, double RelTol = 1e-09, double AbsTol = 1e-12);

void CheckNear(double ValueA, double ValueB, double RelTol = 1e-09, double AbsTol = 1e-12);

void CheckNear(const Matrix& rA, const Matrix& rB, double RelTol = 1e-09, double AbsTol = 1e-12);

void CheckNear(const Vector& rA, const Vector& rB, double RelTol = 1e-09, double AbsTol = 1e-12);

template <class TClassType>
void CalculateResidual(Vector& residual, TClassType& rClassTypeObject, ProcessInfo& rProcessInfo);

void InitializeVariableWithValues(ModelPart& rModelPart,
                                  const Variable<double>& rVariable,
                                  const double Value,
                                  const std::size_t TimeStep = 0);

void InitializeVariableWithRandomValues(ModelPart& rModelPart,
                                        const Variable<double>& rVariable,
                                        const double MinValue,
                                        const double MaxValue,
                                        const std::size_t TimeSteps);

void InitializeVariableWithRandomValues(ModelPart& rModelPart,
                                        const Variable<array_1d<double, 3>>& rVariable,
                                        const double MinValue,
                                        const double MaxValue,
                                        const std::size_t TimeSteps);

template <typename TContainer>
void RunResidualSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    std::vector<Process*>& rPrimalProcesses,
    std::vector<Process*>& rAdjointProcesses,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<void(Matrix&, typename TContainer::data_type&, ProcessInfo&)> CalculateElementResidualScalarSensitivity,
    std::function<double&(NodeType&)> PerturbVariable,
    const double Delta,
    const double Tolerance,
    const int DerivativesOffset = 0,
    const int EquationOffset = 0);

template <typename TContainer>
void RunResidualSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    std::vector<Process*>& rPrimalProcesses,
    std::vector<Process*>& rAdjointProcesses,
    std::function<void(ModelPart&)> UpdateVariablesInModelPart,
    std::function<void(Matrix&, typename TContainer::data_type&, ProcessInfo&)> CalculateElementResidualVectorSensitivity,
    std::function<double&(NodeType&, const int)> PerturbVariable,
    const double Delta,
    const double Tolerance,
    const int DerivativesOffset = 0,
    const int EquationOffset = 0);

template <typename TContainer>
TContainer& GetContainerItems(ModelPart& rModelPart);

std::function<double&(NodeType&)> GetPerturbationMethod(const Variable<double>&);

std::function<double&(NodeType&, const int)> GetPerturbationMethod(
    const Variable<array_1d<double, 3>>&);

enum ElementDataMethods
{
    CalculateEffectiveKinematicViscosity,
    CalculateReactionTerm,
    CalculateSourceTerm
};

template <typename TElementDataType>
class ElementDataWrapper
{
public:
    ElementDataWrapper(TElementDataType& rElementData, const ProcessInfo& rProcessInfo)
        : mrElementData(rElementData), mrProcessInfo(rProcessInfo)
    {
    }

    int CalculateGeometryData()
    {
        RansCalculationUtilities::CalculateGeometryData(
            mrElementData.GetGeometry(), TElementDataType::GetIntegrationMethod(),
            mGaussWeights, mShapeFunctions, mShapeFunctionDerivatives);
        return mGaussWeights.size();
    }

    void CalculateGaussPointData(const int GaussPointIndex)
    {
        const Vector& r_gauss_point_shape_functions = row(mShapeFunctions, GaussPointIndex);
        const Matrix& r_gauss_point_shape_function_derivatives =
            mShapeFunctionDerivatives[GaussPointIndex];
        mrElementData.CalculateGaussPointData(r_gauss_point_shape_functions,
                                              r_gauss_point_shape_function_derivatives,
                                              mrProcessInfo);
    }

    double CalculatePrimalMethod(const ElementDataMethods& rMethod, const int GaussPointIndex) const
    {
        const Vector& r_gauss_point_shape_functions = row(mShapeFunctions, GaussPointIndex);
        const Matrix& r_gauss_point_shape_function_derivatives =
            mShapeFunctionDerivatives[GaussPointIndex];

        switch (rMethod)
        {
        case ElementDataMethods::CalculateEffectiveKinematicViscosity:
            return mrElementData.CalculateEffectiveKinematicViscosity(
                r_gauss_point_shape_functions,
                r_gauss_point_shape_function_derivatives, mrProcessInfo);
        case ElementDataMethods::CalculateReactionTerm:
            return mrElementData.CalculateReactionTerm(
                r_gauss_point_shape_functions,
                r_gauss_point_shape_function_derivatives, mrProcessInfo);
        case ElementDataMethods::CalculateSourceTerm:
            return mrElementData.CalculateSourceTerm(
                r_gauss_point_shape_functions,
                r_gauss_point_shape_function_derivatives, mrProcessInfo);
        }

        return 0.0;
    }

    template <typename TOutputType, typename TDataType>
    void CalculateAdjointMethod(TOutputType& rOutput,
                                const ElementDataMethods& rMethod,
                                const int GaussPointIndex,
                                const Variable<TDataType>& rVariable) const
    {
        const Vector& r_gauss_point_shape_functions = row(mShapeFunctions, GaussPointIndex);
        const Matrix& r_gauss_point_shape_function_derivatives =
            mShapeFunctionDerivatives[GaussPointIndex];

        switch (rMethod)
        {
        case ElementDataMethods::CalculateEffectiveKinematicViscosity:
            mrElementData.CalculateEffectiveKinematicViscosityDerivatives(
                rOutput, rVariable, r_gauss_point_shape_functions,
                r_gauss_point_shape_function_derivatives, mrProcessInfo);
            break;
        case ElementDataMethods::CalculateReactionTerm:
            mrElementData.CalculateReactionTermDerivatives(
                rOutput, rVariable, r_gauss_point_shape_functions,
                r_gauss_point_shape_function_derivatives, mrProcessInfo);
            break;
        case ElementDataMethods::CalculateSourceTerm:
            mrElementData.CalculateSourceTermDerivatives(
                rOutput, rVariable, r_gauss_point_shape_functions,
                r_gauss_point_shape_function_derivatives, mrProcessInfo);
            break;
        }
    }

    const TElementDataType& GetElementData() const
    {
        return mrElementData;
    }

    const Vector GetShapeFunctions(const int GaussPointIndex) const
    {
        return row(mShapeFunctions, GaussPointIndex);
    }

    const Matrix& GetShapeFunctionDerivatives(const int GaussPointIndex) const
    {
        return mShapeFunctionDerivatives[GaussPointIndex];
    }

    const ProcessInfo& GetProcessInfo() const
    {
        return mrProcessInfo;
    }

private:
    TElementDataType& mrElementData;
    const ProcessInfo& mrProcessInfo;
    Vector mGaussWeights;
    Matrix mShapeFunctions;
    ShapeFunctionDerivativesArrayType mShapeFunctionDerivatives;
};

template <typename ElementDataType>
std::function<double(const ElementDataWrapper<ElementDataType>&, const int)> GetPrimalMethod(
    const ElementDataMethods& rMethod)
{
    return [rMethod](const ElementDataWrapper<ElementDataType>& rElementDataWrapper,
                     const int GaussPointIndex) {
        return rElementDataWrapper.CalculatePrimalMethod(rMethod, GaussPointIndex);
    };
}

template <unsigned int TDim, unsigned int TNumNodes, typename ElementDataType>
std::function<void(BoundedVector<double, TNumNodes>&, const ElementDataWrapper<ElementDataType>&, const int)> GetAdjointMethod(
    const ElementDataMethods& rMethod, const Variable<double>& rVariable)
{
    return [rMethod, rVariable](BoundedVector<double, TNumNodes>& rOutput,
                                const ElementDataWrapper<ElementDataType>& rElementDataWrapper,
                                const int GaussPointIndex) {
        rElementDataWrapper.CalculateAdjointMethod(rOutput, rMethod,
                                                   GaussPointIndex, rVariable);
    };
}

template <unsigned int TDim, unsigned int TNumNodes, typename ElementDataType>
std::function<void(BoundedMatrix<double, TNumNodes, TDim>& rOutput, const ElementDataWrapper<ElementDataType>& rElementDataWrapper, const int GaussPointIndex)> GetAdjointMethod(
    const ElementDataMethods& rMethod, const Variable<array_1d<double, 3>>& rVariable)
{
    if (rVariable == SHAPE_SENSITIVITY)
    {
        return [rMethod, rVariable](BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                    const ElementDataWrapper<ElementDataType>& rElementDataWrapper,
                                    const int GaussPointIndex) {
            ShapeParameter deriv;
            Geometry<Point>::JacobiansType J;
            const GeometryType& r_geometry =
                rElementDataWrapper.GetElementData().GetGeometry();
            r_geometry.Jacobian(J, ElementDataType::GetIntegrationMethod());

            const Geometry<Point>::ShapeFunctionsGradientsType& DN_De =
                r_geometry.ShapeFunctionsLocalGradients(
                    ElementDataType::GetIntegrationMethod());

            GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_deriv;

            const Matrix& rJ = J[GaussPointIndex];
            const Matrix& rDN_De = DN_De[GaussPointIndex];
            GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

            std::function<double(const ShapeParameter&, const double, const GeometricalSensitivityUtility::ShapeFunctionsGradientType&)> sensitivity_method;

            switch (rMethod)
            {
            case ElementDataMethods::CalculateEffectiveKinematicViscosity:
                sensitivity_method =
                    [rElementDataWrapper, GaussPointIndex](
                        const ShapeParameter& rDeriv, const double detJ_deriv,
                        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_dx_deriv) {
                        return rElementDataWrapper.GetElementData().CalculateEffectiveKinematicViscosityShapeSensitivity(
                            rElementDataWrapper.GetShapeFunctions(GaussPointIndex),
                            rElementDataWrapper.GetShapeFunctionDerivatives(GaussPointIndex),
                            rDeriv, detJ_deriv, rDN_dx_deriv,
                            rElementDataWrapper.GetProcessInfo());
                    };
                break;
            case ElementDataMethods::CalculateReactionTerm:
                sensitivity_method =
                    [rElementDataWrapper, GaussPointIndex](
                        const ShapeParameter& rDeriv, const double detJ_deriv,
                        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_dx_deriv) {
                        return rElementDataWrapper.GetElementData().CalculateReactionTermShapeSensitivity(
                            rElementDataWrapper.GetShapeFunctions(GaussPointIndex),
                            rElementDataWrapper.GetShapeFunctionDerivatives(GaussPointIndex),
                            rDeriv, detJ_deriv, rDN_dx_deriv,
                            rElementDataWrapper.GetProcessInfo());
                    };
                break;
            case ElementDataMethods::CalculateSourceTerm:
                sensitivity_method =
                    [rElementDataWrapper, GaussPointIndex](
                        const ShapeParameter& rDeriv, const double detJ_deriv,
                        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_dx_deriv) {
                        return rElementDataWrapper.GetElementData().CalculateSourceTermShapeSensitivity(
                            rElementDataWrapper.GetShapeFunctions(GaussPointIndex),
                            rElementDataWrapper.GetShapeFunctionDerivatives(GaussPointIndex),
                            rDeriv, detJ_deriv, rDN_dx_deriv,
                            rElementDataWrapper.GetProcessInfo());
                    };
                break;
            }

            for (unsigned int c = 0; c < TNumNodes; ++c)
            {
                for (unsigned int k = 0; k < TDim; ++k)
                {
                    deriv.NodeIndex = c;
                    deriv.Direction = k;

                    double detJ_deriv;
                    geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, DN_DX_deriv);

                    rOutput(c, k) = sensitivity_method(deriv, detJ_deriv, DN_DX_deriv);
                }
            }
        };
    }

    return [rMethod, rVariable](BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                const ElementDataWrapper<ElementDataType>& rElementDataWrapper,
                                const int GaussPointIndex) {
        rElementDataWrapper.CalculateAdjointMethod(rOutput, rMethod,
                                                   GaussPointIndex, rVariable);
    };
}

template <unsigned int TDim, unsigned int TNumNodes, typename PrimalElementDataType, typename AdjointElementDataType>
void RunAdjointElementDataSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    std::vector<Process*>& rPrimalProcesses,
    std::vector<Process*>& rAdjointProcesses,
    std::function<double&(NodeType&)> PerturbVariable,
    std::function<double(const ElementDataWrapper<PrimalElementDataType>&, const int)> PrimalMethod,
    std::function<void(BoundedVector<double, TNumNodes>& rOutput, const ElementDataWrapper<AdjointElementDataType>&, const int GaussPointIndex)> AdjointMethod,
    const double Delta,
    const double Tolerance)
{
    ModelPart::ElementsContainerType& r_primal_container = rPrimalModelPart.Elements();
    ModelPart::ElementsContainerType& r_adjoint_container = rAdjointModelPart.Elements();

    std::size_t number_of_elements = r_primal_container.size();

    KRATOS_ERROR_IF(number_of_elements != r_adjoint_container.size())
        << "Number of elements mismatch.";

    rAdjointModelPart.GetProcessInfo()[DELTA_TIME] =
        -1.0 * rPrimalModelPart.GetProcessInfo()[DELTA_TIME];

    for (auto process : rPrimalProcesses)
        process->Check();
    for (auto process : rAdjointProcesses)
        process->Check();

    for (auto process : rAdjointProcesses)
        process->Execute();

    ProcessInfo& r_primal_process_info = rPrimalModelPart.GetProcessInfo();
    ProcessInfo& r_adjoint_process_info = rAdjointModelPart.GetProcessInfo();

    const int domain_size = r_primal_process_info[DOMAIN_SIZE];
    KRATOS_ERROR_IF(domain_size != r_adjoint_process_info[DOMAIN_SIZE])
        << "Domain size mismatch.";

    KRATOS_ERROR_IF(AdjointElementDataType::GetIntegrationMethod() != PrimalElementDataType::GetIntegrationMethod())
        << "Integration method mismatch.";

    Matrix adjoint_total_element_residual_sensitivity, damping_matrix, mass_matrix;

    for (std::size_t i_element = 0; i_element < number_of_elements; ++i_element)
    {
        ElementType& r_adjoint_element = *(r_adjoint_container.begin() + i_element);
        const GeometryType& r_adjoint_geometry = r_adjoint_element.GetGeometry();
        AdjointElementDataType::Check(r_adjoint_geometry, r_adjoint_process_info);
        AdjointElementDataType adjoint_element_data(r_adjoint_geometry);

        ElementDataWrapper<AdjointElementDataType> adjoint_element_data_wrapper(
            adjoint_element_data, r_adjoint_process_info);
        const int number_of_gauss_points =
            adjoint_element_data_wrapper.CalculateGeometryData();

        ElementType& r_primal_element = *(r_primal_container.begin() + i_element);
        GeometryType& r_primal_geometry = r_primal_element.GetGeometry();
        PrimalElementDataType::Check(r_primal_geometry, r_primal_process_info);
        PrimalElementDataType primal_elment_data(r_primal_geometry);

        ElementDataWrapper<PrimalElementDataType> primal_element_data_wrapper(
            primal_elment_data, r_primal_process_info);

        for (int g = 0; g < number_of_gauss_points; ++g)
        {
            adjoint_element_data_wrapper.CalculateGaussPointData(g);

            BoundedVector<double, TNumNodes> adjoint_sensitivities;
            AdjointMethod(adjoint_sensitivities, adjoint_element_data_wrapper, g);

            primal_element_data_wrapper.CalculateGeometryData();

            for (auto process : rPrimalProcesses)
                process->Execute();

            primal_element_data_wrapper.CalculateGaussPointData(g);

            const double primal_value_0 = PrimalMethod(primal_element_data_wrapper, g);

            BoundedVector<double, TNumNodes> fd_sensitivities;
            for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
            {
                NodeType& r_node = r_primal_geometry[i_node];
                PerturbVariable(r_node) += Delta;

                primal_element_data_wrapper.CalculateGeometryData();

                for (auto process : rPrimalProcesses)
                    process->Execute();

                primal_element_data_wrapper.CalculateGaussPointData(g);

                const double primal_value_1 =
                    PrimalMethod(primal_element_data_wrapper, g);

                fd_sensitivities[i_node] = (primal_value_1 - primal_value_0) / Delta;

                PerturbVariable(r_node) -= Delta;
            }

            KRATOS_CHECK_VECTOR_NEAR(adjoint_sensitivities, fd_sensitivities, Tolerance);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes, typename PrimalElementDataType, typename AdjointElementDataType>
void RunAdjointElementDataSensitivityTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    std::vector<Process*>& rPrimalProcesses,
    std::vector<Process*>& rAdjointProcesses,
    std::function<double&(NodeType&, const int)> PerturbVariable,
    std::function<double(const ElementDataWrapper<PrimalElementDataType>&, const int)> PrimalMethod,
    std::function<void(BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                       const ElementDataWrapper<AdjointElementDataType>&,
                       const int GaussPointIndex)> AdjointMethod,
    const double Delta,
    const double Tolerance)
{
    ProcessInfo& r_primal_process_info = rPrimalModelPart.GetProcessInfo();

    const int domain_size = r_primal_process_info[DOMAIN_SIZE];

    KRATOS_ERROR_IF(TDim != domain_size) << "Domain size mismatch";

    for (int i_dim = 0; i_dim < domain_size; ++i_dim)
    {
        auto calculate_sensitivities =
            [AdjointMethod, i_dim](BoundedVector<double, TNumNodes>& rOutput,
                                   const ElementDataWrapper<AdjointElementDataType>& rElementDataWrapper,
                                   const int GaussPointIndex) {
                BoundedMatrix<double, TNumNodes, TDim> sensitivities;
                AdjointMethod(sensitivities, rElementDataWrapper, GaussPointIndex);

                for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
                {
                    rOutput[i_node] = sensitivities(i_node, i_dim);
                }
            };

        std::function<double&(NodeType&)> perturb_variable =
            [PerturbVariable, i_dim](NodeType& rNode) -> double& {
            return PerturbVariable(rNode, i_dim);
        };

        RunAdjointElementDataSensitivityTest<TDim, TNumNodes, PrimalElementDataType, AdjointElementDataType>(
            rPrimalModelPart, rAdjointModelPart, rPrimalProcesses, rAdjointProcesses,
            perturb_variable, PrimalMethod, calculate_sensitivities, Delta, Tolerance);
    }
}

} // namespace RansModellingApplicationTestUtilities
} // namespace Kratos

#endif