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

// System includes
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <random>
#include <sstream>

// External includes

// Project includes
#include "containers/model.h"
#include "custom_utilities/rans_variable_utilities.h"
#include "includes/checks.h"

#include "custom_utilities/rans_calculation_utilities.h"
#include "test_utilities.h"

namespace Kratos
{
namespace RansModellingApplicationTestUtilities
{
typedef ModelPart::NodeType NodeType;
typedef ModelPart::ElementType ElementType;
typedef Geometry<NodeType> GeometryType;
typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

/**
 * @brief based on python's math.IsNear()
 */
bool IsNear(const double ValueA, const double ValueB, const double RelTol, const double AbsTol)
{
    if (std::isnan(ValueA) || std::isnan(ValueB))
    {
        return false;
    }
    // Special cases inf, -inf or exact:
    if (ValueA == ValueB)
    {
        return true;
    }
    // Regular floating point number:
    return std::abs(ValueA - ValueB) <=
           std::max(RelTol * std::max(std::abs(ValueA), std::abs(ValueB)), AbsTol);
}

void CheckNear(const double ValueA, const double ValueB, const double RelTol, const double AbsTol)
{
    if (IsNear(ValueA, ValueB, RelTol, AbsTol))
    {
        if (IsNear(ValueA, 0.0, 0.0, AbsTol) &&
            IsNear(ValueB, 0.0, 0.0, AbsTol) && ValueA != 0.0 && ValueB != 0.0)
        {
            // Warn if ValueA and ValueB are non-zero but below the absolute tolerance threshold.
            KRATOS_WARNING("CheckNear")
                << "Comparing values smaller than Tolerance. ValueA / ValueB < "
                   "Tolerance [ "
                << ValueA << " / " << ValueB << " < " << AbsTol << " ]\n";
        }
    }
    else
    {
        // Currently KRATOS_ERROR doesn't handle I/O formatting so stringstream
        // is used here to create the message.
        std::stringstream msg;
        msg << std::fixed << std::scientific << std::setprecision(20);
        msg << "Check failed because\n"
            << "\t ValueA = " << ValueA << " is not close to\n"
            << "\t ValueB = " << ValueB << std::setprecision(2) << '\n'
            << "\t with rel. tol. = " << RelTol << " and abs. tol. = " << AbsTol;
        KRATOS_ERROR << msg.str();
    }
}

void CheckNear(const Matrix& rA, const Matrix& rB, const double RelTol, const double AbsTol)
{
    KRATOS_ERROR_IF(rA.size1() != rB.size1())
        << "rA.size1() = " << rA.size1()
        << " is not equal to rB.size1() = " << rB.size1();

    KRATOS_ERROR_IF(rA.size2() != rB.size2())
        << "rA.size2() = " << rA.size2()
        << " is not equal to rB.size2() = " << rB.size2();

    for (std::size_t i = 0; i < rA.size1(); ++i)
        for (std::size_t j = 0; j < rA.size2(); ++j)
            CheckNear(rA(i, j), rB(i, j), RelTol, AbsTol);
}

void CheckNear(const Vector& rA, const Vector& rB, const double RelTol, const double AbsTol)
{
    KRATOS_ERROR_IF(rA.size() != rB.size())
        << "rA.size() = " << rA.size() << " is not equal to rB.size() = " << rB.size();

    for (std::size_t i = 0; i < rA.size(); ++i)
        CheckNear(rA(i), rB(i), RelTol, AbsTol);
}

template <class TClassType>
void CalculateResidual(Vector& residual, TClassType& rClassTypeObject, ProcessInfo& rProcessInfo)
{
    const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];

    Vector nodal_scalar_values, current_nodal_scalar_rate_values, old_nodal_scalar_rate_values;
    Matrix lhs, damping_matrix, mass_matrix;

    rClassTypeObject.CalculateLocalSystem(lhs, residual, rProcessInfo);
    rClassTypeObject.CalculateMassMatrix(mass_matrix, rProcessInfo);
    rClassTypeObject.CalculateLocalVelocityContribution(damping_matrix, residual, rProcessInfo);

    rClassTypeObject.GetFirstDerivativesVector(nodal_scalar_values);
    rClassTypeObject.GetSecondDerivativesVector(current_nodal_scalar_rate_values);
    rClassTypeObject.GetSecondDerivativesVector(old_nodal_scalar_rate_values, 1);

    noalias(current_nodal_scalar_rate_values) =
        current_nodal_scalar_rate_values * (1 - bossak_alpha) +
        old_nodal_scalar_rate_values * bossak_alpha;

    IndexType residual_equations_size = std::max(
        {damping_matrix.size1(), mass_matrix.size1(), nodal_scalar_values.size(),
         current_nodal_scalar_rate_values.size(), old_nodal_scalar_rate_values.size()});

    if (residual.size() != 0)
    {
        KRATOS_ERROR_IF(residual.size() != residual_equations_size)
            << rClassTypeObject.Info() << "::CalculateRightHandSide RHS vector size doesn't match with max residual_equations_size "
            << residual_equations_size << " [ RHS_size = " << residual.size()
            << " != " << residual_equations_size << " ].\n";
    }

    if (damping_matrix.size1() != 0 && nodal_scalar_values.size() != 0)
    {
        KRATOS_ERROR_IF(damping_matrix.size1() != residual_equations_size)
            << rClassTypeObject.Info() << "::CalculateDampingMatrix damping matrix size1 doesn't match with max residual_equations_size "
            << residual_equations_size
            << " [ DampingMatrix_size1 = " << damping_matrix.size1()
            << " != " << residual_equations_size << " ].\n";
        KRATOS_ERROR_IF(damping_matrix.size1() != residual_equations_size)
            << rClassTypeObject.Info() << "::CalculateDampingMatrix damping matrix size2 doesn't match with max residual_equations_size "
            << residual_equations_size
            << " [ DampingMatrix_size2 = " << damping_matrix.size1()
            << " != " << residual_equations_size << " ].\n";
        KRATOS_ERROR_IF(nodal_scalar_values.size() != residual_equations_size)
            << rClassTypeObject.Info() << "::GetFirstDerivativesVector values vector size doesn't match with max residual_equations_size "
            << residual_equations_size
            << " [ Values_size = " << nodal_scalar_values.size()
            << " != " << residual_equations_size << " ].\n";
    }

    if (mass_matrix.size1() != 0 && current_nodal_scalar_rate_values.size() != 0)
    {
        KRATOS_ERROR_IF(mass_matrix.size1() != residual_equations_size)
            << rClassTypeObject.Info() << "::CalculateMassMatrix damping matrix size1 doesn't match with max residual_equations_size "
            << residual_equations_size
            << " [ MassMatrix_size1 = " << mass_matrix.size1()
            << " != " << residual_equations_size << " ].\n";
        KRATOS_ERROR_IF(mass_matrix.size1() != residual_equations_size)
            << rClassTypeObject.Info() << "::CalculateMassMatrix damping matrix size2 doesn't match with max residual_equations_size "
            << residual_equations_size
            << " [ MassMatrix_size2 = " << mass_matrix.size1()
            << " != " << residual_equations_size << " ].\n";
        KRATOS_ERROR_IF(current_nodal_scalar_rate_values.size() != residual_equations_size)
            << rClassTypeObject.Info() << "::GetSecondDerivativesVector values vector size doesn't match with max residual_equations_size "
            << residual_equations_size
            << " [ Values_size = " << current_nodal_scalar_rate_values.size()
            << " != " << residual_equations_size << " ].\n";
        noalias(residual) -= prod(mass_matrix, current_nodal_scalar_rate_values);
    }
}

void InitializeVariableWithValues(ModelPart& rModelPart,
                                  const Variable<double>& rVariable,
                                  const double Value,
                                  const std::size_t TimeStep)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
        r_node.FastGetSolutionStepValue(rVariable, TimeStep) = Value;
    }
}

void InitializeVariableWithRandomValues(ModelPart& rModelPart,
                                        const Variable<double>& rVariable,
                                        const double MinValue,
                                        const double MaxValue,
                                        const std::size_t TimeSteps)
{
    std::string seed_str = "KratosRANSModellingTestSeed_" + rVariable.Name();
    std::seed_seq seed(seed_str.begin(), seed_str.end());
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(MinValue, MaxValue);

    for (std::size_t i = 0; i < rModelPart.NumberOfNodes(); ++i)
    {
        NodeType& r_node = *(rModelPart.NodesBegin() + i);
        for (std::size_t i_step = 0; i_step < TimeSteps; ++i_step)
            r_node.FastGetSolutionStepValue(rVariable, i_step) = distribution(generator);
    }
}

void InitializeVariableWithRandomValues(ModelPart& rModelPart,
                                        const Variable<array_1d<double, 3>>& rVariable,
                                        const double MinValue,
                                        const double MaxValue,
                                        const std::size_t TimeSteps)
{
    std::string seed_str = "KratosRANSModellingTestSeed_" + rVariable.Name();
    std::seed_seq seed(seed_str.begin(), seed_str.end());
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(MinValue, MaxValue);

    for (std::size_t i = 0; i < rModelPart.NumberOfNodes(); ++i)
    {
        NodeType& r_node = *(rModelPart.NodesBegin() + i);
        for (std::size_t i_step = 0; i_step < TimeSteps; ++i_step)
        {
            array_1d<double, 3>& r_vector =
                r_node.FastGetSolutionStepValue(rVariable, i_step);
            r_vector[0] = distribution(generator);
            r_vector[1] = distribution(generator);
            r_vector[2] = distribution(generator);
        }
    }
}

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
    const int DerivativesOffset,
    const int EquationOffset)
{
    using TElementType = typename TContainer::data_type;

    TContainer& r_primal_container = GetContainerItems<TContainer>(rPrimalModelPart);
    TContainer& r_adjoint_container = GetContainerItems<TContainer>(rAdjointModelPart);

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
    UpdateVariablesInModelPart(rAdjointModelPart);

    ProcessInfo& r_primal_process_info = rPrimalModelPart.GetProcessInfo();
    ProcessInfo& r_adjoint_process_info = rAdjointModelPart.GetProcessInfo();

    const int domain_size = r_primal_process_info[DOMAIN_SIZE];
    KRATOS_ERROR_IF(domain_size != r_adjoint_process_info[DOMAIN_SIZE])
        << "Domain size mismatch.";

    Matrix adjoint_total_element_residual_sensitivity, damping_matrix, mass_matrix;

    for (std::size_t i_element = 0; i_element < number_of_elements; ++i_element)
    {
        TElementType& r_adjoint_element = *(r_adjoint_container.begin() + i_element);
        r_adjoint_element.Check(r_adjoint_process_info);

        CalculateElementResidualScalarSensitivity(
            adjoint_total_element_residual_sensitivity, r_adjoint_element, r_adjoint_process_info);

        TElementType& r_primal_element = *(r_primal_container.begin() + i_element);
        GeometryType& r_primal_geometry = r_primal_element.GetGeometry();
        r_primal_element.Check(r_primal_process_info);

        Vector residual, residual_0, residual_sensitivity;

        for (auto process : rPrimalProcesses)
            process->Execute();
        UpdateVariablesInModelPart(rPrimalModelPart);
        CalculateResidual(residual_0, r_primal_element, r_primal_process_info);

        const std::size_t number_of_nodes = r_primal_geometry.PointsNumber();
        const std::size_t number_of_equations = residual_0.size();
        const std::size_t residual_equation_size = number_of_equations / number_of_nodes;
        const int local_derivative_size =
            adjoint_total_element_residual_sensitivity.size1() / number_of_nodes;
        const int local_equation_size =
            adjoint_total_element_residual_sensitivity.size2() / number_of_nodes;

        residual.resize(number_of_equations);
        residual_sensitivity.resize(number_of_equations);

        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = r_primal_geometry[i_node];
            PerturbVariable(r_node) += Delta;

            for (auto process : rPrimalProcesses)
                process->Execute();

            UpdateVariablesInModelPart(rPrimalModelPart);
            CalculateResidual(residual, r_primal_element, r_primal_process_info);
            noalias(residual_sensitivity) = (residual - residual_0) / Delta;

            for (std::size_t i_check_equation = 0;
                 i_check_equation < number_of_equations; ++i_check_equation)
            {
                const std::size_t i_check_eq_node = i_check_equation / residual_equation_size;
                const std::size_t i_check_eq_dim = i_check_equation % residual_equation_size;

                const double current_adjoint_shape_sensitivity =
                    adjoint_total_element_residual_sensitivity(
                        i_node * local_derivative_size + DerivativesOffset,
                        i_check_eq_node * local_equation_size + EquationOffset + i_check_eq_dim);

                CheckNear(residual_sensitivity[i_check_equation],
                          current_adjoint_shape_sensitivity, Tolerance, 1e-12);
            }

            PerturbVariable(r_node) -= Delta;
        }
    }
}

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
    const int DerivativesOffset,
    const int EquationOffset)
{
    using TElementType = typename TContainer::data_type;

    ProcessInfo& r_primal_process_info = rPrimalModelPart.GetProcessInfo();

    const int domain_size = r_primal_process_info[DOMAIN_SIZE];

    for (int i_dim = 0; i_dim < domain_size; ++i_dim)
    {
        std::function<void(Matrix&, TElementType&, ProcessInfo&)> calculate_sensitivities =
            [CalculateElementResidualVectorSensitivity, i_dim, domain_size,
             DerivativesOffset](Matrix& rDimSensitivities, TElementType& rElement,
                                ProcessInfo& rCurrentProcessInfo) {
                Matrix sensitivities;
                CalculateElementResidualVectorSensitivity(
                    sensitivities, rElement, rCurrentProcessInfo);

                const int number_of_equations = sensitivities.size2();
                const int number_of_nodes = rElement.GetGeometry().PointsNumber();
                const int local_size = sensitivities.size1() / number_of_nodes;
                rDimSensitivities.resize(number_of_nodes, number_of_equations);

                for (int i = 0; i < number_of_nodes; ++i)
                    for (int j = 0; j < number_of_equations; ++j)
                        rDimSensitivities(i, j) = sensitivities(
                            i * local_size + i_dim + DerivativesOffset, j);
            };

        std::function<double&(NodeType&)> perturb_variable =
            [PerturbVariable, i_dim](NodeType& rNode) -> double& {
            return PerturbVariable(rNode, i_dim);
        };

        RunResidualSensitivityTest<TContainer>(
            rPrimalModelPart, rAdjointModelPart, rPrimalProcesses,
            rAdjointProcesses, UpdateVariablesInModelPart, calculate_sensitivities,
            perturb_variable, Delta, Tolerance, DerivativesOffset, EquationOffset);
    }
}

template <>
ModelPart::ElementsContainerType& GetContainerItems<ModelPart::ElementsContainerType>(ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <>
ModelPart::ConditionsContainerType& GetContainerItems<ModelPart::ConditionsContainerType>(ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

std::function<double&(NodeType&)> GetPerturbationMethod(const Variable<double>& rPerturbationVariable)
{
    std::function<double&(NodeType&)> perturbation_method =
        [rPerturbationVariable](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(rPerturbationVariable);
    };

    return perturbation_method;
}

std::function<double&(NodeType&, const int)> GetPerturbationMethod(
    const Variable<array_1d<double, 3>>& rPerturbationVariable)
{
    if (rPerturbationVariable == SHAPE_SENSITIVITY)
    {
        std::function<double&(NodeType&, const int)> perturbation_method =
            [rPerturbationVariable](NodeType& rNode, const int iDim) -> double& {
            array_1d<double, 3>& r_coordinates = rNode.Coordinates();
            return r_coordinates[iDim];
        };
        return perturbation_method;
    }
    else
    {
        std::function<double&(NodeType&, const int)> perturbation_method =
            [rPerturbationVariable](NodeType& rNode, const int iDim) -> double& {
            array_1d<double, 3>& r_vector =
                rNode.FastGetSolutionStepValue(rPerturbationVariable);
            return r_vector[iDim];
        };
        return perturbation_method;
    }
}

// templated method instantiation
template void CalculateResidual(Vector&, ElementType&, ProcessInfo&);
template void CalculateResidual(Vector&, ConditionType&, ProcessInfo&);

template void RunResidualSensitivityTest<ModelPart::ElementsContainerType>(
    ModelPart&,
    ModelPart&,
    std::vector<Process*>&,
    std::vector<Process*>&,
    std::function<void(ModelPart&)>,
    std::function<void(Matrix&, ElementType&, ProcessInfo&)>,
    std::function<double&(NodeType&)>,
    const double,
    const double,
    const int,
    const int);

template void RunResidualSensitivityTest<ModelPart::ConditionsContainerType>(
    ModelPart&,
    ModelPart&,
    std::vector<Process*>&,
    std::vector<Process*>&,
    std::function<void(ModelPart&)>,
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)>,
    std::function<double&(NodeType&)>,
    const double,
    const double,
    const int,
    const int);

template void RunResidualSensitivityTest<ModelPart::ElementsContainerType>(
    ModelPart&,
    ModelPart&,
    std::vector<Process*>&,
    std::vector<Process*>&,
    std::function<void(ModelPart&)>,
    std::function<void(Matrix&, ElementType&, ProcessInfo&)>,
    std::function<double&(NodeType&, const int)>,
    const double,
    const double,
    const int,
    const int);

template void RunResidualSensitivityTest<ModelPart::ConditionsContainerType>(
    ModelPart&,
    ModelPart&,
    std::vector<Process*>&,
    std::vector<Process*>&,
    std::function<void(ModelPart&)>,
    std::function<void(Matrix&, ConditionType&, ProcessInfo&)>,
    std::function<double&(NodeType&, const int)>,
    const double,
    const double,
    const int,
    const int);

} // namespace RansModellingApplicationTestUtilities
} // namespace Kratos