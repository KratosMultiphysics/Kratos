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

// System includes

// Project includes
#include "includes/dof.h"
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/fluid_test_utilities.h"

// Include base h
#include "fluid_adjoint_test_utilities.h"

namespace Kratos
{

namespace FluidAdjointTestUtilitiesHelper
{

using IndexType = std::size_t;

using NodeType = ModelPart::NodeType;

template<class TDataType>
std::function<double&(NodeType&, const IndexType)> GetPerturbationMethod(
    const Variable<TDataType>& rPerturbationVariable);

template<>
std::function<double&(ModelPart::NodeType&, const IndexType)> GetPerturbationMethod(
    const Variable<double>& rPerturbationVariable)
{
    auto perturbation_method = [&rPerturbationVariable](
                                   NodeType& rNode, const IndexType Dimension) -> double& {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(Dimension != 0)
            << "Dimension should be always 0 for scalar perturbation "
               "variables. [ Perturbation variable = "
            << rPerturbationVariable.Name() << ", Dimension = " << Dimension << " ].\n";

        return rNode.FastGetSolutionStepValue(rPerturbationVariable);

        KRATOS_CATCH("");
    };

    return perturbation_method;
}

template <>
std::function<double&(ModelPart::NodeType&, const IndexType)> GetPerturbationMethod(
    const Variable<array_1d<double, 3>>& rPerturbationVariable)
{
    if (rPerturbationVariable == SHAPE_SENSITIVITY) {
        auto perturbation_method = [](NodeType& rNode, const IndexType iDim) -> double& {
            array_1d<double, 3>& r_coordinates = rNode.Coordinates();
            return r_coordinates[iDim];
        };
        return perturbation_method;
    } else {
        auto perturbation_method = [&rPerturbationVariable](
                                       NodeType& rNode, const IndexType Dimension) -> double& {
            KRATOS_TRY

            KRATOS_DEBUG_ERROR_IF(Dimension > 2)
                << "Dimension should be always less than 3 for scalar "
                   "perturbation "
                   "variables. [ Perturbation variable = "
                << rPerturbationVariable.Name() << ", Dimension = " << Dimension
                << " ].\n";

            array_1d<double, 3>& r_vector =
                rNode.FastGetSolutionStepValue(rPerturbationVariable);
            return r_vector[Dimension];

            KRATOS_CATCH("");
        };
        return perturbation_method;
    }
}

template<class TDataType>
IndexType GetVariableDimension(
    const Variable<TDataType>& rVariable,
    const ProcessInfo& rProcessInfo);

template <>
IndexType GetVariableDimension(
    const Variable<double>& rVariable,
    const ProcessInfo& rProcessInfo)
{
    return 1;
}

template <>
IndexType GetVariableDimension(
    const Variable<array_1d<double, 3>>& rVariable,
    const ProcessInfo& rProcessInfo)
{
    return rProcessInfo[DOMAIN_SIZE];
}

template<class TEntityType>
void CalculateResidual(
    Vector& residual,
    TEntityType& rEntity,
    const ProcessInfo& rProcessInfo)
{
    // QSVMS also supports BDF2 time scheme
    // TODO: check for the time scheme.
    const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];

    Vector nodal_scalar_values, current_nodal_scalar_rate_values, old_nodal_scalar_rate_values;
    Matrix lhs, damping_matrix, mass_matrix;

    rEntity.Initialize(rProcessInfo);
    rEntity.CalculateLocalSystem(lhs, residual, rProcessInfo);
    rEntity.CalculateMassMatrix(mass_matrix, rProcessInfo);
    rEntity.CalculateLocalVelocityContribution(damping_matrix, residual, rProcessInfo);

    rEntity.GetFirstDerivativesVector(nodal_scalar_values);
    rEntity.GetSecondDerivativesVector(current_nodal_scalar_rate_values);
    rEntity.GetSecondDerivativesVector(old_nodal_scalar_rate_values, 1);

    noalias(current_nodal_scalar_rate_values) =
        current_nodal_scalar_rate_values * (1 - bossak_alpha) +
        old_nodal_scalar_rate_values * bossak_alpha;

    IndexType residual_equations_size = std::max(
        {damping_matrix.size1(), mass_matrix.size1(), nodal_scalar_values.size(),
         current_nodal_scalar_rate_values.size(), old_nodal_scalar_rate_values.size()});

    if (residual.size() != 0) {
        KRATOS_ERROR_IF(residual.size() != residual_equations_size)
            << rEntity.Info() << "::CalculateRightHandSide RHS vector size doesn't match with max residual_equations_size "
            << residual_equations_size << " [ RHS_size = " << residual.size()
            << " != " << residual_equations_size << " ].\n";
    }

    if (damping_matrix.size1() != 0 && nodal_scalar_values.size() != 0) {
        KRATOS_ERROR_IF(damping_matrix.size1() != residual_equations_size)
            << rEntity.Info() << "::CalculateDampingMatrix damping matrix size1 doesn't match with max residual_equations_size "
            << residual_equations_size
            << " [ DampingMatrix_size1 = " << damping_matrix.size1()
            << " != " << residual_equations_size << " ].\n";
        KRATOS_ERROR_IF(damping_matrix.size1() != residual_equations_size)
            << rEntity.Info() << "::CalculateDampingMatrix damping matrix size2 doesn't match with max residual_equations_size "
            << residual_equations_size
            << " [ DampingMatrix_size2 = " << damping_matrix.size1()
            << " != " << residual_equations_size << " ].\n";
        KRATOS_ERROR_IF(nodal_scalar_values.size() != residual_equations_size)
            << rEntity.Info() << "::GetFirstDerivativesVector values vector size doesn't match with max residual_equations_size "
            << residual_equations_size
            << " [ Values_size = " << nodal_scalar_values.size()
            << " != " << residual_equations_size << " ].\n";
    }

    if (mass_matrix.size1() != 0 && current_nodal_scalar_rate_values.size() != 0) {
        KRATOS_ERROR_IF(mass_matrix.size1() != residual_equations_size)
            << rEntity.Info() << "::CalculateMassMatrix damping matrix size1 doesn't match with max residual_equations_size "
            << residual_equations_size
            << " [ MassMatrix_size1 = " << mass_matrix.size1()
            << " != " << residual_equations_size << " ].\n";
        KRATOS_ERROR_IF(mass_matrix.size1() != residual_equations_size)
            << rEntity.Info() << "::CalculateMassMatrix damping matrix size2 doesn't match with max residual_equations_size "
            << residual_equations_size
            << " [ MassMatrix_size2 = " << mass_matrix.size1()
            << " != " << residual_equations_size << " ].\n";
        KRATOS_ERROR_IF(current_nodal_scalar_rate_values.size() != residual_equations_size)
            << rEntity.Info() << "::GetSecondDerivativesVector values vector size doesn't match with max residual_equations_size "
            << residual_equations_size
            << " [ Values_size = " << current_nodal_scalar_rate_values.size()
            << " != " << residual_equations_size << " ].\n";
        noalias(residual) -= prod(mass_matrix, current_nodal_scalar_rate_values);
    }
}

template<class TContainerType, class TDataType>
void RunAdjointEntityDerivativesTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    const std::function<void(ModelPart&)>& rUpdateModelPart,
    const Variable<TDataType>& rVariable,
    const std::function<void(Matrix&, typename TContainerType::data_type&, const ProcessInfo&)>& rCalculateElementResidualDerivatives,
    const IndexType EquationOffset,
    const IndexType DerivativeOffset,
    const double Delta,
    const double Tolerance)
{
    KRATOS_TRY

    auto& r_primal_container = FluidTestUtilities::GetContainer<TContainerType>()(rPrimalModelPart);
    auto& r_adjoint_container = FluidTestUtilities::GetContainer<TContainerType>()(rAdjointModelPart);
    rAdjointModelPart.GetProcessInfo()[DELTA_TIME] =
        rPrimalModelPart.GetProcessInfo()[DELTA_TIME] * -1.0;

    const IndexType number_of_elements = r_primal_container.size();
    KRATOS_ERROR_IF(number_of_elements != r_adjoint_container.size())
        << "Mismatching number of items (i.e. elements/conditions) in "
           "primal "
           "and adjoint model parts.\n";

    const auto& r_primal_process_info = rPrimalModelPart.GetProcessInfo();
    const auto& r_adjoint_process_info = rAdjointModelPart.GetProcessInfo();

    KRATOS_ERROR_IF(r_primal_process_info[DOMAIN_SIZE] != r_adjoint_process_info[DOMAIN_SIZE])
        << "Domain size mismatch in primal and adjoint model parts.\n";

    Matrix adjoint_residual_derivatives, adjoint_primal_lhs;
    Vector residual_ref, residual, fd_derivatives, adjoint_primal_rhs;

    const auto& perturbation_method = GetPerturbationMethod(rVariable);
    const auto derivative_dimension = GetVariableDimension(rVariable, r_primal_process_info);

    for (IndexType i = 0; i < number_of_elements; ++i) {
        auto& r_primal_element = *(r_primal_container.begin() + i);
        auto& r_adjoint_element = *(r_adjoint_container.begin() + i);

        const IndexType number_of_nodes = r_primal_element.GetGeometry().PointsNumber();
        KRATOS_ERROR_IF(number_of_nodes != r_adjoint_element.GetGeometry().PointsNumber()) << "Number of nodes mismatch between primal and adjoint elements.\n";

        r_primal_element.Initialize(r_primal_process_info);
        r_adjoint_element.Initialize(r_adjoint_process_info);

        // calculate adjoint sensitivities
        rUpdateModelPart(rAdjointModelPart);
        rCalculateElementResidualDerivatives(
            adjoint_residual_derivatives, r_adjoint_element, r_adjoint_process_info);

        // compute the dofs vector
        std::vector<Dof<double>::Pointer> dofs;
        static_cast<const typename TContainerType::data_type&>(r_adjoint_element).GetDofList(dofs, r_adjoint_process_info);

        // get derivative node ids
        std::vector<int> derivative_node_ids;
        for (IndexType i = 0; i < dofs.size(); ++i) {
            if (std::find(derivative_node_ids.begin(), derivative_node_ids.end(), dofs[i]->Id()) == derivative_node_ids.end()) {
                derivative_node_ids.push_back(dofs[i]->Id());
            }
        }

        const IndexType number_of_derivative_nodes = derivative_node_ids.size();

        // calculate primal reference residuals
        rUpdateModelPart(rPrimalModelPart);
        CalculateResidual(residual_ref, r_primal_element, r_primal_process_info);

        r_adjoint_element.CalculateLocalSystem(
            adjoint_primal_lhs, adjoint_primal_rhs, r_adjoint_process_info);
        r_adjoint_element.CalculateLocalVelocityContribution(
            adjoint_primal_lhs, adjoint_primal_rhs, r_adjoint_process_info);

        // TODO : Remove these casts once non const versions of the followings are removed from
        // the element and condition base classes.
        static_cast<const typename TContainerType::data_type&>(r_primal_element).Check(r_primal_process_info);
        static_cast<const typename TContainerType::data_type&>(r_adjoint_element).Check(r_adjoint_process_info);

        const IndexType residual_block_size = residual_ref.size() / number_of_nodes;
        const IndexType adjoint_equation_block_size =
            adjoint_residual_derivatives.size2() / number_of_nodes;
        const IndexType adjoint_derivatives_block_size =
            adjoint_residual_derivatives.size1() / number_of_derivative_nodes;

        // check residuals from adjoint and primal the same
        for (IndexType c = 0; c < number_of_nodes; ++c) {
            for (IndexType k = 0; k < residual_block_size; ++k) {
                KRATOS_CHECK_RELATIVE_NEAR(
                    residual_ref[c * residual_block_size + k],
                    adjoint_primal_rhs[c * adjoint_equation_block_size + k + EquationOffset],
                    Tolerance);
            }
        }

        for (IndexType c = 0; c < number_of_derivative_nodes; ++c) {
            auto& r_node = rPrimalModelPart.GetNode(derivative_node_ids[c]);
            for (IndexType k = 0; k < derivative_dimension; ++k) {
                perturbation_method(r_node, k) += Delta;

                // calculate perturbed residual
                rUpdateModelPart(rPrimalModelPart);
                CalculateResidual(residual, r_primal_element, r_primal_process_info);
                fd_derivatives = (residual - residual_ref) / Delta;

                Vector analytical_derivatives(number_of_nodes * residual_block_size);

                // checking fd derivatives and adjoint derivatives
                for (IndexType a = 0; a < number_of_nodes; ++a) {
                    for (IndexType b = 0; b < residual_block_size; ++b) {
                        analytical_derivatives[a * residual_block_size + b] = adjoint_residual_derivatives(
                            c * adjoint_derivatives_block_size + DerivativeOffset + k,
                            a * adjoint_equation_block_size + EquationOffset + b);
                    }
                }

                KRATOS_CHECK_VECTOR_RELATIVE_NEAR(analytical_derivatives, fd_derivatives, Tolerance);

                perturbation_method(r_node, k) -= Delta;
            }
        }
    }

    KRATOS_CATCH("");
}

}

template<class TDataType>
TDataType FluidAdjointTestUtilities::CalculateRelaxedVariableRate(
    const double BossakAlpha,
    const Variable<TDataType>& rVariable,
    const NodeType& rNode)
{
    return (1 - BossakAlpha) * rNode.FastGetSolutionStepValue(rVariable) +
           BossakAlpha * rNode.FastGetSolutionStepValue(rVariable, 1);
}

template<class TDataType>
void FluidAdjointTestUtilities::RunAdjointEntityDerivativesTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    const std::function<void(ModelPart&)>& rUpdateModelPart,
    const Variable<TDataType>& rVariable,
    const std::function<void(Matrix&, ConditionType&, const ProcessInfo&)>& rCalculateElementResidualDerivatives,
    const IndexType EquationOffset,
    const IndexType DerivativeOffset,
    const double Delta,
    const double Tolerance)
{
    FluidAdjointTestUtilitiesHelper::RunAdjointEntityDerivativesTest<ModelPart::ConditionsContainerType>(
        rPrimalModelPart, rAdjointModelPart, rUpdateModelPart, rVariable,
        rCalculateElementResidualDerivatives, EquationOffset, DerivativeOffset,
        Delta, Tolerance);
}

template<class TDataType>
void FluidAdjointTestUtilities::RunAdjointEntityDerivativesTest(
    ModelPart& rPrimalModelPart,
    ModelPart& rAdjointModelPart,
    const std::function<void(ModelPart&)>& rUpdateModelPart,
    const Variable<TDataType>& rVariable,
    const std::function<void(Matrix&, ElementType&, const ProcessInfo&)>& rCalculateElementResidualDerivatives,
    const IndexType EquationOffset,
    const IndexType DerivativeOffset,
    const double Delta,
    const double Tolerance)
{
    FluidAdjointTestUtilitiesHelper::RunAdjointEntityDerivativesTest<ModelPart::ElementsContainerType>(
        rPrimalModelPart, rAdjointModelPart, rUpdateModelPart, rVariable,
        rCalculateElementResidualDerivatives, EquationOffset, DerivativeOffset,
        Delta, Tolerance);
}


// template instantiations
template double FluidAdjointTestUtilities::CalculateRelaxedVariableRate(const double, const Variable<double>&, const NodeType&);
template array_1d<double, 3> FluidAdjointTestUtilities::CalculateRelaxedVariableRate(const double, const Variable<array_1d<double, 3>>&, const NodeType&);

template void FluidAdjointTestUtilities::RunAdjointEntityDerivativesTest(ModelPart&, ModelPart&, const std::function<void(ModelPart&)>&, const Variable<double>&, const std::function<void(Matrix&, Condition&, const ProcessInfo&)>&, const IndexType, const IndexType, const double, const double);
template void FluidAdjointTestUtilities::RunAdjointEntityDerivativesTest(ModelPart&, ModelPart&, const std::function<void(ModelPart&)>&, const Variable<array_1d<double, 3>>&, const std::function<void(Matrix&, Condition&, const ProcessInfo&)>&, const IndexType, const IndexType, const double, const double);

template void FluidAdjointTestUtilities::RunAdjointEntityDerivativesTest(ModelPart&, ModelPart&, const std::function<void(ModelPart&)>&, const Variable<double>&, const std::function<void(Matrix&, Element&, const ProcessInfo&)>&, const IndexType, const IndexType, const double, const double);
template void FluidAdjointTestUtilities::RunAdjointEntityDerivativesTest(ModelPart&, ModelPart&, const std::function<void(ModelPart&)>&, const Variable<array_1d<double, 3>>&, const std::function<void(Matrix&, Element&, const ProcessInfo&)>&, const IndexType, const IndexType, const double, const double);

} // namespace Kratos
