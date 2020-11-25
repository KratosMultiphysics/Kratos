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
#include <functional>
#include <random>
#include <sstream>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/model_part.h"

// Application includes

// Include base h
#include "fluid_adjoint_test_utilities.h"

namespace Kratos
{

//// Static Operations ///////////////////////////////////////////////////////////////////

template <class TClassType>
void FluidAdjointTestUtilities::CalculateResidual(
    Vector& residual,
    TClassType& rClassTypeObject,
    const ProcessInfo& rProcessInfo)
{
    const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];

    Vector nodal_scalar_values, current_nodal_scalar_rate_values, old_nodal_scalar_rate_values;
    Matrix lhs, damping_matrix, mass_matrix;

    rClassTypeObject.CalculateLocalSystem(lhs, residual, rProcessInfo);
    rClassTypeObject.CalculateMassMatrix(mass_matrix, rProcessInfo);
    rClassTypeObject.CalculateLocalVelocityContribution(damping_matrix, residual, rProcessInfo);

    // TODO : Remove these casts once non const versions of the followings are removed from
    // the element and condition base classes.
    static_cast<const TClassType&>(rClassTypeObject).GetFirstDerivativesVector(nodal_scalar_values);
    static_cast<const TClassType&>(rClassTypeObject).GetSecondDerivativesVector(current_nodal_scalar_rate_values);
    static_cast<const TClassType&>(rClassTypeObject).GetSecondDerivativesVector(old_nodal_scalar_rate_values, 1);

    noalias(current_nodal_scalar_rate_values) =
        current_nodal_scalar_rate_values * (1 - bossak_alpha) +
        old_nodal_scalar_rate_values * bossak_alpha;

    IndexType residual_equations_size = std::max(
        {damping_matrix.size1(), mass_matrix.size1(), nodal_scalar_values.size(),
         current_nodal_scalar_rate_values.size(), old_nodal_scalar_rate_values.size()});

    if (residual.size() != 0) {
        KRATOS_ERROR_IF(residual.size() != residual_equations_size)
            << rClassTypeObject.Info() << "::CalculateRightHandSide RHS vector size doesn't match with max residual_equations_size "
            << residual_equations_size << " [ RHS_size = " << residual.size()
            << " != " << residual_equations_size << " ].\n";
    }

    if (damping_matrix.size1() != 0 && nodal_scalar_values.size() != 0) {
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

    if (mass_matrix.size1() != 0 && current_nodal_scalar_rate_values.size() != 0) {
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

template <>
ModelPart::ElementsContainerType& FluidAdjointTestUtilities::GetContainer(ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <>
ModelPart::ConditionsContainerType& FluidAdjointTestUtilities::GetContainer(ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

ModelPart& FluidAdjointTestUtilities::CreateTestModelPart(
    Model& rModel,
    const std::string& rElementName,
    const std::string& rConditionName,
    const std::function<ModelPart::PropertiesType::Pointer(ModelPart&)>& rGetElementProperties,
    const std::function<ModelPart::PropertiesType::Pointer(ModelPart&)>& rGetConditionProperties,
    const std::function<void(ModelPart&)>& rAddNodalSolutionStepVariablesFuncion,
    const std::function<void(NodeType&)>& rAddDofsFunction,
    const int BufferSize)
{
    auto& r_model_part = rModel.CreateModelPart("test" + std::to_string(rModel.GetModelPartNames().size()), BufferSize);
    rAddNodalSolutionStepVariablesFuncion(r_model_part);

    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 0.0, 1.0, 0.0);
    r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);

    for (auto& r_node : r_model_part.Nodes()) {
        rAddDofsFunction(r_node);
    }

    auto p_elem_prop = rGetElementProperties(r_model_part);

    using nid_list = std::vector<ModelPart::IndexType>;

    r_model_part.CreateNewElement(rElementName, 1, nid_list{3, 2, 1}, p_elem_prop);

    auto p_cond_prop = rGetConditionProperties(r_model_part);
    r_model_part.CreateNewCondition(rConditionName, 1, nid_list{1, 2}, p_cond_prop);
    r_model_part.CreateNewCondition(rConditionName, 2, nid_list{2, 3}, p_cond_prop);
    r_model_part.CreateNewCondition(rConditionName, 3, nid_list{3, 1}, p_cond_prop);

    return r_model_part;
}

//// Static Operations Template Instantiations ///////////////////////////////////////////

template void FluidAdjointTestUtilities::CalculateResidual<ModelPart::ElementType>(
    Vector&,
    ModelPart::ElementType&,
    const ProcessInfo&);

template void FluidAdjointTestUtilities::CalculateResidual<ModelPart::ConditionType>(
    Vector&,
    ModelPart::ConditionType&,
    const ProcessInfo&);

//// DataTypeUtilities Static Operations /////////////////////////////////////////////////

template <class TDataType>
TDataType FluidAdjointTestUtilities::DataTypeUtilities<TDataType>::ComputeRelaxedVariableRate(
    const double BossakAlpha,
    const Variable<TDataType>& rVariable,
    const NodeType& rNode)
{
    return (1 - BossakAlpha) * rNode.FastGetSolutionStepValue(rVariable) +
           BossakAlpha * rNode.FastGetSolutionStepValue(rVariable, 1);
}

template <>
std::function<double&(ModelPart::NodeType&, const IndexType)> FluidAdjointTestUtilities::DataTypeUtilities<double>::GetPerturbationMethod(
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
std::function<double&(ModelPart::NodeType&, const IndexType)> FluidAdjointTestUtilities::DataTypeUtilities<
    array_1d<double, 3>>::GetPerturbationMethod(const Variable<array_1d<double, 3>>& rPerturbationVariable)
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

template <>
IndexType FluidAdjointTestUtilities::DataTypeUtilities<double>::GetVariableDimension(
    const Variable<double>& rVariable,
    const ProcessInfo& rProcessInfo)
{
    return 1;
}

template <>
IndexType FluidAdjointTestUtilities::DataTypeUtilities<array_1d<double, 3>>::GetVariableDimension(
    const Variable<array_1d<double, 3>>& rVariable,
    const ProcessInfo& rProcessInfo)
{
    return rProcessInfo[DOMAIN_SIZE];
}

template <>
void FluidAdjointTestUtilities::DataTypeUtilities<double>::AssignRandomValues(
    double& rValue,
    const std::string& rSeed,
    const IndexType DomainSize,
    const double MinValue,
    const double MaxValue)
{
    std::seed_seq seed(rSeed.begin(), rSeed.end());
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(MinValue, MaxValue);

    rValue = distribution(generator);
}

template <>
void FluidAdjointTestUtilities::DataTypeUtilities<array_1d<double, 3>>::AssignRandomValues(
    array_1d<double, 3>& rValue,
    const std::string& rSeed,
    const IndexType DomainSize,
    const double MinValue,
    const double MaxValue)
{
    std::seed_seq seed(rSeed.begin(), rSeed.end());
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(MinValue, MaxValue);

    rValue[0] = distribution(generator);
    rValue[1] = distribution(generator);

    if (DomainSize == 3) {
        rValue[2] = distribution(generator);
    } else {
        rValue[2] = 0.0;
    }
}

//// DataTypeUtilities Static Operations Template Instantiations //////////////////////////

template class FluidAdjointTestUtilities::DataTypeUtilities<double>;
template class FluidAdjointTestUtilities::DataTypeUtilities<array_1d<double, 3>>;

//// ContainerDataTypeUtilities Static Operations /////////////////////////////////////////

template <class TContainerType, class TDataType>
void FluidAdjointTestUtilities::ContainerDataTypeUtilities<TContainerType, TDataType>::RunAdjointElementTest(
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

    auto& r_primal_container = GetContainer<TContainerType>(rPrimalModelPart);
    auto& r_adjoint_container = GetContainer<TContainerType>(rAdjointModelPart);
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

    Matrix adjoint_residual_derivatives;
    Vector residual_ref, residual, fd_derivatives;

    const auto& perturbation_method = DataTypeUtilities<TDataType>::GetPerturbationMethod(rVariable);
    const auto derivative_dimension = DataTypeUtilities<TDataType>::GetVariableDimension(rVariable, r_primal_process_info);

    for (IndexType i = 0; i < number_of_elements; ++i) {
        auto& r_primal_element = *(r_primal_container.begin() + i);
        auto& r_adjoint_element = *(r_adjoint_container.begin() + i);

        const IndexType number_of_nodes = r_primal_element.GetGeometry().PointsNumber();
        KRATOS_ERROR_IF(number_of_nodes != r_adjoint_element.GetGeometry().PointsNumber()) << "Number of nodes mismatch between primal and adjoint elements.\n";

        r_primal_element.Initialize(r_primal_process_info);
        r_adjoint_element.Initialize(r_adjoint_process_info);

        // TODO : Remove these casts once non const versions of the followings are removed from
        // the element and condition base classes.
        static_cast<const typename TContainerType::data_type&>(r_primal_element).Check(r_primal_process_info);
        static_cast<const typename TContainerType::data_type&>(r_adjoint_element).Check(r_adjoint_process_info);

        // calculate adjoint sensitivities
        rUpdateModelPart(rAdjointModelPart);
        rCalculateElementResidualDerivatives(
            adjoint_residual_derivatives, r_adjoint_element, r_adjoint_process_info);

        // calculate primal reference residuals
        rUpdateModelPart(rPrimalModelPart);
        CalculateResidual(residual_ref, r_primal_element, r_primal_process_info);

        const IndexType residual_block_size = residual_ref.size() / number_of_nodes;
        const IndexType adjoint_equation_block_size =
            adjoint_residual_derivatives.size2() / number_of_nodes;
        const IndexType adjoint_derivatives_block_size =
            adjoint_residual_derivatives.size1() / number_of_nodes;

        for (IndexType c = 0; c < number_of_nodes; ++c) {
            auto& r_node = r_primal_element.GetGeometry()[c];
            for (IndexType k = 0; k < derivative_dimension; ++k) {
                perturbation_method(r_node, k) += Delta;

                // calculate perturbed residual
                rUpdateModelPart(rPrimalModelPart);
                CalculateResidual(residual, r_primal_element, r_primal_process_info);
                fd_derivatives = (residual - residual_ref) / Delta;

                // checking fd derivatives and adjoint derivatives
                for (IndexType a = 0; a < number_of_nodes; ++a) {
                    for (IndexType b = 0; b < residual_block_size; ++b) {
                        const double adjoint_derivative_value = adjoint_residual_derivatives(
                            c * adjoint_derivatives_block_size + DerivativeOffset + k,
                            a * adjoint_equation_block_size + EquationOffset + b);
                        const double fd_derivative_value =
                            fd_derivatives[a * residual_block_size + b];

                        KRATOS_CHECK_RELATIVE_NEAR(
                            fd_derivative_value, adjoint_derivative_value, Tolerance);
                    }
                }

                perturbation_method(r_node, k) -= Delta;
            }
        }
    }

    KRATOS_CATCH("");
}

template <class TContainerType, class TDataType>
void FluidAdjointTestUtilities::ContainerDataTypeUtilities<TContainerType, TDataType>::RandomFillNodalHistoricalVariable(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const double MinValue,
    const double MaxValue,
    const int Step)
{
    for (auto& node : rModelPart.Nodes()) {
        std::stringstream seed;
        seed << node.Id() << "_HistoricalV_" << rVariable.Name();
        TDataType& value = node.FastGetSolutionStepValue(rVariable, Step);
        DataTypeUtilities<TDataType>::AssignRandomValues(value, seed.str(), rModelPart.GetProcessInfo()[DOMAIN_SIZE], MinValue, MaxValue);
    }
}

template <class TContainerType, class TDataType>
void FluidAdjointTestUtilities::ContainerDataTypeUtilities<TContainerType, TDataType>::RandomFillContainerVariable(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const double MinValue,
    const double MaxValue)
{
    auto& container = GetContainer<TContainerType>(rModelPart);
    for (auto& item : container) {
        std::stringstream seed;
        seed << item.Id() << "_NonHistoricalV_" << rVariable.Name();
        TDataType value = rVariable.Zero();
        DataTypeUtilities<TDataType>::AssignRandomValues(value, seed.str(), rModelPart.GetProcessInfo()[DOMAIN_SIZE], MinValue, MaxValue);
        item.SetValue(rVariable, value);
    }
}

//// ContainerDataTypeUtilities Static Operations Template Instantiations /////////////////

template class FluidAdjointTestUtilities::ContainerDataTypeUtilities<ModelPart::ElementsContainerType, double>;
template class FluidAdjointTestUtilities::ContainerDataTypeUtilities<ModelPart::ElementsContainerType, array_1d<double, 3>>;

template class FluidAdjointTestUtilities::ContainerDataTypeUtilities<ModelPart::ConditionsContainerType, double>;
template class FluidAdjointTestUtilities::ContainerDataTypeUtilities<ModelPart::ConditionsContainerType, array_1d<double, 3>>;

} // namespace Kratos