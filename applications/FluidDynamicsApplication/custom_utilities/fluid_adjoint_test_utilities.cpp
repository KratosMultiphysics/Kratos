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
#include "includes/dof.h"
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
ModelPart::ElementsContainerType& FluidAdjointTestUtilities::ContainerUtilities<ModelPart::ElementsContainerType>::GetContainer(ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <>
ModelPart::ConditionsContainerType& FluidAdjointTestUtilities::ContainerUtilities<ModelPart::ConditionsContainerType>::GetContainer(ModelPart& rModelPart)
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

    IndexType eq_id = 1;
    for (auto& r_node : r_model_part.Nodes()) {
        for (auto& p_dof : r_node.GetDofs()) {
            p_dof->SetEquationId(eq_id++);
        }
    }

    auto p_elem_prop = rGetElementProperties(r_model_part);

    using nid_list = std::vector<ModelPart::IndexType>;

    auto p_element = r_model_part.CreateNewElement(rElementName, 1, nid_list{3, 2, 1}, p_elem_prop);

    auto p_cond_prop = rGetConditionProperties(r_model_part);
    auto p_condition_1 = r_model_part.CreateNewCondition(rConditionName, 1, nid_list{2, 1}, p_cond_prop);
    auto p_condition_2 = r_model_part.CreateNewCondition(rConditionName, 2, nid_list{2, 3}, p_cond_prop);
    auto p_condition_3 = r_model_part.CreateNewCondition(rConditionName, 3, nid_list{3, 1}, p_cond_prop);

    p_condition_1->SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>({p_element}));
    p_condition_2->SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>({p_element}));
    p_condition_3->SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>({p_element}));

    return r_model_part;
}

template <class TContainerType>
void FluidAdjointTestUtilities::ContainerUtilities<TContainerType>::RunAdjointEntityGetDofListTest(
    ModelPart& rModelPart,
    const std::vector<const Variable<double>*>& rDofVariablesList)
{
    KRATOS_TRY

    auto& r_container = GetContainer(rModelPart);
    const auto& r_process_info = rModelPart.GetProcessInfo();

    for (IndexType i = 0; i < r_container.size(); ++i) {
        const auto& r_entity = *(r_container.begin() + i);

        DofsVectorType dofs;
        r_entity.GetDofList(dofs, r_process_info);

        auto& r_geometry = r_entity.GetGeometry();

        KRATOS_CHECK_EQUAL(dofs.size(), r_geometry.PointsNumber() * rDofVariablesList.size());

        IndexType index = 0;
        for (const auto& r_node : r_geometry) {
            for (auto p_variable : rDofVariablesList) {
                KRATOS_CHECK_EQUAL(r_node.pGetDof(*p_variable), dofs[index++]);
            }
        }
    }

    KRATOS_CATCH("");
}

template <class TContainerType>
void FluidAdjointTestUtilities::ContainerUtilities<TContainerType>::RunAdjointEntityEquationIdVectorTest(
    ModelPart& rModelPart,
    const std::vector<const Variable<double>*>& rDofVariablesList)
{
    KRATOS_TRY

    auto& r_container = GetContainer(rModelPart);
    const auto& r_process_info = rModelPart.GetProcessInfo();

    for (IndexType i = 0; i < r_container.size(); ++i) {
        const auto& r_entity = *(r_container.begin() + i);

        EquationIdVectorType equation_ids;
        r_entity.EquationIdVector(equation_ids, r_process_info);

        auto& r_geometry = r_entity.GetGeometry();

        KRATOS_CHECK_EQUAL(equation_ids.size(),r_geometry.PointsNumber() * rDofVariablesList.size());

        IndexType index = 0;
        for (const auto& r_node : r_geometry) {
            for (auto p_variable : rDofVariablesList) {
                KRATOS_CHECK_EQUAL(r_node.GetDof(*p_variable).EquationId(), equation_ids[index++]);
            }
        }
    }

    KRATOS_CATCH("");
}

template <class TContainerType>
void FluidAdjointTestUtilities::ContainerUtilities<TContainerType>::RunAdjointEntityGetValuesVectorTest(
    ModelPart& rModelPart,
    const std::vector<const Variable<double>*>& rDofVariablesList)
{
    KRATOS_TRY

    auto& r_container = GetContainer(rModelPart);

    for (IndexType i = 0; i < r_container.size(); ++i) {
        const auto& r_entity = *(r_container.begin() + i);

        Vector values;
        r_entity.GetValuesVector(values);

        auto& r_geometry = r_entity.GetGeometry();

        KRATOS_CHECK_EQUAL(values.size(),r_geometry.PointsNumber() * rDofVariablesList.size());

        IndexType index = 0;
        for (const auto& r_node : r_geometry) {
            for (auto p_variable : rDofVariablesList) {
                KRATOS_CHECK_EQUAL(r_node.FastGetSolutionStepValue(*p_variable), values[index++]);
            }
        }
    }

    KRATOS_CATCH("");
}

template <class TContainerType>
void FluidAdjointTestUtilities::ContainerUtilities<TContainerType>::RunAdjointEntityGetFirstDerivativesVectorTest(
    ModelPart& rModelPart,
    const std::function<Vector(const ModelPart::NodeType&)>& rValueRetrievalMethod)
{
    KRATOS_TRY

    auto& r_container = GetContainer(rModelPart);

    for (IndexType i = 0; i < r_container.size(); ++i) {
        const auto& r_entity = *(r_container.begin() + i);

        Vector values;
        r_entity.GetFirstDerivativesVector(values);

        IndexType index = 0;
        for (const auto& r_node : r_entity.GetGeometry()) {
            const Vector& ref_values = rValueRetrievalMethod(r_node);
            for (IndexType i = 0; i < ref_values.size(); ++i) {
                KRATOS_CHECK_EQUAL(ref_values[i], values[index++]);
            }
        }
    }

    KRATOS_CATCH("");
}

template <class TContainerType>
void FluidAdjointTestUtilities::ContainerUtilities<TContainerType>::RunAdjointEntityGetSecondDerivativesVectorTest(
    ModelPart& rModelPart,
    const std::function<Vector(const ModelPart::NodeType&)>& rValueRetrievalMethod)
{
    KRATOS_TRY

    auto& r_container = GetContainer(rModelPart);

    for (IndexType i = 0; i < r_container.size(); ++i) {
        const auto& r_entity = *(r_container.begin() + i);

        Vector values;
        r_entity.GetSecondDerivativesVector(values);

        IndexType index = 0;
        for (const auto& r_node : r_entity.GetGeometry()) {
            const Vector& ref_values = rValueRetrievalMethod(r_node);
            for (IndexType i = 0; i < ref_values.size(); ++i) {
                KRATOS_CHECK_EQUAL(ref_values[i], values[index++]);
            }
        }
    }

    KRATOS_CATCH("");
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

template class FluidAdjointTestUtilities::ContainerUtilities<ModelPart::ElementsContainerType>;

template class FluidAdjointTestUtilities::ContainerUtilities<ModelPart::ConditionsContainerType>;

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
void FluidAdjointTestUtilities::ContainerDataTypeUtilities<TContainerType, TDataType>::RunAdjointEntityTest(
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

    auto& r_primal_container = ContainerUtilities<TContainerType>::GetContainer(rPrimalModelPart);
    auto& r_adjoint_container = ContainerUtilities<TContainerType>::GetContainer(rAdjointModelPart);
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

    const auto& perturbation_method = DataTypeUtilities<TDataType>::GetPerturbationMethod(rVariable);
    const auto derivative_dimension = DataTypeUtilities<TDataType>::GetVariableDimension(rVariable, r_primal_process_info);

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

                // checking fd derivatives and adjoint derivatives
                for (IndexType a = 0; a < number_of_nodes; ++a) {
                    for (IndexType b = 0; b < residual_block_size; ++b) {
                        const double adjoint_derivative_value = adjoint_residual_derivatives(
                            c * adjoint_derivatives_block_size + DerivativeOffset + k,
                            a * adjoint_equation_block_size + EquationOffset + b);
                        const double fd_derivative_value =
                            fd_derivatives[a * residual_block_size + b];

                        KRATOS_CHECK_RELATIVE_NEAR(
                            adjoint_derivative_value, fd_derivative_value, Tolerance);
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
    auto& container = ContainerUtilities<TContainerType>::GetContainer(rModelPart);
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