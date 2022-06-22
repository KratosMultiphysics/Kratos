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

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"

// Application includes

// Include h
#include "fluid_least_squares_shadowing_utilities.h"

namespace Kratos
{

FluidLeastSquaresShadowingUtilities::FluidLeastSquaresShadowingUtilities(
    const std::vector<const Variable<double>*>& rPrimalVariablePointersList,
    const std::vector<const Variable<double>*>& rPrimalFirstDerivativeVariablePointersList,
    const std::vector<const Variable<double>*>& rAdjointVariablePointersList,
    const std::vector<const Variable<double>*>& rAdjointFirstDerivativeVariablePointersList,
    const std::vector<const Variable<double>*>& rLSSVariablePointersList,
    const std::vector<const Variable<double>*>& rLSSFirstDerivativeVariablePointersList,
    const Variable<Vector>& rAuxVariable)
    : mPrimalVariablePointersList(rPrimalVariablePointersList),
      mPrimalFirstDerivativeVariablePointersList(rPrimalFirstDerivativeVariablePointersList),
      mAdjointVariablePointersList(rAdjointVariablePointersList),
      mAdjointFirstDerivativeVariablePointersList(rAdjointFirstDerivativeVariablePointersList),
      mLSSVariablePointersList(rLSSVariablePointersList),
      mLSSFirstDerivativeVariablePointersList(rLSSFirstDerivativeVariablePointersList),
      mrAuxVariable(rAuxVariable)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mPrimalVariablePointersList.size() == mPrimalFirstDerivativeVariablePointersList.size())
        << "Primal first derivative variables list size is not matching with primal variables list size [ primal first derivative variables list size = "
        << mPrimalFirstDerivativeVariablePointersList.size() << ", primal variables list size = " << mPrimalVariablePointersList.size() << " ].\n";

    KRATOS_ERROR_IF_NOT(mPrimalVariablePointersList.size() == mAdjointVariablePointersList.size())
        << "Adjoint variables list size is not matching with primal variables list size [ adjoint variables list size = "
        << mAdjointVariablePointersList.size() << ", primal variables list size = " << mPrimalVariablePointersList.size() << " ].\n";

    KRATOS_ERROR_IF_NOT(mPrimalVariablePointersList.size() == mAdjointFirstDerivativeVariablePointersList.size())
        << "Adjoint first derivative variables list size is not matching with primal variables list size [ adjoint first derivative variables list size = "
        << mAdjointFirstDerivativeVariablePointersList.size() << ", primal variables list size = " << mPrimalVariablePointersList.size() << " ].\n";

    KRATOS_ERROR_IF_NOT(mPrimalVariablePointersList.size() == mLSSVariablePointersList.size())
        << "LSS variables list size is not matching with primal variables list size [ lss variables list size = "
        << mLSSVariablePointersList.size() << ", primal variables list size = " << mPrimalVariablePointersList.size() << " ].\n";

    KRATOS_ERROR_IF_NOT(mPrimalVariablePointersList.size() == mLSSFirstDerivativeVariablePointersList.size())
        << "LSS first derivative variables list size is not matching with primal variables list size [ lss first derivative variables list size = "
        << mLSSFirstDerivativeVariablePointersList.size() << ", primal variables list size = " << mPrimalVariablePointersList.size() << " ].\n";

    // check for common variables
    std::vector<const Variable<double>*> temp;
    const auto& variable_adder = [&](const std::vector<const Variable<double>*>& rInput) {
        for (const auto& r_variable : rInput) {
            for (const auto& r_temp_var : temp) {
                KRATOS_ERROR_IF(*r_temp_var == *r_variable) << r_variable->Name() << " is already used in this utilities. Please use unique list of variables.\n";
            }
            temp.push_back(r_variable);
        }
    };

    variable_adder(mPrimalVariablePointersList);
    variable_adder(mPrimalFirstDerivativeVariablePointersList);
    variable_adder(mAdjointVariablePointersList);
    variable_adder(mAdjointFirstDerivativeVariablePointersList);
    variable_adder(mLSSVariablePointersList);
    variable_adder(mLSSFirstDerivativeVariablePointersList);

    KRATOS_CATCH("");
}

template<class TEntityType>
void FluidLeastSquaresShadowingUtilities::CheckVariables(const TEntityType& rEntity) const
{
    KRATOS_TRY

    for (const auto& r_node : rEntity.GetGeometry()) {
        for (const auto& p_variable : mPrimalVariablePointersList) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*p_variable), r_node);
        }

        for (const auto& p_variable : mPrimalFirstDerivativeVariablePointersList) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*p_variable), r_node);
        }

        for (const auto& p_variable : mAdjointVariablePointersList) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*p_variable), r_node);
        }

        for (const auto& p_variable : mAdjointFirstDerivativeVariablePointersList) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*p_variable), r_node);
        }

        for (const auto& p_variable : mLSSVariablePointersList) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*p_variable), r_node);
        }

        for (const auto& p_variable : mLSSFirstDerivativeVariablePointersList) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*p_variable), r_node);
        }
    }


    KRATOS_CATCH("");
}

template<class TEntityType>
void FluidLeastSquaresShadowingUtilities::GetValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step,
    const std::vector<const Variable<double>*>& rVariablePointersList)
{
    KRATOS_TRY

    const auto& r_geometry = rEntity.GetGeometry();
    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType number_of_variables = rVariablePointersList.size();
    const IndexType local_size = number_of_nodes * number_of_variables;

    if (rOutput.size() != local_size) {
        rOutput.resize(local_size);
    }

    IndexType local_index = 0;
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        for (IndexType j = 0; j < number_of_variables; ++j) {
            rOutput[local_index++] = r_geometry[i].FastGetSolutionStepValue(*rVariablePointersList[j], Step);
        }
    }

    KRATOS_CATCH("");
}

template<class TEntityType>
void FluidLeastSquaresShadowingUtilities::GetPrimalValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step) const
{
    FluidLeastSquaresShadowingUtilities::GetValues(rOutput, rEntity, Step, mPrimalVariablePointersList);
}

template<class TEntityType>
void FluidLeastSquaresShadowingUtilities::GetPrimalFirstDerivativeValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step) const
{
    FluidLeastSquaresShadowingUtilities::GetValues(rOutput, rEntity, Step, mPrimalFirstDerivativeVariablePointersList);
}

template<class TEntityType>
void FluidLeastSquaresShadowingUtilities::GetAdjointValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step) const
{
    FluidLeastSquaresShadowingUtilities::GetValues(rOutput, rEntity, Step, mAdjointVariablePointersList);
}

template<class TEntityType>
void FluidLeastSquaresShadowingUtilities::GetAdjointFirstDerivativeValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step) const
{
    FluidLeastSquaresShadowingUtilities::GetValues(rOutput, rEntity, Step, mAdjointFirstDerivativeVariablePointersList);
}

template<class TEntityType>
void FluidLeastSquaresShadowingUtilities::GetLSSValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step) const
{
    FluidLeastSquaresShadowingUtilities::GetValues(rOutput, rEntity, Step, mLSSVariablePointersList);
}

template<class TEntityType>
void FluidLeastSquaresShadowingUtilities::GetLSSFirstDerivativeValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step) const
{
    FluidLeastSquaresShadowingUtilities::GetValues(rOutput, rEntity, Step, mLSSFirstDerivativeVariablePointersList);
}

const std::vector<const Variable<double>*>& FluidLeastSquaresShadowingUtilities::GetPrimalVariablePointersList() const
{
    return mPrimalVariablePointersList;
}

const std::vector<const Variable<double>*>& FluidLeastSquaresShadowingUtilities::GetPrimalFirstDerivativeVariablePointersList() const
{
    return mPrimalFirstDerivativeVariablePointersList;
}

const std::vector<const Variable<double>*>& FluidLeastSquaresShadowingUtilities::GetAdjointVariablePointersList() const
{
    return mAdjointVariablePointersList;
}

const std::vector<const Variable<double>*>& FluidLeastSquaresShadowingUtilities::GetAdjointFirstDerivativeVariablePointersList() const
{
    return mAdjointFirstDerivativeVariablePointersList;
}

const std::vector<const Variable<double>*>& FluidLeastSquaresShadowingUtilities::GetLSSVariablePointersList() const
{
    return mLSSVariablePointersList;
}

const std::vector<const Variable<double>*>& FluidLeastSquaresShadowingUtilities::GetLSSFirstDerivativeVariablePointersList() const
{
    return mLSSFirstDerivativeVariablePointersList;
}

const Variable<Vector>& FluidLeastSquaresShadowingUtilities::GetAuxVariable() const
{
    return mrAuxVariable;
}

// template instantiations
template void FluidLeastSquaresShadowingUtilities::GetValues(Vector&, const ConditionType&, const IndexType, const std::vector<const Variable<double>*>&);
template void FluidLeastSquaresShadowingUtilities::CheckVariables(const ConditionType&) const;
template void FluidLeastSquaresShadowingUtilities::GetPrimalValues(Vector&, const ConditionType&, const IndexType) const;
template void FluidLeastSquaresShadowingUtilities::GetPrimalFirstDerivativeValues(Vector&, const ConditionType&, const IndexType) const;
template void FluidLeastSquaresShadowingUtilities::GetAdjointValues(Vector&, const ConditionType&, const IndexType) const;
template void FluidLeastSquaresShadowingUtilities::GetAdjointFirstDerivativeValues(Vector&, const ConditionType&, const IndexType) const;
template void FluidLeastSquaresShadowingUtilities::GetLSSValues(Vector&, const ConditionType&, const IndexType) const;
template void FluidLeastSquaresShadowingUtilities::GetLSSFirstDerivativeValues(Vector&, const ConditionType&, const IndexType) const;

template void FluidLeastSquaresShadowingUtilities::GetValues(Vector&, const ElementType&, const IndexType, const std::vector<const Variable<double>*>&);
template void FluidLeastSquaresShadowingUtilities::CheckVariables(const ElementType&) const;
template void FluidLeastSquaresShadowingUtilities::GetPrimalValues(Vector&, const ElementType&, const IndexType) const;
template void FluidLeastSquaresShadowingUtilities::GetPrimalFirstDerivativeValues(Vector&, const ElementType&, const IndexType) const;
template void FluidLeastSquaresShadowingUtilities::GetAdjointValues(Vector&, const ElementType&, const IndexType) const;
template void FluidLeastSquaresShadowingUtilities::GetAdjointFirstDerivativeValues(Vector&, const ElementType&, const IndexType) const;
template void FluidLeastSquaresShadowingUtilities::GetLSSValues(Vector&, const ElementType&, const IndexType) const;
template void FluidLeastSquaresShadowingUtilities::GetLSSFirstDerivativeValues(Vector&, const ElementType&, const IndexType) const;

} // namespace Kratos

