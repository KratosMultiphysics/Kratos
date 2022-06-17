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
    const std::vector<const Variable<double>*>& rLSSFirstDerivativeVariablePointersList)
    : mPrimalVariablePointersList(rPrimalVariablePointersList),
      mPrimalFirstDerivativeVariablePointersList(rPrimalFirstDerivativeVariablePointersList),
      mAdjointVariablePointersList(rAdjointVariablePointersList),
      mAdjointFirstDerivativeVariablePointersList(rAdjointFirstDerivativeVariablePointersList),
      mLSSVariablePointersList(rLSSVariablePointersList),
      mLSSFirstDerivativeVariablePointersList(rLSSFirstDerivativeVariablePointersList)
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

    KRATOS_CATCH("");
}

void FluidLeastSquaresShadowingUtilities::CheckVariables(const ElementType& rElement) const
{
    KRATOS_TRY

    for (const auto& r_node : rElement.GetGeometry()) {
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

void FluidLeastSquaresShadowingUtilities::GetValues(
    Vector& rOutput,
    const ElementType& rElement,
    const IndexType Step,
    const std::vector<const Variable<double>*>& rVariablePointersList)
{
    KRATOS_TRY

    const auto& r_geometry = rElement.GetGeometry();
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

void FluidLeastSquaresShadowingUtilities::GetPrimalValues(
    Vector& rOutput,
    const ElementType& rElement,
    const IndexType Step) const
{
    FluidLeastSquaresShadowingUtilities::GetValues(rOutput, rElement, Step, mPrimalVariablePointersList);
}

void FluidLeastSquaresShadowingUtilities::GetPrimalFirstDerivativeValues(
    Vector& rOutput,
    const ElementType& rElement,
    const IndexType Step) const
{
    FluidLeastSquaresShadowingUtilities::GetValues(rOutput, rElement, Step, mPrimalFirstDerivativeVariablePointersList);
}

void FluidLeastSquaresShadowingUtilities::GetAdjointValues(
    Vector& rOutput,
    const ElementType& rElement,
    const IndexType Step) const
{
    FluidLeastSquaresShadowingUtilities::GetValues(rOutput, rElement, Step, mAdjointVariablePointersList);
}

void FluidLeastSquaresShadowingUtilities::GetAdjointFirstDerivativeValues(
    Vector& rOutput,
    const ElementType& rElement,
    const IndexType Step) const
{
    FluidLeastSquaresShadowingUtilities::GetValues(rOutput, rElement, Step, mAdjointFirstDerivativeVariablePointersList);
}

void FluidLeastSquaresShadowingUtilities::GetLSSValues(
    Vector& rOutput,
    const ElementType& rElement,
    const IndexType Step) const
{
    FluidLeastSquaresShadowingUtilities::GetValues(rOutput, rElement, Step, mLSSVariablePointersList);
}

void FluidLeastSquaresShadowingUtilities::GetLSSFirstDerivativeValues(
    Vector& rOutput,
    const ElementType& rElement,
    const IndexType Step) const
{
    FluidLeastSquaresShadowingUtilities::GetValues(rOutput, rElement, Step, mLSSFirstDerivativeVariablePointersList);
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

} // namespace Kratos

