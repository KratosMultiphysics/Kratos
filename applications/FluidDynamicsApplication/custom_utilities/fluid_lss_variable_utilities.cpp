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
#include "fluid_lss_variable_utilities.h"

namespace Kratos
{

template<class TEntityType>
void FluidLSSVariableUtilities::GetValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step,
    const std::vector<IndirectVariableType>& rIndirectVariablesList)
{
    KRATOS_TRY

    const auto& r_geometry = rEntity.GetGeometry();
    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType number_of_variables = rIndirectVariablesList.size();
    const IndexType local_size = number_of_nodes * number_of_variables;

    if (rOutput.size() != local_size) {
        rOutput.resize(local_size);
    }

    IndexType local_index = 0;
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const auto& r_node = r_geometry[i];
        for (IndexType j = 0; j < number_of_variables; ++j) {
            rOutput[local_index++] = rIndirectVariablesList[j](r_node, Step);
        }
    }

    KRATOS_CATCH("");
}

void FluidLSSVariableUtilities::CheckVariables(
    std::vector<const Variable<double>*>& rAllVariablesList,
    const std::vector<const Variable<double>*>& rCurrentVariablesList)
{
    KRATOS_TRY

    for (const auto& r_variable : rCurrentVariablesList) {
        if (r_variable) {
            const auto& p_itr = std::find_if(rAllVariablesList.begin(), rAllVariablesList.end(), [&r_variable](const Variable<double>* pVariable) {
                return *pVariable == *r_variable;
            });
            KRATOS_ERROR_IF(p_itr != rAllVariablesList.end()) << r_variable->Name() << " is already used in this utilities. Please use unique list of variables.\n";
            rAllVariablesList.push_back(r_variable);
        }
    }

    KRATOS_CATCH("");
}

void FluidLSSVariableUtilities::AddIndirectVariables(
    std::vector<IndirectVariableType>& rIndirectVariablesList,
    const std::vector<const Variable<double>*>& rVariablesList)
{
    KRATOS_TRY

    rIndirectVariablesList.clear();
    for (const auto& p_itr : rVariablesList) {
        if (p_itr) {
            IndirectVariableType temp(*p_itr);
            rIndirectVariablesList.push_back(temp);
        } else {
            IndirectVariableType temp;
            rIndirectVariablesList.push_back(temp);
        }
    }

    KRATOS_CATCH("");
}

FluidLSSVariableUtilities::FluidLSSVariableUtilities(
    const std::vector<const Variable<double>*>& rPrimalVariablePointersList,
    const std::vector<const Variable<double>*>& rPrimalFirstDerivativeVariablePointersList,
    const std::vector<const Variable<double>*>& rAdjointVariablePointersList,
    const std::vector<const Variable<double>*>& rAdjointFirstDerivativeVariablePointersList,
    const std::vector<const Variable<double>*>& rLSSVariablePointersList,
    const std::vector<const Variable<double>*>& rLSSFirstDerivativeVariablePointersList)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rPrimalVariablePointersList.size() == rPrimalFirstDerivativeVariablePointersList.size())
        << "Primal first derivative variables list size is not matching with primal variables list size [ primal first derivative variables list size = "
        << rPrimalFirstDerivativeVariablePointersList.size() << ", primal variables list size = " << rPrimalVariablePointersList.size() << " ].\n";

    KRATOS_ERROR_IF_NOT(rPrimalVariablePointersList.size() == rAdjointVariablePointersList.size())
        << "Adjoint variables list size is not matching with primal variables list size [ adjoint variables list size = "
        << rAdjointVariablePointersList.size() << ", primal variables list size = " << rPrimalVariablePointersList.size() << " ].\n";

    KRATOS_ERROR_IF_NOT(rPrimalVariablePointersList.size() == rAdjointFirstDerivativeVariablePointersList.size())
        << "Adjoint first derivative variables list size is not matching with primal variables list size [ adjoint first derivative variables list size = "
        << rAdjointFirstDerivativeVariablePointersList.size() << ", primal variables list size = " << rPrimalVariablePointersList.size() << " ].\n";

    KRATOS_ERROR_IF_NOT(rPrimalVariablePointersList.size() == rLSSVariablePointersList.size())
        << "LSS variables list size is not matching with primal variables list size [ lss variables list size = "
        << rLSSVariablePointersList.size() << ", primal variables list size = " << rPrimalVariablePointersList.size() << " ].\n";

    KRATOS_ERROR_IF_NOT(rPrimalVariablePointersList.size() == rLSSFirstDerivativeVariablePointersList.size())
        << "LSS first derivative variables list size is not matching with primal variables list size [ lss first derivative variables list size = "
        << rLSSFirstDerivativeVariablePointersList.size() << ", primal variables list size = " << rPrimalVariablePointersList.size() << " ].\n";

    // check for common variables
    std::vector<const Variable<double>*> temp;
    CheckVariables(temp, rPrimalVariablePointersList);
    CheckVariables(temp, rPrimalFirstDerivativeVariablePointersList);
    CheckVariables(temp, rAdjointVariablePointersList);
    CheckVariables(temp, rAdjointFirstDerivativeVariablePointersList);
    CheckVariables(temp, rLSSVariablePointersList);
    CheckVariables(temp, rLSSFirstDerivativeVariablePointersList);

    AddIndirectVariables(mPrimalIndirectVariablesList, rPrimalVariablePointersList);
    AddIndirectVariables(mPrimalFirstDerivativeIndirectVariablesList, rPrimalFirstDerivativeVariablePointersList);
    AddIndirectVariables(mAdjointIndirectVariablesList, rAdjointVariablePointersList);
    AddIndirectVariables(mAdjointFirstDerivativeIndirectVariablesList, rAdjointFirstDerivativeVariablePointersList);
    AddIndirectVariables(mLSSIndirectVariablesList, rLSSVariablePointersList);
    AddIndirectVariables(mLSSFirstDerivativeIndirectVariablesList, rLSSFirstDerivativeVariablePointersList);

    KRATOS_CATCH("");
}

template<class TEntityType>
void FluidLSSVariableUtilities::GetPrimalValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step) const
{
    FluidLSSVariableUtilities::GetValues(rOutput, rEntity, Step, mPrimalIndirectVariablesList);
}

template<class TEntityType>
void FluidLSSVariableUtilities::GetPrimalFirstDerivativeValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step) const
{
    FluidLSSVariableUtilities::GetValues(rOutput, rEntity, Step, mPrimalFirstDerivativeIndirectVariablesList);
}

template<class TEntityType>
void FluidLSSVariableUtilities::GetAdjointValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step) const
{
    FluidLSSVariableUtilities::GetValues(rOutput, rEntity, Step, mAdjointIndirectVariablesList);
}

template<class TEntityType>
void FluidLSSVariableUtilities::GetAdjointFirstDerivativeValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step) const
{
    FluidLSSVariableUtilities::GetValues(rOutput, rEntity, Step, mAdjointFirstDerivativeIndirectVariablesList);
}

template<class TEntityType>
void FluidLSSVariableUtilities::GetLSSValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step) const
{
    FluidLSSVariableUtilities::GetValues(rOutput, rEntity, Step, mLSSIndirectVariablesList);
}

template<class TEntityType>
void FluidLSSVariableUtilities::GetLSSFirstDerivativeValues(
    Vector& rOutput,
    const TEntityType& rEntity,
    const IndexType Step) const
{
    FluidLSSVariableUtilities::GetValues(rOutput, rEntity, Step, mLSSFirstDerivativeIndirectVariablesList);
}

const std::vector<FluidLSSVariableUtilities::IndirectVariableType>& FluidLSSVariableUtilities::GetPrimalIndirectVariablesList() const
{
    return mPrimalIndirectVariablesList;
}

const std::vector<FluidLSSVariableUtilities::IndirectVariableType>& FluidLSSVariableUtilities::GetPrimalFirstDerivativeIndirectVariablesList() const
{
    return mPrimalFirstDerivativeIndirectVariablesList;
}

const std::vector<FluidLSSVariableUtilities::IndirectVariableType>& FluidLSSVariableUtilities::GetAdjointIndirectVariablesList() const
{
    return mAdjointIndirectVariablesList;
}

const std::vector<FluidLSSVariableUtilities::IndirectVariableType>& FluidLSSVariableUtilities::GetAdjointFirstDerivativeIndirectVariablesList() const
{
    return mAdjointFirstDerivativeIndirectVariablesList;
}

const std::vector<FluidLSSVariableUtilities::IndirectVariableType>& FluidLSSVariableUtilities::GetLSSIndirectVariablesList() const
{
    return mLSSIndirectVariablesList;
}

const std::vector<FluidLSSVariableUtilities::IndirectVariableType>& FluidLSSVariableUtilities::GetLSSFirstDerivativeIndirectVariablesList() const
{
    return mLSSFirstDerivativeIndirectVariablesList;
}

// template instantiations
template void FluidLSSVariableUtilities::GetValues(Vector&, const ConditionType&, const IndexType, const std::vector<IndirectVariableType>&);
template void FluidLSSVariableUtilities::GetPrimalValues(Vector&, const ConditionType&, const IndexType) const;
template void FluidLSSVariableUtilities::GetPrimalFirstDerivativeValues(Vector&, const ConditionType&, const IndexType) const;
template void FluidLSSVariableUtilities::GetAdjointValues(Vector&, const ConditionType&, const IndexType) const;
template void FluidLSSVariableUtilities::GetAdjointFirstDerivativeValues(Vector&, const ConditionType&, const IndexType) const;
template void FluidLSSVariableUtilities::GetLSSValues(Vector&, const ConditionType&, const IndexType) const;
template void FluidLSSVariableUtilities::GetLSSFirstDerivativeValues(Vector&, const ConditionType&, const IndexType) const;

template void FluidLSSVariableUtilities::GetValues(Vector&, const ElementType&, const IndexType, const std::vector<IndirectVariableType>&);
template void FluidLSSVariableUtilities::GetPrimalValues(Vector&, const ElementType&, const IndexType) const;
template void FluidLSSVariableUtilities::GetPrimalFirstDerivativeValues(Vector&, const ElementType&, const IndexType) const;
template void FluidLSSVariableUtilities::GetAdjointValues(Vector&, const ElementType&, const IndexType) const;
template void FluidLSSVariableUtilities::GetAdjointFirstDerivativeValues(Vector&, const ElementType&, const IndexType) const;
template void FluidLSSVariableUtilities::GetLSSValues(Vector&, const ElementType&, const IndexType) const;
template void FluidLSSVariableUtilities::GetLSSFirstDerivativeValues(Vector&, const ElementType&, const IndexType) const;

} // namespace Kratos

