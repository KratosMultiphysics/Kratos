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

#if !defined(KRATOS_FLUID_LSS_UTILITIES_H)
#define KRATOS_FLUID_LSS_UTILITIES_H

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/element.h"

// Application includes
#include "custom_utilities/indirect_variable.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidLSSVariableUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using ConditionType = Condition;

    using ElementType = Element;

    using IndirectVariableType = IndirectVariable<double>;

    KRATOS_CLASS_POINTER_DEFINITION(FluidLSSVariableUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    FluidLSSVariableUtilities(
        const std::vector<const Variable<double>*>& rPrimalVariablePointersList,
        const std::vector<const Variable<double>*>& rPrimalFirstDerivativeVariablePointersList,
        const std::vector<const Variable<double>*>& rAdjointVariablePointersList,
        const std::vector<const Variable<double>*>& rAdjointFirstDerivativeVariablePointersList,
        const std::vector<const Variable<double>*>& rLSSVariablePointersList,
        const std::vector<const Variable<double>*>& rLSSFirstDerivativeVariablePointersList);

    ///@}
    ///@name Operations
    ///@{

    template<class TEntityType>
    void GetPrimalValues(
        Vector& rOutput,
        const TEntityType& rEntity,
        const IndexType Step = 0) const;

    template<class TEntityType>
    void GetPrimalFirstDerivativeValues(
        Vector& rOutput,
        const TEntityType& rEntity,
        const IndexType Step = 0) const;

    template<class TEntityType>
    void GetAdjointValues(
        Vector& rOutput,
        const TEntityType& rEntity,
        const IndexType Step = 0) const;

    template<class TEntityType>
    void GetAdjointFirstDerivativeValues(
        Vector& rOutput,
        const TEntityType& rEntity,
        const IndexType Step = 0) const;

    template<class TEntityType>
    void GetLSSValues(
        Vector& rOutput,
        const TEntityType& rEntity,
        const IndexType Step = 0) const;

    template<class TEntityType>
    void GetLSSFirstDerivativeValues(
        Vector& rOutput,
        const TEntityType& rEntity,
        const IndexType Step = 0) const;

    const std::vector<IndirectVariableType>& GetPrimalIndirectVariablesList() const;

    const std::vector<IndirectVariableType>& GetPrimalFirstDerivativeIndirectVariablesList() const;

    const std::vector<IndirectVariableType>& GetAdjointIndirectVariablesList() const;

    const std::vector<IndirectVariableType>& GetAdjointFirstDerivativeIndirectVariablesList() const;

    const std::vector<IndirectVariableType>& GetLSSIndirectVariablesList() const;

    const std::vector<IndirectVariableType>& GetLSSFirstDerivativeIndirectVariablesList() const;

    ///@}
private:
    ///@name Private Members
    ///@{

    std::vector<IndirectVariableType> mPrimalIndirectVariablesList;
    std::vector<IndirectVariableType> mPrimalFirstDerivativeIndirectVariablesList;
    std::vector<IndirectVariableType> mAdjointIndirectVariablesList;
    std::vector<IndirectVariableType> mAdjointFirstDerivativeIndirectVariablesList;
    std::vector<IndirectVariableType> mLSSIndirectVariablesList;
    std::vector<IndirectVariableType> mLSSFirstDerivativeIndirectVariablesList;

    ///@}
    ///@name Private Static Operations
    ///@{

    template<class TEntityType>
    static void GetValues(
        Vector& rOutput,
        const TEntityType& rEntity,
        const IndexType Step,
        const std::vector<IndirectVariableType>& rIndirectVariablesList);

    static void CheckVariables(
        std::vector<const Variable<double>*>& rAllVariablesList,
        const std::vector<const Variable<double>*>& rCurrentVariablesList);

    static void AddIndirectVariables(
        std::vector<IndirectVariableType>& rIndirectVariablesList,
        const std::vector<const Variable<double>*>& rVariablesList);

    ///@}
};

///@}

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_LSS_UTILITIES_H
