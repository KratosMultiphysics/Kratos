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

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidLeastSquaresShadowingUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using ConditionType = Condition;

    using ElementType = Element;

    KRATOS_CLASS_POINTER_DEFINITION(FluidLeastSquaresShadowingUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    FluidLeastSquaresShadowingUtilities(
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
    void CheckVariables(const TEntityType& rEntity) const;

    template<class TEntityType>
    void GetPrimalValues(Vector& rOutput, const TEntityType& rEntity, const IndexType Step = 0) const;

    template<class TEntityType>
    void GetPrimalFirstDerivativeValues(Vector& rOutput, const TEntityType& rEntity, const IndexType Step = 0) const;

    template<class TEntityType>
    void GetAdjointValues(Vector& rOutput, const TEntityType& rEntity, const IndexType Step = 0) const;

    template<class TEntityType>
    void GetAdjointFirstDerivativeValues(Vector& rOutput, const TEntityType& rEntity, const IndexType Step = 0) const;

    template<class TEntityType>
    void GetLSSValues(Vector& rOutput, const TEntityType& rEntity, const IndexType Step = 0) const;

    template<class TEntityType>
    void GetLSSFirstDerivativeValues(Vector& rOutput, const TEntityType& rEntity, const IndexType Step = 0) const;

    const std::vector<const Variable<double>*>& GetPrimalVariablePointersList() const;

    const std::vector<const Variable<double>*>& GetPrimalFirstDerivativeVariablePointersList() const;

    const std::vector<const Variable<double>*>& GetAdjointVariablePointersList() const;

    const std::vector<const Variable<double>*>& GetAdjointFirstDerivativeVariablePointersList() const;

    const std::vector<const Variable<double>*>& GetLSSVariablePointersList() const;

    const std::vector<const Variable<double>*>& GetLSSFirstDerivativeVariablePointersList() const;

    ///@}
private:
    ///@name Private Members
    ///@{

    const std::vector<const Variable<double>*> mPrimalVariablePointersList;
    const std::vector<const Variable<double>*> mPrimalFirstDerivativeVariablePointersList;
    const std::vector<const Variable<double>*> mAdjointVariablePointersList;
    const std::vector<const Variable<double>*> mAdjointFirstDerivativeVariablePointersList;
    const std::vector<const Variable<double>*> mLSSVariablePointersList;
    const std::vector<const Variable<double>*> mLSSFirstDerivativeVariablePointersList;

    ///@}
    ///@name Private Static Operations
    ///@{

    template<class TEntityType>
    static void GetValues(Vector& rOutput, const TEntityType& rEntity, const IndexType Step, const std::vector<const Variable<double>*>& rVariablePointersList);

    ///@}
};

///@}

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_LSS_UTILITIES_H
