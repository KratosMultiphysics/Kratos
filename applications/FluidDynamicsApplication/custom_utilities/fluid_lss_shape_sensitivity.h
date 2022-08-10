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

#if !defined(KRATOS_FLUID_LSS_SHAPE_SENSITIVITY_H)
#define KRATOS_FLUID_LSS_SHAPE_SENSITIVITY_H

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/kratos_parameters.h"
#include "response_functions/adjoint_response_function.h"

// Application includes
#include "custom_utilities/fluid_lss_sensitivity.h"
#include "custom_utilities/fluid_adjoint_slip_utilities.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidLSSShapeSensitivity : public FluidLSSSensitivity
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = FluidLSSSensitivity;

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(FluidLSSShapeSensitivity);

    ///@}
    ///@name Life Cycle
    ///@{

    FluidLSSShapeSensitivity(
        Parameters Settings,
        const IndexType Dimension
    );

    ~FluidLSSShapeSensitivity() override = default;

    ///@}
    ///@name Operations
    ///@{

    const Variable<double>& GetDerivativeVariable() const override;

    void CalculateResidualSensitivity(
        Vector& rOutput,
        Element& rElement,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo) override;

    void CalculateResidualSensitivity(
        Vector& rOutput,
        Condition& rCondition,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo) override;

    double CalculateResponseSensitivity(
        Element& rElement,
        AdjointResponseFunction& rResponse,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo) override;

    double CalculateResponseSensitivity(
        Condition& rCondition,
        AdjointResponseFunction& rResponse,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo) override;

    ///@}
private:
    ///@name Private Classes
    ///@{

    struct TLS
    {
        Vector mResiduals;
        Matrix mAuxMatrix;
        Matrix mResidualShapeDerivatives;
        Matrix mRotatedResidualShapeDerivatives;
        Vector mResponseShapeDerivatives;
        std::vector<IndexType> mDerivativeNodeIds;
    };

    ///@}
    ///@name Private Members
    ///@{

    const IndexType mDimension;

    IndexType mShapeDerivativeNodeId;

    IndexType mShapeDerivativeDirection;

    std::vector<TLS> mTLS;

    ///@}
    ///@name Private Operations
    ///@{

    template<class TEntityType>
    void CalculateResidualSensitivity(
        Vector& rOutput,
        TEntityType& rEntity,
        Vector& rResiduals,
        Matrix& rAuxMatrix,
        Matrix& rResidualShapeDerivatives,
        Matrix& rRotatedResidualShapeDerivatives,
        std::vector<IndexType>& rDerivativeNodeIds,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo) const;

    template<class TEntityType>
    double CalculateResponseSensitivity(
        TEntityType& rEntity,
        Matrix& rResidualShapeDerivatives,
        Vector& rResponseShapeDerivatives,
        AdjointResponseFunction& rResponseFunction,
        const ProcessInfo& rProcessInfo) const;


    ///@}
};

///@}

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_LSS_SHAPE_SENSITIVITY_H
