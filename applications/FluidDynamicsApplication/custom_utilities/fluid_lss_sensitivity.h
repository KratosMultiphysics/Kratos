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

#if !defined(KRATOS_FLUID_LSS_SENSITIVITY_H)
#define KRATOS_FLUID_LSS_SENSITIVITY_H

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "response_functions/adjoint_response_function.h"

// Application includes
#include "custom_utilities/fluid_adjoint_slip_utilities.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidLSSSensitivity
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FluidLSSSensitivity);

    ///@}
    ///@name Life Cycle
    ///@{

    FluidLSSSensitivity() {}

    virtual ~FluidLSSSensitivity() = default;

    ///@}
    ///@name Operations
    ///@{

    virtual const Variable<double>& GetDerivativeVariable() const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling the base class GetDerivativeVariable. Please call the derrived class same method.\n";

        return Variable<double>::StaticObject();
        KRATOS_CATCH("");
    }

    virtual void CalculateResidualSensitivity(
        Vector& rOutput,
        Element& rElement,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling the base class CalculateResidualSensitivity. Please call the derrived class same method.\n";

        KRATOS_CATCH("");
    }

    virtual void CalculateResidualSensitivity(
        Vector& rOutput,
        Condition& rCondition,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling the base class CalculateResidualSensitivity. Please call the derrived class same method.\n";

        KRATOS_CATCH("");
    }

    virtual double CalculateResponseSensitivity(
        Element& rElement,
        AdjointResponseFunction& rResponseFunction,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling the base class CalculateResponseSensitivity. Please call the derrived class same method.\n";

        return 0.0;

        KRATOS_CATCH("");
    }

    virtual double CalculateResponseSensitivity(
        Condition& rCondition,
        AdjointResponseFunction& rResponseFunction,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling the base class CalculateResponseSensitivity. Please call the derrived class same method.\n";

        return 0.0;

        KRATOS_CATCH("");
    }

    ///@}
};

///@}

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_LSS_SENSITIVITY_H
