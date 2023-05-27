//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"

namespace Kratos {

template <class TDataType>
class KRATOS_API(KRATOS_CORE) FixedIntervalRecurringEventUtility
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FixedIntervalRecurringEventUtility);

    ///@}
    ///@name Life cycle
    ///@{

    FixedIntervalRecurringEventUtility(
        const TDataType InitialValue,
        const TDataType Interval);

    ///@}
    ///@name Public operations
    ///@{

    bool IsEventExpected(const TDataType CurrentValue) const;

    void ScheduleNextEvent(const TDataType CurrentValue);

    ///@}
    ///@name Public input output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const TDataType mInterval;

    TDataType mNextEventValue;

    ///@}
};

/// output stream function
template <class TDataType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const FixedIntervalRecurringEventUtility<TDataType>& rThis)
{
    rOStream << rThis.Info();
    return rOStream;
}

} // namespace Kratos
