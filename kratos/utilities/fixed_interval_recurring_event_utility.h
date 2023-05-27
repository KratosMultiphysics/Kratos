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

/**
 * @class FixedIntervalRecurringEventUtility
 * @ingroup KratosCore
 * @brief Utility class to handle fixed interval recurring events.
 * @tparam TDataType The data type used for event values.
 */
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

    /**
     * @brief Constructor.
     * @param InitialValue The initial value of the event.
     * @param Interval The time interval between events.
     */
    FixedIntervalRecurringEventUtility(
        const TDataType InitialValue,
        const TDataType Interval);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Checks if an event is expected based on the current value.
     * @param CurrentValue The current value.
     * @return true if an event is expected, false otherwise.
     */
    bool IsEventExpected(const TDataType CurrentValue) const;

    /**
     * @brief Schedules the next event based on the current value.
     * @param CurrentValue The current value.
     */
    void ScheduleNextEvent(const TDataType CurrentValue);

    ///@}
    ///@name Public input output
    ///@{

    /**
     * @brief Returns the information of the utility as a string.
     * @return The utility information as a string.
     */
    std::string Info() const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const TDataType mInterval; /// The time interval between events.
    TDataType mNextEventValue; /// The value of the next event.

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
