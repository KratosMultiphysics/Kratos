//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_INTERVAL_UTILITY_H_INCLUDED)
#define  KRATOS_INTERVAL_UTILITY_H_INCLUDED

#include "includes/define.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
namespace Detail
{

/**
 *  @brief Utility class for membership tests on a 1D interval.
 *
 *  @note This class template has specializations for @a double and @a int but is not implemented for other types.
 *  @ingroup KratosCore
 */
template <class TValue>
class IntervalUtility
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(IntervalUtility);

    /// Default constructor initializing boundaries to "Begin" and "End".
    IntervalUtility();

    /**
     *  @brief Construct from parameters containing "interval".
     *
     *  @details "interval" is expected as an array with exactly 2 items, defining
     *           the begin and end of the interval respectively. The first item can
     *           either be a numeric value or "Begin" (setting the smallest representable
     *           value of @a TValue), while the second one can be a numeric value or
     *           "End" (setting the highest representable value of @a TValue).
     *
     *  @note String values ("Begin" and "End") are replaced with their numeric counterparts
     *        in the input @a Settings.
     *
     *  @note If "interval" is not in @a Settings, an "interval" with values corresponding
     *        to "Begin" and "End" are added to it.
     *
     *  @note Other parameters in @a Settings are not checked.
     */
    IntervalUtility(Parameters Settings);

    IntervalUtility(IntervalUtility&& rOther) = default;

    IntervalUtility(const IntervalUtility& rOther) = default;

    TValue GetIntervalBegin() const noexcept;

    TValue GetIntervalEnd() const noexcept;

    /**
     *  @brief Check whether the input value is within the defined interval [Begin, End).
     *
     *  @details This member has explicit specializations for different types
     *           that have slight variations in behaviour around the interval
     *           boundaries. Check the individual specializations for the exact
     *           behaviour.
     */
    bool IsInInterval(TValue Value) const noexcept;

    static Parameters GetDefaultParameters();

private:
    TValue mBegin;

    TValue mEnd;
};

} // namespace Detail

/// A class providing membership tests on 1D rational intervals (eg.: time intervals).
using IntervalUtility = Detail::IntervalUtility<double>;

/// A class providing membership tests on 1D integer intervals (eg.: step intervals).
using DiscreteIntervalUtility = Detail::IntervalUtility<int>;

} // namespace Kratos

#endif // KRATOS_INTERVAL_UTILITY_H_INCLUDED
