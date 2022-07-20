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

template <class TValue>
class IntervalUtility
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(IntervalUtility);

    IntervalUtility(Parameters Settings);

    TValue GetIntervalBegin() const noexcept;

    TValue GetIntervalEnd() const noexcept;

    bool IsInInterval(TValue Value) const noexcept;

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
