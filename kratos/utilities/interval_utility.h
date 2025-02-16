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

#pragma once

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/schema.hpp" // Schema

namespace Kratos {

/**this function manages intervals. It aims at being used within processes
*
*/
class KRATOS_API(KRATOS_CORE) IntervalUtility final
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(IntervalUtility);

    /// @brief Constructor with parameters
    /// @param Settings
    IntervalUtility(Parameters Settings);

    /// @brief Get the initial time of the interval
    double GetIntervalBegin() const;

    /// @brief Get the final time of the interval
    double GetIntervalEnd() const;

    /// @brief Check if the time is in interval
    /// @param Time
    bool IsInInterval(double Time);

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

    Schema GetSchema() const;

private:
    double mIntervalBegin;
    double mIntervalEnd;
};


std::ostream& operator << (std::ostream& rOStream, const IntervalUtility& rThis);

} // namespace Kratos
