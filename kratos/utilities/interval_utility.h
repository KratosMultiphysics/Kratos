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

#include <cmath>
#include "includes/define.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

/**this function manages intervals. It aims at being used within processes
*
*/
class KRATOS_API(KRATOS_CORE) IntervalUtility
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

private:
    double mIntervalBegin;
    double mIntervalEnd;
};


/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const IntervalUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : ";
    rThis.PrintData(rOStream);
    return rOStream;
}


}

#endif // KRATOS_INTERVAL_UTILITY_H_INCLUDED
