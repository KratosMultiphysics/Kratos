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
class IntervalUtility
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(IntervalUtility);

    IntervalUtility(  Parameters settings )
    {
        KRATOS_TRY

        if(settings.Has("interval"))
        {
            if(settings["interval"][1].IsString() )
            {
                if(settings["interval"][1].GetString() == std::string("End"))
                    settings["interval"][1].SetDouble(1e30);
                else
                    KRATOS_ERROR << "the second value of interval can be \"End\" or a number, interval currently: \n"+settings["interval"].PrettyPrintJsonString();
            }
        }
        else
        {
            Parameters defaults(R"( {"default_interval": [0.0, 1e30]} )");
            settings.AddValue("interval", defaults["default_interval"]);
        }

        mIntervalBegin = settings["interval"][0].GetDouble();
        mIntervalEnd = settings["interval"][1].GetDouble();

        KRATOS_CATCH("");
    }

    const double GetIntervalBegin() const
    {
        return mIntervalBegin;
    }

    const double GetIntervalEnd() const
    {
        return mIntervalEnd;
    }

    bool IsInInterval(double time )
    {
        const double eps = std::max(1e-14*mIntervalBegin, 1e-30);
        if(time > mIntervalBegin-eps && time < mIntervalEnd+eps)
            return true;
        else
            return false;
    }

    /// Turn back information as a string.
    std::string Info() const
    {
        return "IntervalUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
        rOStream << "[" << GetIntervalBegin() << ", " << GetIntervalEnd() << "]";
    }

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
