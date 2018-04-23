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
            const unsigned int interval_size = settings["interval"].size()-1;
            if(settings["interval"][interval_size].IsString() )
            {
                if(settings["interval"][interval_size].GetString() == std::string("End"))
                    settings["interval"][interval_size].SetDouble(1e30);
                else
                    KRATOS_ERROR << "the second value of interval can be \"End\" or a number, interval currently: \n"+settings["interval"].PrettyPrintJsonString();
            }
        }
        else
        {
            Parameters defaults(R"( {"default_interval": [0.0, 1e30]} )");
            settings.AddValue("interval", defaults["default_interval"]);
        }

        mInterval = settings["interval"].GetVector();

        KRATOS_CATCH("");
    }

    double GetIntervalBegin()
    {
        return mInterval[0];
    }

    double GetIntervalEnd()
    {
        return mInterval[mInterval.size() - 1];
    }

    double GetSubIntervalBegin(int n)
    {
        KRATOS_DEBUG_ERROR_IF((n<1)||(n>mInterval.size()-1)) << "wrong interval number: " << n << std::endl;
        return mInterval[n-1];
    }

    double GetSubIntervalEnd(int n)
    {
        KRATOS_DEBUG_ERROR_IF((n<1) || (n>mInterval.size() - 1)) << "wrong interval number: " << n << std::endl;
        return mInterval[n];
    }

    bool IsInInterval(double time )
    {
        const size_t interval_size = mInterval.size() - 1;
        const double eps = std::max(1e-14*mInterval[0], 1e-30);
        if(time > mInterval[0] -eps && time < mInterval[interval_size] +eps)
            return true;
        else
            return false;
    }

    int GetSubInterval(double time)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(IsInInterval(time)) << "the current time: " << time << " is not inside the interval - Check it with the function 'IsInInterval'" << std::endl;
        for (unsigned int i = 1; i < mInterval.size(); ++i)
        {
            if (mInterval[i] > time)
            {
                return i;
            }
        }
        return mInterval.size()-1;
    }

    int GetSubIntervalNb()
    {
        return mInterval.size()-1;
    }

    void RemoveSubInterval()
    {
        double interval_min = mInterval[0];
        double interval_max = mInterval[mInterval.size() - 1];
        mInterval.resize(2);
        for (unsigned int i = 1; i < mInterval.size()-1; ++i)
        {
            //mInterval.
        }
        mInterval[0] = interval_min;
        mInterval[1] = interval_max;
    }

private:
    Vector mInterval;
};


}

#endif // KRATOS_INTERVAL_UTILITY_H_INCLUDED
