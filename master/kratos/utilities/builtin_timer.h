//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre
//
//

#pragma once

// System includes
#include <chrono>

namespace Kratos
{

///@name Kratos Classes
///@{

class BuiltinTimer
{
public:
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BuiltinTimer() : mStartTime{std::chrono::steady_clock::now()}
    {
    }

    ///@}
    ///@name Operations
    ///@{

    inline double ElapsedSeconds() const
    {
        using namespace std::chrono;
        const auto current_time = steady_clock::now();
        return duration_cast<duration<double>>(current_time - mStartTime).count();
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    const std::chrono::steady_clock::time_point mStartTime;

    ///@}
};

///@} // Kratos Classes

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const BuiltinTimer& rThis)
{
    const double elapsed_secs = rThis.ElapsedSeconds();

    int minutes = elapsed_secs / 60;
    int hours = minutes / 60;
    minutes = minutes % 60;
    double seconds = elapsed_secs - hours*3600.0 - minutes*60.0;

    if (hours > 0) {
        rOStream << hours << " [h] ";
    }
    if (minutes > 0) {
        rOStream << minutes << " [m] ";
        // remove decimals if over a minute
        seconds = int(seconds);
    }
    rOStream << seconds << " [s]";
    return rOStream;

}

}  // namespace Kratos.
