//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre
//
//

#if !defined(KRATOS_BUILTIN_TIMER_H_INCLUDED)
#define  KRATOS_BUILTIN_TIMER_H_INCLUDED

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

}  // namespace Kratos.

#endif // KRATOS_BUILTIN_TIMER_H_INCLUDED  defined 
