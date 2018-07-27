//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//                    
//

#if !defined(KRATOS_TIMER_H_INCLUDED )
#define  KRATOS_TIMER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

// External includes
#include <boost/timer.hpp> // to be removed after replacing the boost timers with Kratos timer.

// Project includes
#include "includes/define.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class Timer
 * @ingroup KratosCore
 * @brief This utility can be used to compute the time employed on computations
 * @author Pooyan Dadvand
 * @author Riccardo Rossi
 * @todo The boost::timer dependency is not in this file but in the files that include it
 * @todo Use logger
 */
class KRATOS_API(KRATOS_CORE) Timer
{
    class TimerData
    {
        int mRepeatNumber;
        double mStartTime;
        double mTotalElapsedTime;
        double mMaximumTime;
        double mMinimumTime;
    public:
        TimerData() : mRepeatNumber(int()), mStartTime(double()), mTotalElapsedTime(double()), mMaximumTime(double()), mMinimumTime(double()) {}
        double GetStartTime()
        {
            return mStartTime;
        }
        void SetStartTime(double StartTime)
        {
            mStartTime = StartTime;
        }
        void Update(double StopTime)
        {
            double elapsed = StopTime - mStartTime;
            if(mRepeatNumber == 0)
                mMinimumTime = elapsed;
            mTotalElapsedTime += elapsed;
            if(mMaximumTime < elapsed)
                mMaximumTime = elapsed;

            if((mMinimumTime > elapsed))
                mMinimumTime = elapsed;

            mRepeatNumber++;

        }
        /// Print object's data.
        void PrintData(std::ostream& rOStream, double GlobalElapsedTime = -1.00) const
        {
            if(mRepeatNumber != 0)
            {
                if(GlobalElapsedTime <= 0.00)
                    rOStream << mRepeatNumber << " \t" << mTotalElapsedTime << "s     \t" << mMaximumTime << "s     \t" << mMinimumTime << "s     \t" << mTotalElapsedTime/static_cast<double>(mRepeatNumber) << "s     \t" ;
                else
                    rOStream << mRepeatNumber << " \t" << mTotalElapsedTime << "s     \t" << mMaximumTime << "s     \t" << mMinimumTime << "s     \t" << mTotalElapsedTime/static_cast<double>(mRepeatNumber) << "s     \t" << (mTotalElapsedTime/GlobalElapsedTime)*100.00 << "%" ;
            }

        }
    };

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Timer
    KRATOS_CLASS_POINTER_DEFINITION(Timer);

    typedef double TimeType;


    typedef std::map<std::string, TimerData> ContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Timer();

    /// Destructor.
    virtual ~Timer()
    {

    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static void Start(std::string const& IntervalName)
    {
        msTimeTable[IntervalName].SetStartTime(GetTime());
    }

    static void Stop(std::string const& IntervalName)
    {
        double stop_time = GetTime();
        ContainerType::iterator i_time_data = msTimeTable.find(IntervalName);

        if(i_time_data == msTimeTable.end())
            return;
        /* 	  KRATOS_THROW_ERROR(std::logical_error, "Stopping a not running time interval: ", IntervalName); */

        i_time_data->second.Update(stop_time);

        PrintIntervalInformation(IntervalName, i_time_data->second.GetStartTime(), stop_time);
    }

    static inline double GetTime()
    {
#ifndef _OPENMP
        return std::clock()/static_cast<double>(CLOCKS_PER_SEC);
#else
        return omp_get_wtime();
#endif
    }

    ///@}
    ///@name Access
    ///@{

    static int SetOuputFile(std::string const& OutputFileName)
    {
        if(msOutputFile.is_open())
            msOutputFile.close();

        msOutputFile.open(OutputFileName.c_str());

        msOutputFile << "                                         Start   \tStop     \tElapsed " << std::endl;

        return msOutputFile.is_open();
    }

    static int CloseOuputFile()
    {
        if(msOutputFile.is_open())
            msOutputFile.close();

        return msOutputFile.is_open();
    }

    static bool GetPrintOnScreen()
    {
        return msPrintOnScreen;
    }

    static void SetPrintOnScreen(bool const PrintOnScreen)
    {
        msPrintOnScreen = PrintOnScreen;
    }


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    static void PrintIntervalInformation(std::string const& IntervalName, double StartTime, double StopTime)
    {
        if(msOutputFile.is_open())
        {
            msOutputFile << IntervalName << " ";

            for(int i = IntervalName.size() + 1 ; i < 40 ; i++)
                msOutputFile << ".";

            msOutputFile << " " << StartTime << "s     \t" << StopTime << "s     \t" << StopTime - StartTime <<"s" << std::endl;
        }
        else if(msPrintOnScreen)
        {
            std::cout << IntervalName << " ";

            for(int i = IntervalName.size() + 1 ; i < 40 ; i++)
                std::cout << ".";

            std::cout << " " << StartTime << "s     \t" << StopTime << "s     \t" << StopTime - StartTime <<"s" << std::endl;
        }
    }

    static void PrintTimingInformation()
    {
        if(msOutputFile.is_open())
            PrintTimingInformation(msOutputFile);
        else if(msPrintOnScreen)
            PrintTimingInformation(std::cout);
    }

    static void PrintTimingInformation(std::ostream& rOStream)
    {
        double global_elapsed_time = GetTime() - msGlobalStart;
        rOStream << "                                 Repeat # \tTotal     \tMax     \tMin     \tAverage     \t%" << std::endl;
        for(ContainerType::iterator i_time_data = msTimeTable.begin() ; i_time_data != msTimeTable.end() ; i_time_data++)
        {
            rOStream << i_time_data->first;
            for(int i =  i_time_data->first.size() + 1 ; i < 40 ; i++)
                rOStream << ".";

            rOStream << " ";
            i_time_data->second.PrintData(rOStream, global_elapsed_time);
            rOStream << std::endl;
        }
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Timer";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        PrintTimingInformation(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:

    ///@name Static Member Variables
    ///@{

    static ContainerType msTimeTable;

    static std::ofstream msOutputFile;

    static bool msPrintOnScreen;

    static double msGlobalStart;



    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    Timer& operator=(Timer const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    /*       Timer(Timer const& rOther); */


    ///@}

}; // Class Timer

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
// 				    Timer& rThis){}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Timer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_TIMER_H_INCLUDED  defined 


