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
//                   Vicente Mataix Ferrandiz
//

// System includes
#include <iomanip>
#include <sstream>

// External includes

// Project includes
#include "utilities/timer.h"
#include "input_output/logger.h"

namespace Kratos
{
void Timer::TimerData::PrintData(
    std::ostream& rOStream,
    double GlobalElapsedTime
    ) const
{
    if(mRepeatNumber != 0) {
        if(GlobalElapsedTime <= 0.0) {
            rOStream.precision(6);
            rOStream
            << mRepeatNumber
            << "\t\t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(6)
            << std::uppercase
            << std::setw(6)
            << mTotalElapsedTime
            << "s    \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(6)
            << std::uppercase
            << std::setw(6)
            << mMaximumTime
            << "s    \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(6)
            << std::uppercase
            << std::setw(6)
            << mMinimumTime
            << "s    \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(6)
            << std::uppercase
            << std::setw(6)
            << mTotalElapsedTime/static_cast<double>(mRepeatNumber)
            << "s    \t" ;
        } else {
            rOStream.precision(6);
            rOStream
            << mRepeatNumber
            << "\t\t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(6)
            << std::uppercase
            << std::setw(6)
            << mTotalElapsedTime
            << "s    \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(6)
            << std::uppercase
            << std::setw(6)
            << mMaximumTime
            << "s    \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(6)
            << std::uppercase
            << std::setw(6)
            << mMinimumTime
            << "s    \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(6)
            << std::uppercase
            << std::setw(6)
            << mTotalElapsedTime/static_cast<double>(mRepeatNumber)
            << "s    \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(6)
            << std::uppercase
            << std::setw(6)
            << (mTotalElapsedTime/GlobalElapsedTime)*100.00 << "%" ;
        }
    }
}

/// Default constructor.
Timer::Timer(){}

void Timer::Start(std::string const& rIntervalName)
{
    const auto it_internal_name = msInternalNameDatabase.find(rIntervalName);
    if(it_internal_name == msInternalNameDatabase.end()) {
        const std::string internal_name = GetInternalName(rIntervalName);
        msInternalNameDatabase.insert(std::pair<std::string, std::string>(rIntervalName, internal_name));
        msTimeTable[internal_name].SetStartTime(GetTime());
        ++msCounter;
    }
}

void Timer::Stop(std::string const& rIntervalName)
{
    const double stop_time = GetTime();
    const std::string& r_name = msInternalNameDatabase[rIntervalName];
    ContainerType::iterator it_time_data = msTimeTable.find(r_name);

    if(it_time_data == msTimeTable.end())
        return;

    it_time_data->second.Update(stop_time);

    if (msPrintIntervalInformation) {
        PrintIntervalInformation(r_name, it_time_data->second.GetStartTime(), stop_time);
    }
}

int Timer::SetOuputFile(std::string const& rOutputFileName)
{
    if(msOutputFile.is_open())
        msOutputFile.close();

    msOutputFile.open(rOutputFileName.c_str());

    if (msPrintIntervalInformation) {
        msOutputFile << "                                         Start      \t\tStop          \t\tElapsed" << std::endl;
    }

    return msOutputFile.is_open();
}

int Timer::CloseOuputFile()
{
    if(msOutputFile.is_open())
        msOutputFile.close();

    return msOutputFile.is_open();
}

bool Timer::GetPrintOnScreen()
{
    return msPrintOnScreen;
}

void Timer::SetPrintOnScreen(bool const PrintOnScreen)
{
    msPrintOnScreen = PrintOnScreen;
}

bool Timer::GetPrintIntervalInformation()
{
    return msPrintIntervalInformation;
}

void Timer::SetPrintIntervalInformation(bool const PrintIntervalInformation)
{
    msPrintIntervalInformation = PrintIntervalInformation;
}

void Timer::PrintIntervalInformation(std::string const& rIntervalName, const double StartTime, const double StopTime)
{
    const std::string& r_name = msInternalNameDatabase[rIntervalName];
    if (msOutputFile.is_open())
        PrintIntervalInformation(msOutputFile, r_name, StartTime, StopTime);
    if(msPrintOnScreen)
        PrintIntervalInformation(std::cout, r_name, StartTime, StopTime);
}

void Timer::PrintIntervalInformation(std::ostream& rOStream, std::string const& rIntervalName, const double StartTime, const double StopTime)
{
    rOStream << rIntervalName << " ";

    for(int i = rIntervalName.size() + 1 ; i < 40 ; i++)
        rOStream << ".";

    rOStream.precision(6);
    rOStream << " "
    << std::setiosflags(std::ios::scientific)
    << std::setprecision(6)
    << std::uppercase
    << std::setw(6)
    << StartTime << "s\t\t"
    << std::setiosflags(std::ios::scientific)
    << std::setprecision(6)
    << std::uppercase
    << std::setw(6)
    << StopTime << "s\t\t"
    << std::setiosflags(std::ios::scientific)
    << std::setprecision(6)
    << std::uppercase
    << std::setw(6)
    << StopTime - StartTime <<"s" << std::endl;
}

void Timer::PrintTimingInformation()
{
    if(msOutputFile.is_open())
        PrintTimingInformation(msOutputFile);
    if(msPrintOnScreen)
        PrintTimingInformation(std::cout);
}

void Timer::PrintTimingInformation(std::ostream& rOStream)
{
    const double global_elapsed_time = ElapsedSeconds(mStartTime);
    rOStream << "                                 Repeat #\t\tTotal           \tMax             \tMin             \tAverage           \tTime%" << std::endl;
    for(auto& r_time_data : msTimeTable) {
        rOStream << r_time_data.first;
        for(int i =  r_time_data.first.size() + 1 ; i < 40 ; i++)
            rOStream << ".";

        rOStream << " ";
        r_time_data.second.PrintData(rOStream, global_elapsed_time);
        rOStream << std::endl;
    }
}

std::string Timer::GetInternalName(const std::string& rName)
{
    std::string initial_string = std::to_string(msCounter);
    const std::size_t size_string = initial_string.size();
    for (std::size_t i = 0; i < (NumberOfZeros - size_string); ++i) {
        initial_string = "0" + initial_string;
    }
    std::string internal_name = initial_string + "." + rName;

    return internal_name;
}

Timer::InternalNameDatabaseType Timer::msInternalNameDatabase;
Timer::ContainerType Timer::msTimeTable;
std::ofstream Timer::msOutputFile;
std::size_t Timer::msCounter = 0;
bool Timer::msPrintOnScreen = false;
bool Timer::msPrintIntervalInformation = true;
const std::chrono::steady_clock::time_point Timer::mStartTime = std::chrono::steady_clock::now();

} /// namespace Kratos

