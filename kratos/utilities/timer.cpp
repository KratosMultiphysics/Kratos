//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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
            << std::setw(4)
            << std::setiosflags(std::ios::fixed)
            << mRepeatNumber
            << std::resetiosflags(std::ios::fixed)
            << "\t\t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(4)
            << std::uppercase
            << std::setw(4)
            << mTotalElapsedTime
            << "     \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(4)
            << std::uppercase
            << std::setw(4)
            << mMaximumTime
            << "     \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(4)
            << std::uppercase
            << std::setw(4)
            << mMinimumTime
            << "     \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(4)
            << std::uppercase
            << std::setw(4)
            << mTotalElapsedTime/static_cast<double>(mRepeatNumber)
            << "     \t" ;
        } else {
            rOStream.precision(6);
            rOStream
            << std::setw(4)
            << std::setiosflags(std::ios::fixed)
            << mRepeatNumber
            << "\t\t"
            << std::resetiosflags(std::ios::fixed)
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(4)
            << std::uppercase
            << std::setw(4)
            << mTotalElapsedTime
            << "     \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(4)
            << std::uppercase
            << std::setw(4)
            << mMaximumTime
            << "     \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(4)
            << std::uppercase
            << std::setw(4)
            << mMinimumTime
            << "     \t"
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(4)
            << std::uppercase
            << std::setw(4)
            << mTotalElapsedTime/static_cast<double>(mRepeatNumber)
            << std::resetiosflags(std::ios::scientific)
            << "     \t"
            << std::setiosflags(std::ios::fixed)
            << std::setprecision(3)
            << std::uppercase
            << std::setw(3)
            << (mTotalElapsedTime/GlobalElapsedTime)*100.00 << "%" ;
        }
    }
}

/// Default constructor.
Timer::Timer(){}

void Timer::Start(std::string const& rIntervalName)
{
    GetLabelsStackInstance().push_back(rIntervalName);
    auto full_name = CreateFullLabel();
    const auto it_internal_name = msInternalNameDatabase.find(full_name);
    if(it_internal_name == msInternalNameDatabase.end()) {
        const std::string internal_name = GetInternalName(full_name);
        msInternalNameDatabase.insert(std::pair<std::string, std::string>(full_name, internal_name));
        msTimeTable[internal_name].SetStartTime(GetTime());
        ++msCounter;
    }
    const std::string& r_name = msInternalNameDatabase[full_name];
    ContainerType::iterator it_time_data = msTimeTable.find(r_name);
    it_time_data->second.SetStartTime(GetTime());
}

void Timer::Stop(std::string const& rIntervalName)
{
    auto full_name = CreateFullLabel();
    GetLabelsStackInstance().pop_back();
    const double stop_time = GetTime();
    const std::string& r_name = msInternalNameDatabase[full_name];
    ContainerType::iterator it_time_data = msTimeTable.find(r_name);

    if(it_time_data == msTimeTable.end())
        return;

    it_time_data->second.Update(stop_time);

    if (msPrintIntervalInformation) {
        PrintIntervalInformation(r_name, it_time_data->second.GetStartTime(), stop_time);
    }
}

int Timer::SetOutputFile(std::string const& rOutputFileName)
{
    if(msOutputFile.is_open())
        msOutputFile.close();

    msOutputFile.open(rOutputFileName.c_str());

    if (msPrintIntervalInformation) {
        msOutputFile << "                                         Start      \t\tStop          \t\tElapsed" << std::endl;
    }

    return msOutputFile.is_open();
}

int Timer::CloseOutputFile()
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
    << std::setprecision(4)
    << std::uppercase
    << std::setw(4)
    << StartTime << "s\t\t"
    << std::setiosflags(std::ios::scientific)
    << std::setprecision(4)
    << std::uppercase
    << std::setw(4)
    << StopTime << "s\t\t"
    << std::setiosflags(std::ios::scientific)
    << std::setprecision(4)
    << std::uppercase
    << std::setw(4)
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
    std::size_t max_string_length = 40;
    for(auto& r_time_data : msTimeTable) {
        if (r_time_data.first.size() > max_string_length) {
            max_string_length = r_time_data.first.size();
        }
    }
    const double global_elapsed_time = ElapsedSeconds(mStartTime);
    std::string header = "   Repeat #\t\tTotal       \tMax         \tMin         \tAverage       \tTime%";
    for (std::size_t i = 0; i < max_string_length - 6; i++) {
        header.insert(0, " ");
    }
    rOStream << header << std::endl;
    for(auto& r_time_data : msTimeTable) {
        rOStream << r_time_data.first;
        for(int i =  r_time_data.first.size() ; i < static_cast<int>(max_string_length) ; i++)
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
