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
//

// System includes

// External includes

// Project includes
#include "utilities/timer.h"
#include "input_output/logger.h"

namespace Kratos
{
/// Default constructor.
Timer::Timer(){}

void Timer::Start(std::string const& rIntervalName)
{
    msTimeTable[rIntervalName].SetStartTime(GetTime());
}

void Timer::Stop(std::string const& rIntervalName)
{
    const double stop_time = GetTime();
    ContainerType::iterator it_time_data = msTimeTable.find(rIntervalName);

    if(it_time_data == msTimeTable.end())
        return;

    it_time_data->second.Update(stop_time);

    PrintIntervalInformation(rIntervalName, it_time_data->second.GetStartTime(), stop_time);
}

int Timer::SetOuputFile(std::string const& rOutputFileName)
{
    if(msOutputFile.is_open())
        msOutputFile.close();

    msOutputFile.open(rOutputFileName.c_str());

    msOutputFile << "                                         Start   \tStop     \tElapsed " << std::endl;

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

void Timer::PrintIntervalInformation(std::string const& rIntervalName, const double StartTime, const double StopTime)
{
    if (msOutputFile.is_open()) {
        msOutputFile << rIntervalName << " ";

        for(int i = rIntervalName.size() + 1 ; i < 40 ; i++)
            msOutputFile << ".";

        msOutputFile << " " << StartTime << "s     \t" << StopTime << "s     \t" << StopTime - StartTime <<"s" << std::endl;
    }

    if(msPrintOnScreen) {
        KRATOS_INFO("") << rIntervalName << " ";

        for(int i = rIntervalName.size() + 1 ; i < 40 ; i++)
            KRATOS_INFO("") << ".";


        KRATOS_INFO("") << " " << StartTime << "s     \t" << StopTime << "s     \t" << StopTime - StartTime <<"s" << std::endl;
    }
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
    rOStream << "                                 Repeat # \tTotal     \tMax     \tMin     \tAverage     \t%" << std::endl;
    for(auto& r_time_data : msTimeTable) {
        rOStream << r_time_data.first;
        for(int i =  r_time_data.first.size() + 1 ; i < 40 ; i++)
            rOStream << ".";

        rOStream << " ";
        r_time_data.second.PrintData(rOStream, global_elapsed_time);
        rOStream << std::endl;
    }
}

Timer::ContainerType Timer::msTimeTable;
std::ofstream Timer::msOutputFile;
bool Timer::msPrintOnScreen = false;
const std::chrono::steady_clock::time_point Timer::mStartTime = std::chrono::steady_clock::now();

} /// namespace Kratos

