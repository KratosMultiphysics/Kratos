//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_DEFINE_H_INCLUDED
#define CO_SIM_IO_DEFINE_H_INCLUDED

// System includes
#include <string>
#include <unordered_map>

namespace CoSimIO {

typedef std::unordered_map<std::string, std::string> SettingsType;

enum class ControlSignal
{
    Dummy,
    BreakSolutionLoop,
    ConvergenceAchieved,

    AdvanceInTime,
    InitializeSolutionStep,
    Predict,
    SolveSolutionStep,
    FinalizeSolutionStep,
    OutputSolutionStep,

    ImportGeometry,
    ExportGeometry,
    ImportMesh,
    ExportMesh,
    ImportData,
    ExportData,
};

enum ConnectionStatus
{
    NotConnected,
    Connected,
    Disconnected,
    ConnectionError,
    DisconnectionError
};

} //namespace CoSimIO

#endif // CO_SIM_IO_DEFINE_H_INCLUDED