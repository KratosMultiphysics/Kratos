// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

#ifndef KRATOS_CO_SIM_IO_IMPL_H_INCLUDED
#define KRATOS_CO_SIM_IO_IMPL_H_INCLUDED

// Optional includes
#ifdef KRATOS_CO_SIM_IO_ENABLE_SOCKETS
#include "co_sim_sockets_comm.h"
#endif /* KRATOS_CO_SIM_IO_ENABLE_SOCKETS */


#ifdef KRATOS_CO_SIM_IO_ENABLE_MPI
#include "co_sim_mpi_comm.h"
#endif /* KRATOS_CO_SIM_IO_ENABLE_MPI */

// System includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

// Project includes
#include "co_sim_file_comm.h"


namespace CoSim {

inline CoSimIO::CoSimIO(const std::string& rName, const std::string& rSettingsFileName, const bool IsConnectionMaster)
    : CoSimIO::CoSimIO(rName, Internals::ReadSettingsFile(rSettingsFileName), IsConnectionMaster) { } // forwarding constructor call

inline CoSimIO::CoSimIO(const std::string& rName, SettingsType rSettings, const bool IsConnectionMaster)
{
    Initialize(rName, rSettings, IsConnectionMaster);
}

inline bool CoSimIO::Connect()
{
    return mpComm->Connect();
}

inline bool CoSimIO::Disconnect()
{
    return mpComm->Disconnect();
}


inline void CoSimIO::SendControlSignal(const Internals::ControlSignal Signal, const std::string& rIdentifier)
{
    mpComm->SendControlSignal(Signal, rIdentifier);

}
inline Internals::ControlSignal CoSimIO::RecvControlSignal(std::string& rIdentifier)
{
    return mpComm->RecvControlSignal(rIdentifier);
}


inline void CoSimIO::RegisterAdvanceInTime(double (*pFuncPtr)(double))
{
    mpAdvInTime = pFuncPtr;
}

inline void CoSimIO::RegisterInitializeSolutionStep(void (*pFuncPtr)())
{
    mpInitSolStep = pFuncPtr;
}

inline void CoSimIO::RegisterSolveSolutionStep(void (*pFuncPtr)())
{
    mpSolSolStep = pFuncPtr;
}

inline void CoSimIO::RegisterFinalizeSolutionStep(void (*pFuncPtr)())
{
    mpFinSolStep = pFuncPtr;
}

inline void CoSimIO::RegisterImportGeometry(void (*pFuncPtr)(const std::string&))
{
    mpImportGeom = pFuncPtr;
}

inline void CoSimIO::RegisterExportGeometry(void (*pFuncPtr)(const std::string&))
{
    mpExportGeom = pFuncPtr;
}

inline void CoSimIO::RegisterImportMesh(void (*pFuncPtr)(const std::string&))
{
    mpImportMesh = pFuncPtr;
}

inline void CoSimIO::RegisterExportMesh(void (*pFuncPtr)(const std::string&))
{
    mpExportMesh = pFuncPtr;
}

inline void CoSimIO::RegisterImportData(void (*pFuncPtr)(const std::string&))
{
    mpImportData = pFuncPtr;
}

inline void CoSimIO::RegisterExportData(void (*pFuncPtr)(const std::string&))
{
    mpExportData = pFuncPtr;
}


inline void CoSimIO::Run()
{
    CoSim::Internals::ControlSignal control_signal;
    std::string identifier;
    while(true) {
        control_signal = RecvControlSignal(identifier);
        if (control_signal == CoSim::Internals::ControlSignal::BreakSolutionLoop) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Received control-signal for: BreakSolutionLoop" << std::endl;
            break; // coupled simulation is done
        } else if (control_signal == CoSim::Internals::ControlSignal::AdvanceInTime) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Received control-signal for: AdvanceInTime" << std::endl;
            KRATOS_CO_SIM_ERROR_IF_NOT(mpAdvInTime) << "No function was registered for \"AdvanceInTime\"!" << std::endl;

            std::vector<double> time_vec(1);
            DataContainers::Data time_data = {time_vec};
            Import(time_data, "time_from_co_sim");
            time_vec[0] = mpAdvInTime(time_vec[0]);
            Export(time_data, "time_to_co_sim");
        } else if (control_signal == CoSim::Internals::ControlSignal::InitializeSolutionStep) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Received control-signal for: InitializeSolutionStep" << std::endl;
            KRATOS_CO_SIM_ERROR_IF_NOT(mpInitSolStep) << "No function was registered for \"InitializeSolutionStep\"!" << std::endl;
            mpInitSolStep();

        } else if (control_signal == CoSim::Internals::ControlSignal::SolveSolutionStep) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Received control-signal for: SolveSolutionStep" << std::endl;
            KRATOS_CO_SIM_ERROR_IF_NOT(mpSolSolStep) << "No function was registered for \"SolveSolutionStep\"!" << std::endl;
            mpSolSolStep();
        } else if (control_signal == CoSim::Internals::ControlSignal::FinalizeSolutionStep) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Received control-signal for: FinalizeSolutionStep" << std::endl;
            KRATOS_CO_SIM_ERROR_IF_NOT(mpFinSolStep) << "No function was registered for \"FinalizeSolutionStep\"!" << std::endl;
            mpFinSolStep();
        } else if (control_signal == CoSim::Internals::ControlSignal::ImportGeometry) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Received control-signal for: ImportGeometry" << std::endl;
            KRATOS_CO_SIM_ERROR_IF_NOT(mpImportGeom) << "No function was registered for \"ImportGeometry\"!" << std::endl;
            mpImportGeom(identifier);
        } else if (control_signal == CoSim::Internals::ControlSignal::ExportGeometry) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Received control-signal for: ExportGeometry" << std::endl;
            KRATOS_CO_SIM_ERROR_IF_NOT(mpExportGeom) << "No function was registered for \"ExportGeometry\"!" << std::endl;
            mpExportGeom(identifier);
        } else if (control_signal == CoSim::Internals::ControlSignal::ImportMesh) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Received control-signal for: ImportMesh" << std::endl;
            KRATOS_CO_SIM_ERROR_IF_NOT(mpImportMesh) << "No function was registered for \"ImportMesh\"!" << std::endl;
            mpImportMesh(identifier);
        } else if (control_signal == CoSim::Internals::ControlSignal::ExportMesh) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Received control-signal for: ExportMesh" << std::endl;
            KRATOS_CO_SIM_ERROR_IF_NOT(mpExportMesh) << "No function was registered for \"ExportMesh\"!" << std::endl;
            mpExportMesh(identifier);
        } else if (control_signal == CoSim::Internals::ControlSignal::ImportData) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Received control-signal for: ImportData" << std::endl;
            KRATOS_CO_SIM_ERROR_IF_NOT(mpImportData) << "No function was registered for \"ImportData\"!" << std::endl;
            mpImportData(identifier);
        } else if (control_signal == CoSim::Internals::ControlSignal::ExportData) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Received control-signal for: ExportData" << std::endl;
            KRATOS_CO_SIM_ERROR_IF_NOT(mpExportData) << "No function was registered for \"ExportData\"!" << std::endl;
            mpExportData(identifier);
        } else {
            KRATOS_CO_SIM_ERROR << "Unknown control signal received: " << static_cast<int>(control_signal) << std::endl;;
        }
    }
}

inline bool CoSimIO::IsConverged()
{
    std::string dummy("");
    return RecvControlSignal(dummy) == CoSim::Internals::ControlSignal::ConvergenceAchieved;
}

template<class DataContainer>
bool CoSimIO::Import(DataContainer& rContainer, const std::string& rIdentifier)
{
    return mpComm->Import(rContainer, rIdentifier);
}

template<class DataContainer>
bool CoSimIO::Export(const DataContainer& rContainer, const std::string& rIdentifier)
{
    return mpComm->Export(rContainer, rIdentifier);
}

inline void CoSimIO::Initialize(const std::string& rName, SettingsType& rSettings, const bool IsConnectionMaster)
{
    std::string comm_format("file"); // default is file-communication
    if (rSettings.count("communication_format") != 0) { // communication format has been specified
        comm_format = rSettings.at("communication_format");
    }

    KRATOS_CO_SIM_INFO("CoSimIO") << "CoSimIO for \"" << rName << "\" uses communication format: " << comm_format << std::endl;

    if (comm_format == "file") {
        mpComm = std::unique_ptr<CoSimComm>(new FileComm(rName, rSettings, IsConnectionMaster)); // make_unique is C++14
    } else if (comm_format == "sockets") {
#ifdef KRATOS_CO_SIM_IO_ENABLE_SOCKETS
        mpComm = std::unique_ptr<CoSimComm>(new SocketsComm(rName, rSettings, IsConnectionMaster)); // make_unique is C++14
#else
        KRATOS_CO_SIM_ERROR << "Support for Sockets was not compiled!" << std::endl;
#endif /* KRATOS_CO_SIM_IO_ENABLE_SOCKETS */
    } else if (comm_format == "mpi") {
#ifdef KRATOS_CO_SIM_IO_ENABLE_MPI
        mpComm = std::unique_ptr<CoSimComm>(new MPIComm(rName, rSettings, IsConnectionMaster)); // make_unique is C++14
#else
        KRATOS_CO_SIM_ERROR << "Support for MPI was not compiled!" << std::endl;
#endif /* KRATOS_CO_SIM_IO_ENABLE_MPI */
    } else {
        KRATOS_CO_SIM_ERROR << "Unsupported communication format: " << comm_format << std::endl;
    }
}

} // namespace CoSim

#endif /* KRATOS_CO_SIM_IO_IMPL_H_INCLUDED */
