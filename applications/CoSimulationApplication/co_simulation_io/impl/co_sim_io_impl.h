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

inline CoSimIO::CoSimIO(const std::string& rName, SettingsType rSettings, const bool IsConnectionMaster) : mIsConnectionMaster(IsConnectionMaster)
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
    KRATOS_CO_SIM_ERROR_IF_NOT(mIsConnectionMaster) << "This function can only be called as the Connection-Master!" << std::endl;
    mpComm->SendControlSignal(Signal, rIdentifier);
}
inline Internals::ControlSignal CoSimIO::RecvControlSignal(std::string& rIdentifier)
{
    KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
    return mpComm->RecvControlSignal(rIdentifier);
}


inline void CoSimIO::RegisterDataExchange(DataExchangeFunctionType pFuncPtr, const std::string& rName)
{
    KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
    KRATOS_CO_SIM_ERROR_IF(mDataExchangeFunctions.count(rName) > 0) << "Function already registered for \"" << rName << "\"!" << std::endl;
    mDataExchangeFunctions[rName] = pFuncPtr;
}


inline void CoSimIO::RegisterAdvanceInTime(double (*pFuncPtr)(double))
{
    KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
    mpAdvInTime = pFuncPtr;
}

inline void CoSimIO::RegisterInitializeSolutionStep(void (*pFuncPtr)())
{
    KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
    mpInitSolStep = pFuncPtr;
}

inline void CoSimIO::RegisterSolveSolutionStep(void (*pFuncPtr)())
{
    KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
    mpSolSolStep = pFuncPtr;
}

inline void CoSimIO::RegisterFinalizeSolutionStep(void (*pFuncPtr)())
{
    KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
    mpFinSolStep = pFuncPtr;
}

inline void CoSimIO::Run()
{
    KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;

    const std::map<const CoSim::Internals::ControlSignal, const std::string> signal_to_name = {
        {CoSim::Internals::ControlSignal::ImportGeometry, "ImportGeometry"},
        {CoSim::Internals::ControlSignal::ExportGeometry, "ExportGeometry"},
        {CoSim::Internals::ControlSignal::ImportMesh,     "ImportMesh"},
        {CoSim::Internals::ControlSignal::ExportMesh,     "ExportMesh"},
        {CoSim::Internals::ControlSignal::ImportData,     "ImportData"},
        {CoSim::Internals::ControlSignal::ExportData,     "ExportData"}
    };

    CoSim::Internals::ControlSignal control_signal;
    std::string identifier;
    while(true) {
        control_signal = RecvControlSignal(identifier);
        if (control_signal == CoSim::Internals::ControlSignal::BreakSolutionLoop) {
            break; // coupled simulation is done
        } else if (control_signal == CoSim::Internals::ControlSignal::AdvanceInTime) {
            KRATOS_CO_SIM_ERROR_IF_NOT(mpAdvInTime) << "No function was registered for \"AdvanceInTime\"!" << std::endl;

            std::vector<double> time_vec(1);
            DataContainers::Data time_data = {time_vec};
            Import(time_data, "time_from_co_sim");
            time_vec[0] = mpAdvInTime(time_vec[0]);
            Export(time_data, "time_to_co_sim");
        } else if (control_signal == CoSim::Internals::ControlSignal::InitializeSolutionStep) {
            KRATOS_CO_SIM_ERROR_IF_NOT(mpInitSolStep) << "No function was registered for \"InitializeSolutionStep\"!" << std::endl;
            mpInitSolStep();

        } else if (control_signal == CoSim::Internals::ControlSignal::SolveSolutionStep) {
            KRATOS_CO_SIM_ERROR_IF_NOT(mpSolSolStep) << "No function was registered for \"SolveSolutionStep\"!" << std::endl;
            mpSolSolStep();
        } else if (control_signal == CoSim::Internals::ControlSignal::FinalizeSolutionStep) {
            KRATOS_CO_SIM_ERROR_IF_NOT(mpFinSolStep) << "No function was registered for \"FinalizeSolutionStep\"!" << std::endl;
            mpFinSolStep();
        } else if (signal_to_name.count(control_signal) > 0) {
            const auto& r_function_name(signal_to_name.at(control_signal));

            KRATOS_CO_SIM_ERROR_IF_NOT((mDataExchangeFunctions.count(r_function_name)>0)) << "No function was registered for \"" << r_function_name << "\"!" << std::endl;
            mDataExchangeFunctions.at(r_function_name)(identifier);
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
