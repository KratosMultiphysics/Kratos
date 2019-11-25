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

#ifndef KRATOS_CO_SIM_CONNECTION_H_INCLUDED
#define KRATOS_CO_SIM_CONNECTION_H_INCLUDED

// Optional includes
#ifdef KRATOS_CO_SIM_IO_ENABLE_SOCKETS
#include "co_sim_sockets_communication.h"
#endif /* KRATOS_CO_SIM_IO_ENABLE_SOCKETS */


#ifdef KRATOS_CO_SIM_IO_ENABLE_MPI
#include "co_sim_mpi_communication.h"
#endif /* KRATOS_CO_SIM_IO_ENABLE_MPI */

// System includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <utility>

// Project includes
#include "co_sim_file_communication.h"

namespace CoSimIO {
namespace Internals {

class CoSimConnection
{

public:
    typedef void (*DataExchangeFunctionType)(const std::string&, const std::string&);
    typedef void (*DataExchangeCFunctionType)(const char*, const char*);

    explicit CoSimConnection(const std::string& rName, const std::string& rSettingsFileName)
    : CoSimConnection(rName, Internals::ReadSettingsFile(rSettingsFileName)) { } // forwarding constructor call

    explicit CoSimConnection(const std::string& rName, SettingsType Settings) : mConnectionName(rName)
    {
        Initialize(Settings);
    }

    bool Connect()
    {
        return mpComm->Connect();
    }

    bool Disconnect()
    {
        return mpComm->Disconnect();
    }

    void SendControlSignal(const std::string& rIdentifier, const CoSimIO::ControlSignal Signal)
    {
        KRATOS_CO_SIM_ERROR_IF_NOT(mIsConnectionMaster) << "This function can only be called as the Connection-Master!" << std::endl;
        mpComm->SendControlSignal(rIdentifier, Signal);
    }

    // Only used for AdvanceInTime
    void Register(
        const std::string& rFunctionName,
        void (*pFunctionPointer)(double*))
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
        KRATOS_CO_SIM_ERROR_IF(rFunctionName != "AdvanceInTime") << "Only \"AdvanceInTime\" can be registered with this function, trying to register \"" << rFunctionName << "\"!" << std::endl;
        KRATOS_CO_SIM_ERROR_IF(mpAdvInTime) << "A function was already registered for " << rFunctionName << "!" << std::endl;

        mpAdvInTime = pFunctionPointer;
    }

    // Used for the solving functions
    void Register(
        const std::string& rFunctionName,
        void (*pFunctionPointer)(void))
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;

        if (rFunctionName == "InitializeSolutionStep") {
            KRATOS_CO_SIM_ERROR_IF(mpInitSolStep) << "A function was already registered for \"InitializeSolutionStep\"!" << std::endl;
            mpInitSolStep = pFunctionPointer;
        } else if (rFunctionName == "SolveSolutionStep") {
            KRATOS_CO_SIM_ERROR_IF(mpSolSolStep) << "A function was already registered for \"SolveSolutionStep\"!" << std::endl;
            mpSolSolStep = pFunctionPointer;
        } else if (rFunctionName == "FinalizeSolutionStep") {
            KRATOS_CO_SIM_ERROR_IF(mpFinSolStep) << "A function was already registered for \"FinalizeSolutionStep\"!" << std::endl;
            mpFinSolStep = pFunctionPointer;
        } else {
            KRATOS_CO_SIM_ERROR << "Only functions for \"InitializeSolutionStep\", \"SolveSolutionStep\" \"FinalizeSolutionStep\" can be registered using this function, tried registering " << rFunctionName << std::endl;
        }
    }

    // Used for the data-exchange functions coming from C or Fortran
    void Register(
        const std::string& rFunctionName,
        void (*pFunctionPointer)(const char*, const char*))
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
        KRATOS_CO_SIM_ERROR_IF(mDataExchangeFunctions.size() > 0) << "Mixing of registering functions with different arguments is not allowed!" << std::endl;
        KRATOS_CO_SIM_ERROR_IF((mDataExchangeCFunctions.count(rFunctionName)>0)) << "A function was already registered for " << rFunctionName << "!" << std::endl;
        // TODO maybe check if the name of the function is allowed?
        mDataExchangeCFunctions[rFunctionName] = pFunctionPointer;
    }

    // Used for the data-exchange functions coming from C++
    void Register(
        const std::string& rFunctionName,
        void (*pFunctionPointer)(const std::string&, const std::string&))
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
        KRATOS_CO_SIM_ERROR_IF(mDataExchangeCFunctions.size() > 0) << "Mixing of registering functions with different arguments is not allowed!" << std::endl;
        KRATOS_CO_SIM_ERROR_IF((mDataExchangeFunctions.count(rFunctionName)>0)) << "A function was already registered for " << rFunctionName << "!" << std::endl;
        // TODO maybe check if the name of the function is allowed?
        mDataExchangeFunctions[rFunctionName] = pFunctionPointer;
    }

    void Run()
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;

        const std::map<const CoSimIO::ControlSignal, const std::string> signal_to_name {
            {CoSimIO::ControlSignal::ImportGeometry, "ImportGeometry"},
            {CoSimIO::ControlSignal::ExportGeometry, "ExportGeometry"},
            {CoSimIO::ControlSignal::ImportMesh,     "ImportMesh"},
            {CoSimIO::ControlSignal::ExportMesh,     "ExportMesh"},
            {CoSimIO::ControlSignal::ImportData,     "ImportData"},
            {CoSimIO::ControlSignal::ExportData,     "ExportData"}
        };

        CoSimIO::ControlSignal control_signal;
        std::string identifier;
        while(true) {
            control_signal = RecvControlSignal(identifier);
            if (control_signal == CoSimIO::ControlSignal::BreakSolutionLoop) {
                break; // coupled simulation is done
            } else if (control_signal == CoSimIO::ControlSignal::AdvanceInTime) {
                KRATOS_CO_SIM_ERROR_IF_NOT(mpAdvInTime) << "No function was registered for \"AdvanceInTime\"!" << std::endl;

                std::vector<double> time_vec(1);
                // DataContainers::Data time_data = {time_vec};
                // Import(time_data, "time_from_co_sim");
                mpAdvInTime(&time_vec[0]);
                // Export(time_data, "time_to_co_sim");
            } else if (control_signal == CoSimIO::ControlSignal::InitializeSolutionStep) {
                KRATOS_CO_SIM_ERROR_IF_NOT(mpInitSolStep) << "No function was registered for \"InitializeSolutionStep\"!" << std::endl;
                mpInitSolStep();

            } else if (control_signal == CoSimIO::ControlSignal::SolveSolutionStep) {
                KRATOS_CO_SIM_ERROR_IF_NOT(mpSolSolStep) << "No function was registered for \"SolveSolutionStep\"!" << std::endl;
                mpSolSolStep();
            } else if (control_signal == CoSimIO::ControlSignal::FinalizeSolutionStep) {
                KRATOS_CO_SIM_ERROR_IF_NOT(mpFinSolStep) << "No function was registered for \"FinalizeSolutionStep\"!" << std::endl;
                mpFinSolStep();
            } else if (signal_to_name.count(control_signal) > 0) {
                const auto& r_function_name(signal_to_name.at(control_signal));

                KRATOS_CO_SIM_ERROR_IF_NOT((mDataExchangeFunctions.count(r_function_name)>0)) << "No function was registered for \"" << r_function_name << "\"!" << std::endl;
                mDataExchangeFunctions.at(r_function_name)(mConnectionName, identifier);
            } else {
                KRATOS_CO_SIM_ERROR << "Unknown control signal received: " << static_cast<int>(control_signal) << std::endl;;
            }
        }
    }

    void IsConverged(int& rConvergenceSignal)
    {
        std::string dummy("");
        rConvergenceSignal = (RecvControlSignal(dummy) == CoSimIO::ControlSignal::ConvergenceAchieved);
    }


    template<class... Args>
    void ImportData(Args&&... args)
    {
        mpComm->ImportData(std::forward<Args>(args)...);
    }

    template<class... Args>
    void ExportData(Args&&... args)
    {
        mpComm->ExportData(std::forward<Args>(args)...);
    }

    template<class... Args>
    void ImportMesh(Args&&... args)
    {
        mpComm->ImportMesh(std::forward<Args>(args)...);
    }

    template<class... Args>
    void ExportMesh(Args&&... args)
    {
        mpComm->ExportMesh(std::forward<Args>(args)...);
    }

    template<class... Args>
    void ImportGeometry(Args&&... args)
    {
        KRATOS_CO_SIM_ERROR << "Importing of Geometry is not yet implemented!" << std::endl;
        mpComm->ImportGeometry(std::forward<Args>(args)...);
    }

    template<class... Args>
    void ExportGeometry(Args&&... args)
    {
        KRATOS_CO_SIM_ERROR << "Exporting of Geometry is not yet implemented!" << std::endl;
        mpComm->ExportGeometry(std::forward<Args>(args)...);
    }

private:
    std::unique_ptr<CoSimCommunication> mpComm; // handles communication (File, Sockets, MPI, ...)

    std::string mConnectionName;

    bool mIsConnectionMaster = false;

    void (*mpAdvInTime)(double*) = nullptr;
    void (*mpInitSolStep)()      = nullptr;
    void (*mpSolSolStep)()       = nullptr;
    void (*mpFinSolStep)()       = nullptr;

    std::map<const std::string, DataExchangeFunctionType> mDataExchangeFunctions;
    std::map<const std::string, DataExchangeCFunctionType> mDataExchangeCFunctions;

    void Initialize(SettingsType& rSettings)
    {
        std::string comm_format("file"); // default is file-communication
        if (rSettings.count("communication_format") != 0) { // communication format has been specified
            comm_format = rSettings.at("communication_format");
        }

        if (rSettings.count("is_connection_master") != 0) { // is_connection_master has been specified
            mIsConnectionMaster = (rSettings.at("is_connection_master") == "1");
        }

        KRATOS_CO_SIM_INFO("CoSimIO") << "CoSimIO for \"" << mConnectionName << "\" uses communication format: " << comm_format << std::endl;

        if (comm_format == "file") {
            mpComm = std::unique_ptr<CoSimCommunication>(new CoSimFileCommunication(mConnectionName, rSettings, mIsConnectionMaster));
        } else if (comm_format == "sockets") {
            #ifdef KRATOS_CO_SIM_IO_ENABLE_SOCKETS
            mpComm = std::unique_ptr<CoSimCommunication>(new CoSimSocketsCommunication(mConnectionName, rSettings, mIsConnectionMaster));
            #else
            KRATOS_CO_SIM_ERROR << "Support for Sockets was not compiled!" << std::endl;
            #endif /* KRATOS_CO_SIM_IO_ENABLE_SOCKETS */
        } else if (comm_format == "mpi") {
            #ifdef KRATOS_CO_SIM_IO_ENABLE_MPI
            mpComm = std::unique_ptr<CoSimCommunication>(new CoSimMPICommunication(mConnectionName, rSettings, mIsConnectionMaster));
            #else
            KRATOS_CO_SIM_ERROR << "Support for MPI was not compiled!" << std::endl;
            #endif /* KRATOS_CO_SIM_IO_ENABLE_MPI */
        } else {
            KRATOS_CO_SIM_ERROR << "Unsupported communication format: " << comm_format << std::endl;
        }
    }

    CoSimIO::ControlSignal RecvControlSignal(std::string& rIdentifier)
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
        return mpComm->RecvControlSignal(rIdentifier);
    }

}; // class CoSimConnection

} // namespace Internals

} // namespace CoSimIO

#endif /* KRATOS_CO_SIM_CONNECTION_H_INCLUDED */
