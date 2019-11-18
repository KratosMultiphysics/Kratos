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
#include <map>
#include <stdexcept>

// Project includes
#include "co_sim_file_comm.h"

namespace CoSimIO {
namespace Internals {

class CoSimIOImpl
{

public:
    typedef CoSimComm::SettingsType SettingsType;

    typedef void (*DataExchangeFunctionType)(const std::string&);

    explicit CoSimIOImpl(const std::string& rName, const std::string& rSettingsFileName, const bool IsConnectionMaster=false)
    : CoSimIOImpl(rName, Internals::ReadSettingsFile(rSettingsFileName), IsConnectionMaster) { } // forwarding constructor call

    explicit CoSimIOImpl(const std::string& rName, SettingsType rSettings, const bool IsConnectionMaster=false) : mIsConnectionMaster(IsConnectionMaster)
    {
        Initialize(rName, rSettings, IsConnectionMaster);
    }

    bool Connect()
    {
        return mpComm->Connect();
    }

    bool Disconnect()
    {
        return mpComm->Connect();
    }

    void SendControlSignal(const Internals::ControlSignal Signal, const std::string& rIdentifier)
    {
        KRATOS_CO_SIM_ERROR_IF_NOT(mIsConnectionMaster) << "This function can only be called as the Connection-Master!" << std::endl;
        mpComm->SendControlSignal(Signal, rIdentifier);
    }
    Internals::ControlSignal RecvControlSignal(std::string& rIdentifier)
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
        return mpComm->RecvControlSignal(rIdentifier);
    }


    void RegisterDataExchange(DataExchangeFunctionType pFuncPtr, const std::string& rName)
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
        KRATOS_CO_SIM_ERROR_IF(mDataExchangeFunctions.count(rName) > 0) << "Function already registered for \"" << rName << "\"!" << std::endl;
        mDataExchangeFunctions[rName] = pFuncPtr;
    }

    void RegisterAdvanceInTime(double (*pFuncPtr)(double))
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
        mpAdvInTime = pFuncPtr;
    }

    void RegisterInitializeSolutionStep(void (*pFuncPtr)())
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
        mpInitSolStep = pFuncPtr;
    }

    void RegisterSolveSolutionStep(void (*pFuncPtr)())
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
        mpSolSolStep = pFuncPtr;
    }

    void RegisterFinalizeSolutionStep(void (*pFuncPtr)())
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;
        mpFinSolStep = pFuncPtr;
    }

    void Run()
    {
        KRATOS_CO_SIM_ERROR_IF(mIsConnectionMaster) << "This function can only be called as the Connection-Slave!" << std::endl;

        const std::map<const CoSimIO::Internals::ControlSignal, const std::string> signal_to_name = {
            {CoSimIO::Internals::ControlSignal::ImportGeometry, "ImportGeometry"},
            {CoSimIO::Internals::ControlSignal::ExportGeometry, "ExportGeometry"},
            {CoSimIO::Internals::ControlSignal::ImportMesh,     "ImportMesh"},
            {CoSimIO::Internals::ControlSignal::ExportMesh,     "ExportMesh"},
            {CoSimIO::Internals::ControlSignal::ImportData,     "ImportData"},
            {CoSimIO::Internals::ControlSignal::ExportData,     "ExportData"}
        };

        CoSimIO::Internals::ControlSignal control_signal;
        std::string identifier;
        while(true) {
            control_signal = RecvControlSignal(identifier);
            if (control_signal == CoSimIO::Internals::ControlSignal::BreakSolutionLoop) {
                break; // coupled simulation is done
            } else if (control_signal == CoSimIO::Internals::ControlSignal::AdvanceInTime) {
                KRATOS_CO_SIM_ERROR_IF_NOT(mpAdvInTime) << "No function was registered for \"AdvanceInTime\"!" << std::endl;

                std::vector<double> time_vec(1);
                // DataContainers::Data time_data = {time_vec};
                // Import(time_data, "time_from_co_sim");
                time_vec[0] = mpAdvInTime(time_vec[0]);
                // Export(time_data, "time_to_co_sim");
            } else if (control_signal == CoSimIO::Internals::ControlSignal::InitializeSolutionStep) {
                KRATOS_CO_SIM_ERROR_IF_NOT(mpInitSolStep) << "No function was registered for \"InitializeSolutionStep\"!" << std::endl;
                mpInitSolStep();

            } else if (control_signal == CoSimIO::Internals::ControlSignal::SolveSolutionStep) {
                KRATOS_CO_SIM_ERROR_IF_NOT(mpSolSolStep) << "No function was registered for \"SolveSolutionStep\"!" << std::endl;
                mpSolSolStep();
            } else if (control_signal == CoSimIO::Internals::ControlSignal::FinalizeSolutionStep) {
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

    bool IsConverged()
    {
        std::string dummy("");
        return RecvControlSignal(dummy) == CoSimIO::Internals::ControlSignal::ConvergenceAchieved;
    }

    void ImportData(
        const char* pIdentifier,
        int* pSize,
        double** ppData)
    {
        mpComm->ImportData(pIdentifier, pSize, ppData);
    }

    void ExportData(
        const char* pIdentifier,
        const int Size,
        const double* pData)
    {
        mpComm->ExportData(pIdentifier, Size, pData);
    }

    void ImportMesh(
        const char* pIdentifier,
        int* pNumberOfNodes,
        int* pNumberOfElements,
        double** ppNodalCoordinates,
        int** ppElementConnectivities,
        int** ppElementTypes)
    {
        mpComm->ImportMesh(pIdentifier, pNumberOfNodes, pNumberOfElements, ppNodalCoordinates, ppElementConnectivities, ppElementTypes);
    }

    void ExportMesh(
        const char* pIdentifier,
        const int NumberOfNodes,
        const int NumberOfElements,
        const double* pNodalCoordinates,
        const int* pElementConnectivities,
        const int* pElementTypes)
    {
        mpComm->ExportMesh(pIdentifier, NumberOfNodes, NumberOfElements, pNodalCoordinates, pElementConnectivities, pElementTypes);
    }

    void ImportGeometry()
    {
        KRATOS_CO_SIM_ERROR << "Importing of Geometry is not yet implemented!" << std::endl;
    }

    void ExportGeometry()
    {
        KRATOS_CO_SIM_ERROR << "Exporting of Geometry is not yet implemented!" << std::endl;
    }

private:
    std::unique_ptr<CoSimComm> mpComm; // handles communication (File, Sockets, MPI, ...)

    bool mIsConnectionMaster = false;

    double (*mpAdvInTime)(double) = nullptr;
    void (*mpInitSolStep)()     = nullptr;
    void (*mpSolSolStep)()      = nullptr;
    void (*mpFinSolStep)()      = nullptr;

    void (*mpImportData)()      = nullptr;
    void (*mpExportData)()      = nullptr;
    void (*mpImportMesh)()      = nullptr;
    void (*mpExportMesh)()      = nullptr;
    void (*mpImportGeometry)()  = nullptr;
    void (*mpExportGeometry)()  = nullptr;

    std::map<const std::string, DataExchangeFunctionType> mDataExchangeFunctions;

    void Initialize(const std::string& rName, SettingsType& rSettings, const bool IsConnectionMaster)
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

}; // class CoSimIOImpl


// TODO make sure this is unique even across compilation units (test somehow)
static std::unordered_map<std::string, std::unique_ptr<CoSimIOImpl>> s_co_sim_ios;

static bool HasIO(const char* pConnectionName)
{
    return s_co_sim_ios.find(std::string(pConnectionName)) != s_co_sim_ios.end();
}

static CoSimIOImpl& GetIO(const char* pConnectionName)
{
    KRATOS_CO_SIM_ERROR_IF_NOT(HasIO(pConnectionName)) << "Trying to use connection " << pConnectionName << " which does not exist!" << std::endl;
    return *s_co_sim_ios.at(std::string(pConnectionName));
}

} // namespace Internals
} // namespace CoSimIO

#endif /* KRATOS_CO_SIM_IO_IMPL_H_INCLUDED */
