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

#ifndef KRATOS_CO_SIM_COMM_H_INCLUDED
#define KRATOS_CO_SIM_COMM_H_INCLUDED

// System includes

// Project includes
#include "co_sim_io_internals.h"

namespace CoSimIO {

class CoSimComm
{
public:
    typedef Internals::SettingsType SettingsType;

    explicit CoSimComm(const std::string& rName, SettingsType& rSettings, const bool IsConnectionMaster) : mrSettings(rSettings),mName(rName), mIsConnectionMaster(IsConnectionMaster)
    {
        const SettingsType default_settings = {
            {"echo_level",   "1"},
            {"print_timing", "0"}
        };
        Internals::AddMissingSettings(default_settings, mrSettings);

        mEchoLevel = std::stoi(mrSettings.at("echo_level"));
        mPrintTiming = (mrSettings.at("print_timing") == "1");
    }

    virtual ~CoSimComm() = default; // impl of disconnect has to be in derived class due to order of class destruction

    bool Connect()
    {
        KRATOS_CO_SIM_INFO("CoSimIO") << "Connecting \"" << mName << "\" as Connection-" << (mIsConnectionMaster ? "MASTER" : "SLAVE") << " ..." << std::endl;

        KRATOS_CO_SIM_ERROR_IF(mIsConnected) << "A connection was already established!" << std::endl;

        mIsConnected = ConnectDetail();

        KRATOS_CO_SIM_ERROR_IF_NOT(mIsConnected) << "Connection was not successful!" << std::endl;

        KRATOS_CO_SIM_INFO("CoSimIO") << "Connection established" << std::endl;

        return mIsConnected;
    }

    bool Disconnect()
    {
        KRATOS_CO_SIM_INFO("CoSimIO") << "Disconnecting \"" << mName << "\" ..." << std::endl;

        if (mIsConnected) {
            mIsConnected = !DisconnectDetail();
            if (mIsConnected) {
                KRATOS_CO_SIM_INFO("CoSimIO") << "Warning: Disconnect was not successful!" << std::endl;
                return false;
            }
        } else {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Warning: Calling Disconnect but there was no active connection!" << std::endl;
            return false;
        }

        KRATOS_CO_SIM_INFO("CoSimIO") << "Disconnecting successful" << std::endl;

        return true;
    }

    void SendControlSignal(Internals::ControlSignal Signal, const std::string& rIdentifier)
    {
        CheckConnection(); return SendControlSignalDetail(Signal, rIdentifier);
    }
    Internals::ControlSignal RecvControlSignal(std::string& rIdentifier)
    {
        CheckConnection(); return RecvControlSignalDetail(rIdentifier);
    }

    void ImportData(
        const std::string& rIdentifier,
        int* pSize,
        double** ppData)
    {
        CheckConnection();
        ImportDataImpl(rIdentifier, pSize, ppData);
    }

    void ExportData(
        const std::string& rIdentifier,
        const int Size,
        const double* pData)
    {
        CheckConnection();
        ExportDataImpl(rIdentifier, Size, pData);
    }

    void ImportMesh(
        const std::string& rIdentifier,
        int* pNumberOfNodes,
        int* pNumberOfElements,
        double** ppNodalCoordinates,
        int** ppElementConnectivities,
        int** ppElementTypes)
    {
        CheckConnection();
        ImportMeshImpl(rIdentifier, pNumberOfNodes, pNumberOfElements, ppNodalCoordinates, ppElementConnectivities, ppElementTypes);
    }

    void ExportMesh(
        const std::string& rIdentifier,
        const int NumberOfNodes,
        const int NumberOfElements,
        const double* pNodalCoordinates,
        const int* pElementConnectivities,
        const int* pElementTypes)
    {
        CheckConnection();
        ExportMeshImpl(rIdentifier, NumberOfNodes, NumberOfElements, pNodalCoordinates, pElementConnectivities, pElementTypes);
    }

protected:
    SettingsType& mrSettings;

    std::string GetName() const        {return mName;}
    int GetEchoLevel() const           {return mEchoLevel;}
    bool GetIsConnectionMaster() const {return mIsConnectionMaster;}
    bool GetPrintTiming() const        {return mPrintTiming;}
    bool GetIsConnected() const        {return mIsConnected;}

private:
    std::string mName;
    int mEchoLevel = 1;
    bool mIsConnectionMaster = false;
    bool mPrintTiming = false;
    bool mIsConnected = false;

    virtual bool ConnectDetail() = 0;
    virtual bool DisconnectDetail() = 0;

    virtual void SendControlSignalDetail(Internals::ControlSignal Signal, const std::string& rIdentifier)
    {
        KRATOS_CO_SIM_ERROR << "SendControlSignalDetail not implemented for this comm-type" << std::endl;
    }
    virtual Internals::ControlSignal RecvControlSignalDetail(std::string& rIdentifier)
    {
        KRATOS_CO_SIM_ERROR << "RecvControlSignalDetail not implemented for this comm-type" << std::endl;
        return Internals::ControlSignal::Dummy;
    }

    virtual void ImportDataImpl(
        const std::string& rIdentifier,
        int* pSize,
        double** ppData)
    {
        KRATOS_CO_SIM_ERROR << "ImportDataImpl not implemented for this comm-type!" << std::endl;
    }

    virtual void ExportDataImpl(
        const std::string& rIdentifier,
        const int Size,
        const double* pData)
    {
        KRATOS_CO_SIM_ERROR << "ImportDataImpl not implemented for this comm-type!" << std::endl;
    }

    virtual void ImportMeshImpl(
        const std::string& rIdentifier,
        int* pNumberOfNodes,
        int* pNumberOfElements,
        double** ppNodalCoordinates,
        int** ppElementConnectivities,
        int** ppElementTypes)
    {
        KRATOS_CO_SIM_ERROR << "ImportDataImpl not implemented for this comm-type!" << std::endl;
    }

    virtual void ExportMeshImpl(
        const std::string& rIdentifier,
        const int NumberOfNodes,
        const int NumberOfElements,
        const double* pNodalCoordinates,
        const int* pElementConnectivities,
        const int* pElementTypes)
    {
        KRATOS_CO_SIM_ERROR << "ImportDataImpl not implemented for this comm-type!" << std::endl;
    }

    void CheckConnection()
    {
        KRATOS_CO_SIM_ERROR_IF_NOT(mIsConnected) << "No active connection exists!" << std::endl;;
    }
};

} // namespace CoSimIO

#endif /* KRATOS_CO_SIM_COMM_H_INCLUDED */
