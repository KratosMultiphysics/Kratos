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
#include <unordered_map>
#include <string>
#include <stdexcept>

// Project includes
#include "co_sim_data_containers.h"
#include "co_sim_api_internals.h"

namespace CoSim {

#define CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE(DataContainerType)                              \
    bool Import(DataContainerType& rDataContainer, const std::string& rIdentifier)               \
        { CheckConnection(); return ImportDetail(rDataContainer, rIdentifier); }                 \
    bool Export(const DataContainerType& rDataContainer, const std::string& rIdentifier)         \
        { CheckConnection(); return ExportDetail(rDataContainer, rIdentifier); }                 \

#define CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE_DETAIL(DataContainerType)                       \
    virtual bool ImportDetail(DataContainerType& rDataContainer, const std::string& rIdentifier) \
        { throw std::runtime_error("Type of data not yet supported"); }                          \
    virtual bool ExportDetail(const DataContainerType& rDataContainer, const std::string& rIdentifier) \
        { throw std::runtime_error("Type of data not yet supported"); }                          \

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
        CS_LOG << "Connecting \"" << mName << "\" as Connection-" << (mIsConnectionMaster ? "MASTER" : "SLAVE") << " ..." << std::endl;

        if (mIsConnected) {
            throw std::runtime_error("A connection was already established!");
        }

        mIsConnected = ConnectDetail();

        if (!mIsConnected) {
            throw std::runtime_error("Connection was not successful!");
        }

        CS_LOG << "Connection established" << std::endl;

        return mIsConnected;
    }

    bool Disconnect()
    {
        CS_LOG << "Disconnecting \"" << mName << "\" ..." << std::endl;

        if (mIsConnected) {
            mIsConnected = !DisconnectDetail();
            if (mIsConnected) {
                CS_LOG << "Warning: Disconnect was not successful!" << std::endl;
                return false;
            }
        } else {
            CS_LOG << "Warning: Calling Disconnect but there was no active connection!" << std::endl;
            return false;
        }

        CS_LOG << "Disconnecting successful" << std::endl;

        return true;
    }

    void SendControlSignal(const int rSignal, const std::string& rIdentifier)
    {
        CheckConnection(); return SendControlSignalDetail(rSignal, rIdentifier);
    }
    int RecvControlSignal(std::string& rIdentifier)
    {
        CheckConnection(); return RecvControlSignalDetail(rIdentifier);
    }

    CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE(DataContainers::Geometry);
    CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE(DataContainers::Mesh);
    CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE(DataContainers::Data);

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

    virtual void SendControlSignalDetail(const int rSignal, const std::string& rIdentifier)
    {
        throw std::runtime_error("SendControlSignalDetail not implemented for this comm-type");
    }
    virtual int RecvControlSignalDetail(std::string& rIdentifier)
    {
        throw std::runtime_error("RecvControlSignalDetail not implemented for this comm-type");
    }

    CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE_DETAIL(DataContainers::Geometry);
    CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE_DETAIL(DataContainers::Mesh);
    CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE_DETAIL(DataContainers::Data);

    void CheckConnection()
    {
        if (!mIsConnected) {
            throw std::runtime_error("No active connection exists!");
        }
    }

};

#undef CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE // this macro should only be used for this class
#undef CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE_DETAIL // this macro should only be used for this class

} // namespace CoSim

#endif /* KRATOS_CO_SIM_COMM_H_INCLUDED */
