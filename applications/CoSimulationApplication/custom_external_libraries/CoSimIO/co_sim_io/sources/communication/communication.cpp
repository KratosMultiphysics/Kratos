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

// System includes
#include <thread>
#include <system_error>

// Project includes
#include "includes/communication/communication.hpp"
#include "includes/file_serializer.hpp"
#include "includes/utilities.hpp"
#include "includes/version.hpp"

namespace CoSimIO {
namespace Internals {

void AddFilePermissions(const fs::path& rPath)
{
    CO_SIM_IO_TRY

    fs::permissions(rPath,
        fs::perms::owner_read  | fs::perms::owner_write |
        fs::perms::group_read  | fs::perms::group_write |
        fs::perms::others_read | fs::perms::others_write,
        fs::perm_options::add);

    CO_SIM_IO_CATCH
}

Communication::Communication(
    const Info& I_Settings,
    std::shared_ptr<DataCommunicator> I_DataComm)
    : mpDataComm(I_DataComm),
      mMyName(I_Settings.Get<std::string>("my_name")),
      mConnectTo(I_Settings.Get<std::string>("connect_to")),
      mAlwaysUseSerializer(I_Settings.Get<bool>("always_use_serializer", false)),
      mWorkingDirectory(I_Settings.Get<std::string>("working_directory", fs::relative(fs::current_path()).string())),
      mEchoLevel(I_Settings.Get<int>("echo_level", 0)),
      mPrintTiming(I_Settings.Get<bool>("print_timing", false))
{
    CO_SIM_IO_TRY

    if (I_Settings.Has("is_primary_connection")) {
        mIsPrimaryConnection = I_Settings.Get<bool>("is_primary_connection");
        mPrimaryWasExplicitlySpecified = true;
    } else {
        // automatically determine the primary connection in case the user didn't specify it
        mIsPrimaryConnection = mMyName < mConnectTo;
        mPrimaryWasExplicitlySpecified = false;
    }
    mConnectionName = Utilities::CreateConnectionName(mMyName, mConnectTo);

    CO_SIM_IO_ERROR_IF_NOT(fs::exists(mWorkingDirectory)) << "The working directory " << mWorkingDirectory << " does not exist!" << std::endl;

    if (I_Settings.Has("serializer_trace_type")) {
        mSerializerTraceType = Serializer::StringToTraceType(I_Settings.Get<std::string>("serializer_trace_type"));
    }

    mCommInFolder = I_Settings.Get<bool>("use_folder_for_communication", true);
    mCommFolder = GetWorkingDirectory();
    if (mCommInFolder) {
        mCommFolder /= ".CoSimIOComm_" + GetConnectionName();
    }

    CO_SIM_IO_CATCH
}

Info Communication::Connect(const Info& I_Info)
{
    CO_SIM_IO_TRY

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>0 && mpDataComm->Rank() == 0)
        << "Establishing connection for \"" << mConnectionName
        << "\"\n    from: \"" << mMyName
        << "\"\n    to:   \"" << mConnectTo
        << "\"\n    as " << (mIsPrimaryConnection ? "PRIMARY" : "SECONDARY")
        << " connection; working directory: " << mWorkingDirectory << " ..." << std::endl;

    CO_SIM_IO_ERROR_IF(mIsConnected) << "A connection was already established!" << std::endl;

    BaseConnectDetail(I_Info);

    PrepareConnection(I_Info);

    HandShake(I_Info);

    Info connect_detail_info = ConnectDetail(I_Info);
    mIsConnected = true;
    connect_detail_info.Set<bool>("is_connected", true);
    connect_detail_info.Set<int>("connection_status", ConnectionStatus::Connected);
    connect_detail_info.Set<std::string>("working_directory", mWorkingDirectory.string());

    CO_SIM_IO_ERROR_IF_NOT(mIsConnected) << "Connection was not successful!" << std::endl;

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>0 && mpDataComm->Rank() == 0) << "Connection established" << std::endl;

    return connect_detail_info;

    CO_SIM_IO_CATCH
}

Info Communication::Disconnect(const Info& I_Info)
{
    CO_SIM_IO_TRY

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>0 && mpDataComm->Rank() == 0) << "Disconnecting \"" << mConnectionName << "\" ..." << std::endl;

    if (mIsConnected) {
        Info disconnect_detail_info = DisconnectDetail(I_Info);
        mIsConnected = false;
        disconnect_detail_info.Set<bool>("is_connected", false);

        BaseDisconnectDetail(I_Info);

        if (mIsConnected) {
            CO_SIM_IO_INFO("CoSimIO") << "Warning: Disconnect was not successful!" << std::endl;
            disconnect_detail_info.Set<int>("connection_status", ConnectionStatus::DisconnectionError);
        } else {
            CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>0 && mpDataComm->Rank() == 0) << "Disconnecting successful" << std::endl;
            disconnect_detail_info.Set<int>("connection_status", ConnectionStatus::Disconnected);
        }
        return disconnect_detail_info;

    } else {
        CO_SIM_IO_INFO("CoSimIO") << "Warning: Calling Disconnect but there was no active connection!" << std::endl;
        Info disconnect_info;
        disconnect_info.Set<bool>("is_connected", false);
        disconnect_info.Set<int>("connection_status", ConnectionStatus::DisconnectionError);
        return disconnect_info;
    }

    CO_SIM_IO_CATCH
}

void Communication::BaseConnectDetail(const Info& I_Info)
{
    CO_SIM_IO_TRY

    if (mCommInFolder && GetIsPrimaryConnection() && mpDataComm->Rank() == 0) {
        // delete and recreate directory to remove potential leftovers
        std::error_code ec;
        fs::remove_all(mCommFolder, ec);
        CO_SIM_IO_INFO_IF("CoSimIO", ec) << "Warning, communication directory (" << mCommFolder << ") could not be deleted!\nError code: " << ec.message() << std::endl;
        if (!fs::exists(mCommFolder)) {
            fs::create_directory(mCommFolder);
            AddFilePermissions(mCommFolder); // otherwise the process with lower rights cannot delete files in it
        }
    }

    SynchronizeAll("conn");

    CO_SIM_IO_CATCH
}

void Communication::BaseDisconnectDetail(const Info& I_Info)
{
    CO_SIM_IO_TRY

    SynchronizeAll("disconn");

    if (mCommInFolder && GetIsPrimaryConnection() && mpDataComm->Rank() == 0) {
        // delete directory to remove potential leftovers
        std::error_code ec;
        fs::remove_all(mCommFolder, ec);
        if (ec) {
            CO_SIM_IO_INFO("CoSimIO") << "Warning, communication directory (" << mCommFolder << ") could not be deleted!\nError code: " << ec.message() << std::endl;
        }
    }

    CO_SIM_IO_CATCH
}

Info Communication::ImportInfoImpl(const Info& I_Info)
{
    CO_SIM_IO_TRY

    Info imported_info;
    Info rec_info = ReceiveObjectWithStreamSerializer(I_Info, imported_info);
    imported_info.Set<double>("elapsed_time", rec_info.Get<double>("elapsed_time"));
    imported_info.Set<double>("elapsed_time_ipc", rec_info.Get<double>("elapsed_time_ipc"));
    imported_info.Set<double>("elapsed_time_serializer", rec_info.Get<double>("elapsed_time_serializer"));
    imported_info.Set<std::size_t>("memory_usage_ipc", rec_info.Get<std::size_t>("memory_usage_ipc"));
    return imported_info;

    CO_SIM_IO_CATCH
}

Info Communication::ExportInfoImpl(const Info& I_Info)
{
    CO_SIM_IO_TRY

    return SendObjectWithStreamSerializer(I_Info, I_Info);

    CO_SIM_IO_CATCH
}

Info Communication::ImportDataImpl(
    const Info& I_Info,
    Internals::DataContainer<double>& rData)
{
    CO_SIM_IO_TRY

    if (mAlwaysUseSerializer) {
        return ReceiveObjectWithStreamSerializer(I_Info, rData);
    } else {
        Info info;
        const double elapsed_time = ReceiveDataContainer(I_Info, rData);
        info.Set<double>("elapsed_time", elapsed_time);
        info.Set<std::size_t>("memory_usage_ipc", rData.size()*sizeof(double));
        return info;
    }

    CO_SIM_IO_CATCH
}

Info Communication::ExportDataImpl(
    const Info& I_Info,
    const Internals::DataContainer<double>& rData)
{
    CO_SIM_IO_TRY

    if (mAlwaysUseSerializer) {
        return SendObjectWithStreamSerializer(I_Info, rData);
    } else {
        Info info;
        const double elapsed_time = SendDataContainer(I_Info, rData);
        info.Set<double>("elapsed_time", elapsed_time);
        info.Set<std::size_t>("memory_usage_ipc", rData.size()*sizeof(double));
        return info;
    }

    CO_SIM_IO_CATCH
}

Info Communication::ImportMeshImpl(
    const Info& I_Info,
    ModelPart& O_ModelPart)
{
    CO_SIM_IO_TRY

    return ReceiveObjectWithStreamSerializer(I_Info, O_ModelPart);

    CO_SIM_IO_CATCH
}

Info Communication::ExportMeshImpl(
    const Info& I_Info,
    const ModelPart& I_ModelPart)
{
    CO_SIM_IO_TRY

    return SendObjectWithStreamSerializer(I_Info, I_ModelPart);

    CO_SIM_IO_CATCH
}

void Communication::CheckConnection(const Info& I_Info)
{
    CO_SIM_IO_ERROR_IF_NOT(mIsConnected) << "No active connection exists!" << std::endl;
    CO_SIM_IO_ERROR_IF_NOT(I_Info.Has("identifier")) << "\"identifier\" must be specified!" << std::endl;
    Utilities::CheckEntry(I_Info.Get<std::string>("identifier"), "identifier");
}

void Communication::PostChecks(const Info& I_Info)
{
    CO_SIM_IO_ERROR_IF_NOT(I_Info.Has("elapsed_time")) << "\"elapsed_time\" must be specified!" << std::endl;
    CO_SIM_IO_ERROR_IF_NOT(I_Info.Has("memory_usage_ipc")) << "\"memory_usage_ipc\" must be specified!" << std::endl;
}

fs::path Communication::GetTmpFileName(
    const fs::path& rPath,
    const bool UseAuxFileForFileAvailability) const
{
    CO_SIM_IO_TRY

    if (!UseAuxFileForFileAvailability) {
        if (mCommInFolder) {
            return rPath.string().insert(mCommFolder.string().length()+1, ".");
        } else {
            return "." + rPath.string();
        }
    } else {
        return rPath;
    }

    CO_SIM_IO_CATCH
}

fs::path Communication::GetFileName(const fs::path& rPath, const std::string& rExtension) const
{
    CO_SIM_IO_TRY

    fs::path local_copy(rPath);
    local_copy += "." + rExtension;

    if (mCommInFolder) {
        return mCommFolder / local_copy;
    } else {
        return local_copy;
    }

    CO_SIM_IO_CATCH
}

fs::path Communication::GetFileName(const fs::path& rPath, const int Rank, const std::string& rExtension) const
{
    CO_SIM_IO_TRY

    fs::path local_copy(rPath);
    local_copy += "_s" + std::to_string(mpDataComm->Rank()) + "_d" + std::to_string(Rank);

    return GetFileName(local_copy, rExtension);

    CO_SIM_IO_CATCH
}

void Communication::WaitForPath(
    const fs::path& rPath,
    const bool UseAuxFileForFileAvailability,
    const int PrintEchoLevel) const
{
    CO_SIM_IO_TRY

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>=PrintEchoLevel) << "Waiting for: " << rPath << std::endl;
    if (!UseAuxFileForFileAvailability) {
        Utilities::WaitUntilPathExists(rPath);
    } else {
        fs::path avail_file = fs::path(rPath.string()+".avail");
        Utilities::WaitUntilPathExists(avail_file);

        // once the file exists it means that the real file was written, hence it can be removed
        RemovePath(avail_file);
    }
    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>=PrintEchoLevel) << "Found: " << rPath << std::endl;

    CO_SIM_IO_CATCH
}

void Communication::WaitUntilFileIsRemoved(
    const fs::path& rPath,
    const int PrintEchoLevel) const
{
    CO_SIM_IO_TRY

    std::error_code ec;
    if (fs::exists(rPath, ec)) { // only issue the wating message if the file exists initially
        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>=PrintEchoLevel) << "Waiting for: " << rPath << " to be removed" << std::endl;
        while(fs::exists(rPath, ec)) {
            std::this_thread::sleep_for(std::chrono::microseconds(10)); // wait 0.001s before next check
        }
        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>=PrintEchoLevel) << rPath << " was removed" << std::endl;
    }

    CO_SIM_IO_CATCH
}

void Communication::MakeFileVisible(
    const fs::path& rPath,
    const bool UseAuxFileForFileAvailability) const
{
    CO_SIM_IO_TRY

    if (!UseAuxFileForFileAvailability) {
        const fs::path tmp_file_name = GetTmpFileName(rPath, UseAuxFileForFileAvailability);
        AddFilePermissions(tmp_file_name);
        std::error_code ec;
        fs::rename(tmp_file_name, rPath, ec);
        CO_SIM_IO_ERROR_IF(ec) << rPath << " could not be made visible!\nError code: " << ec.message() << std::endl;
    } else {
        AddFilePermissions(rPath);
        std::ofstream avail_file;
        avail_file.open(rPath.string() + ".avail");
        avail_file.close();
    }

    CO_SIM_IO_CATCH
}

void Communication::RemovePath(const fs::path& rPath) const
{
    CO_SIM_IO_TRY

    // In windows the file cannot be removed if another file handle is using it
    // this can be the case here if the partner checks if the file (still) exists
    // hence we try multiple times to delete it
    std::error_code ec;
    for (std::size_t i=0; i<5; ++i) {
        if (fs::remove(rPath, ec)) {
            return; // if file could be removed succesfully then return
        }
    }
    CO_SIM_IO_ERROR << rPath << " could not be deleted!\nError code: " << ec.message() << std::endl;

    CO_SIM_IO_CATCH
}

void Communication::SynchronizeAll(const std::string& rTag) const
{
    CO_SIM_IO_TRY

    // first synchronize among the partitions
    mpDataComm->Barrier();

    // then synchronize among the partners
    if (mpDataComm->Rank() == 0) {
        const fs::path file_name_primary(GetFileName("CoSimIO_primary_" + GetConnectionName() + "_" + rTag, "sync"));
        const fs::path file_name_secondary(GetFileName("CoSimIO_secondary_" + GetConnectionName() + "_" + rTag, "sync"));

        if (GetIsPrimaryConnection()) {
            std::ofstream sync_file;
            sync_file.open(GetTmpFileName(file_name_primary));
            sync_file.close();
            CO_SIM_IO_ERROR_IF_NOT(fs::exists(GetTmpFileName(file_name_primary))) << "Primary sync file " << file_name_primary << " could not be created!" << std::endl;
            MakeFileVisible(file_name_primary);

            WaitForPath(file_name_secondary, true, 2);
            RemovePath(file_name_secondary);

            WaitUntilFileIsRemoved(file_name_primary, 2);
        } else {
            WaitForPath(file_name_primary, true, 2);
            RemovePath(file_name_primary);

            std::ofstream sync_file;
            sync_file.open(GetTmpFileName(file_name_secondary));
            sync_file.close();
            CO_SIM_IO_ERROR_IF_NOT(fs::exists(GetTmpFileName(file_name_secondary))) << "Secondary sync file " << file_name_secondary << " could not be created!" << std::endl;
            MakeFileVisible(file_name_secondary);

            WaitUntilFileIsRemoved(file_name_secondary, 2);
        }
    }

    // and finally synchronize again among the partitions
    mpDataComm->Barrier();

    CO_SIM_IO_CATCH
}

Info Communication::GetMyInfo() const
{
    CoSimIO::Info my_info;

    my_info.Set<int>("version_major", GetMajorVersion());
    my_info.Set<int>("version_minor", GetMinorVersion());
    my_info.Set<std::string>("version_patch", GetPatchVersion());

    my_info.Set<bool>("primary_was_explicitly_specified", mPrimaryWasExplicitlySpecified);
    my_info.Set<std::string>("communication_format", GetCommunicationName());
    my_info.Set<std::string>("operating_system", GetOsName());

    my_info.Set<int>("num_processes",   GetDataCommunicator().Size());

    my_info.Set<bool>("is_big_endian", Utilities::IsBigEndian());

    my_info.Set<bool>("always_use_serializer", mAlwaysUseSerializer);
    my_info.Set<std::string>("serializer_trace_type", Serializer::TraceTypeToString(mSerializerTraceType));

    my_info.Set<Info>("communication_settings", GetCommunicationSettings());

    return my_info;
}

void Communication::HandShake(const Info& I_Info)
{
    CO_SIM_IO_TRY

    if (mpDataComm->Rank() == 0) {
        const fs::path file_name_p2s(GetFileName("CoSimIO_" + GetConnectionName() + "_compatibility_check_primary_to_secondary", "dat"));
        const fs::path file_name_s2p(GetFileName("CoSimIO_" + GetConnectionName() + "_compatibility_check_secondary_to_primary", "dat"));

        auto exchange_data_for_handshake = [this](
            const fs::path& rMyFileName, const fs::path& rOtherFileName){

            // first export my info
            WaitUntilFileIsRemoved(rMyFileName,1); // in case of leftovers

            { // necessary as FileSerializer releases resources on destruction!
                FileSerializer serializer_save(GetTmpFileName(rMyFileName).string(), mSerializerTraceType);
                serializer_save.save("info", GetMyInfo());
            }

            MakeFileVisible(rMyFileName);

            // now get the info from the other
            WaitForPath(rOtherFileName, true, 1);

            { // necessary as FileSerializer releases resources on destruction!
                FileSerializer serializer_load(rOtherFileName.string(), mSerializerTraceType);
                serializer_load.load("info", mPartnerInfo);
            }

            RemovePath(rOtherFileName);
        };

        if (GetIsPrimaryConnection()) {
            exchange_data_for_handshake(file_name_p2s, file_name_s2p);
        } else {
            exchange_data_for_handshake(file_name_s2p, file_name_p2s);
        }

        // perform checks for compatibility
        CO_SIM_IO_ERROR_IF(GetMajorVersion() != mPartnerInfo.Get<int>("version_major")) << "Major version mismatch! My version: " << GetMajorVersion() << "; partner version: " << mPartnerInfo.Get<int>("version_major") << std::endl;
        CO_SIM_IO_ERROR_IF(GetMinorVersion() != mPartnerInfo.Get<int>("version_minor")) << "Minor version mismatch! My version: " << GetMinorVersion() << "; partner version: " << mPartnerInfo.Get<int>("version_minor") << std::endl;

        CO_SIM_IO_ERROR_IF(mPrimaryWasExplicitlySpecified != mPartnerInfo.Get<bool>("primary_was_explicitly_specified")) << std::boolalpha << "Mismatch in how the primary connection was specified!\nPrimary connection was explicitly specified for me: " << mPrimaryWasExplicitlySpecified << "\nPrimary connection was explicitly specified for partner: " << mPartnerInfo.Get<bool>("primary_was_explicitly_specified") << std::noboolalpha << std::endl;

        CO_SIM_IO_ERROR_IF(GetCommunicationName() != mPartnerInfo.Get<std::string>("communication_format")) << "Mismatch in communication_format!\nMy communication_format: " << GetCommunicationName() << "\nPartner communication_format: " << mPartnerInfo.Get<std::string>("communication_format") << std::endl;

        CO_SIM_IO_ERROR_IF(GetDataCommunicator().Size() != mPartnerInfo.Get<int>("num_processes")) << "Mismatch in num_processes!\nMy num_processes: " << GetDataCommunicator().Size() << "\nPartner num_processes: " << mPartnerInfo.Get<int>("num_processes") << std::endl;

        CO_SIM_IO_ERROR_IF(mAlwaysUseSerializer != mPartnerInfo.Get<bool>("always_use_serializer")) << std::boolalpha << "Mismatch in always_use_serializer!\nMy always_use_serializer: " << mAlwaysUseSerializer << "\nPartner always_use_serializer: " << mPartnerInfo.Get<bool>("always_use_serializer") << std::noboolalpha << std::endl;

        CO_SIM_IO_ERROR_IF(Serializer::TraceTypeToString(mSerializerTraceType) != mPartnerInfo.Get<std::string>("serializer_trace_type")) << "Mismatch in serializer_trace_type!\nMy serializer_trace_type: " << Serializer::TraceTypeToString(mSerializerTraceType) << "\nPartner serializer_trace_type: " << mPartnerInfo.Get<std::string>("serializer_trace_type") << std::endl;

        auto print_endianness = [](const bool IsBigEndian){return IsBigEndian ? "big endian" : "small endian";};

        CO_SIM_IO_INFO_IF("CoSimIO", Utilities::IsBigEndian() != mPartnerInfo.Get<bool>("is_big_endian")) << "WARNING: Parnters have different endianness, check results carefully! It is recommended to use serialized ascii communication.\n    My endianness:      " << print_endianness(Utilities::IsBigEndian()) << "\n    Partner endianness: " << print_endianness(mPartnerInfo.Get<bool>("is_big_endian")) << std::endl;

        // more things can be done in derived class if necessary
        DerivedHandShake();
    }

    // sync the partner info among the partitions
    mpDataComm->Broadcast(mPartnerInfo, 0);

    CO_SIM_IO_CATCH
}

void Communication::PrintElapsedTime(
    const Info& I_Info,
    const Info& O_Info,
    const std::string& rLabel)
{
    const std::string identifier =I_Info.Get<std::string>("identifier");
    const double dur = O_Info.Get<double>("elapsed_time");
    CO_SIM_IO_INFO_IF("CoSimIO-Timing", GetPrintTiming() && GetDataCommunicator().Rank()==0) << rLabel << " \"" << identifier << "\" took " << dur << " [s]" << std::endl;
}

} // namespace Internals
} // namespace CoSimIO
