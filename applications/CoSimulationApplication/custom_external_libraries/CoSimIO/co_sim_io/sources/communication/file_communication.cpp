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
#include <chrono>
#include <iomanip>

// Project includes
#include "includes/communication/file_communication.hpp"
#include "includes/utilities.hpp"
#include "includes/file_serializer.hpp"

namespace CoSimIO {
namespace Internals {

namespace {

// Important: having this in a function also makes sure that the FileSerializer releases its resources (i.e. the file) at destruction
template<class TObject>
void SerializeToFile(
    const fs::path& rPath,
    const TObject& rObject,
    const Serializer::TraceType SerializerTrace)
{
    CO_SIM_IO_TRY

    FileSerializer serializer(rPath.string(), SerializerTrace);
    serializer.save("obj", rObject);

    CO_SIM_IO_CATCH
}

// important: having this in a function also makes sure that the FileSerializer releases its resources (i.e. the file) at destruction
template<class TObject>
void SerializeFromFile(
    const fs::path& rPath,
    TObject& rObject,
    const Serializer::TraceType SerializerTrace)
{
    CO_SIM_IO_TRY

    FileSerializer serializer(rPath.string(), SerializerTrace);
    serializer.load("obj", rObject);

    CO_SIM_IO_CATCH
}

}

FileCommunication::FileCommunication(
    const Info& I_Settings,
    std::shared_ptr<DataCommunicator> I_DataComm)
    : Communication(I_Settings, I_DataComm),
      mUseAuxFileForFileAvailability(I_Settings.Get<bool>("use_aux_file_for_file_availability", USE_AUX_FILE_FOR_FILE_AVAILABILITY)),
      mUseFileSerializer(I_Settings.Get<bool>("use_file_serializer", true))
{
#ifdef CO_SIM_IO_COMPILED_IN_WINDOWS
    CO_SIM_IO_INFO_IF("CoSimIO", !mUseAuxFileForFileAvailability) << "WARNING: Using rename for making files available can cause race conditions as it is not atomic in Windows! Use \"use_aux_file_for_file_availability\" = false to avoid this" << std::endl;
#endif
}

FileCommunication::~FileCommunication()
{
    CO_SIM_IO_TRY

    if (GetIsConnected()) {
        CO_SIM_IO_INFO("CoSimIO") << "Warning: Disconnect was not performed, attempting automatic disconnection!" << std::endl;
        Info tmp;
        Disconnect(tmp);
    }

    CO_SIM_IO_CATCH
}

void FileCommunication::DerivedHandShake() const
{
    CO_SIM_IO_TRY

    const bool my_use_aux_file_for_file_availability = GetMyInfo().Get<Info>("communication_settings").Get<bool>("use_aux_file_for_file_availability");
    const bool partner_use_aux_file_for_file_availability = GetPartnerInfo().Get<Info>("communication_settings").Get<bool>("use_aux_file_for_file_availability");
    CO_SIM_IO_ERROR_IF(my_use_aux_file_for_file_availability != partner_use_aux_file_for_file_availability) << std::boolalpha << "Mismatch in use_aux_file_for_file_availability!\nMy use_aux_file_for_file_availability: " << my_use_aux_file_for_file_availability << "\nPartner use_aux_file_for_file_availability: " << partner_use_aux_file_for_file_availability << std::noboolalpha << "\nNote that the default on unix is false, while the default on windows is true!"<< std::endl;

    const bool my_use_file_serializer = GetMyInfo().Get<Info>("communication_settings").Get<bool>("use_file_serializer");
    const bool partner_use_file_serializer = GetPartnerInfo().Get<Info>("communication_settings").Get<bool>("use_file_serializer");
    CO_SIM_IO_ERROR_IF(my_use_file_serializer != partner_use_file_serializer) << std::boolalpha << "Mismatch in use_file_serializer!\nMy use_file_serializer: " << my_use_file_serializer << "\nPartner use_file_serializer: " << partner_use_file_serializer << std::noboolalpha << std::endl;

    CO_SIM_IO_CATCH
}

Info FileCommunication::GetCommunicationSettings() const
{
    CO_SIM_IO_TRY

    Info info;
    info.Set("use_aux_file_for_file_availability", mUseAuxFileForFileAvailability);
    info.Set("use_file_serializer", mUseFileSerializer);

    return info;

    CO_SIM_IO_CATCH
}

Info FileCommunication::ImportInfoImpl(const Info& I_Info)
{
    CO_SIM_IO_TRY

    if (mUseFileSerializer) {
        Info received_info;
        const Info rec_info = GenericReceiveWithFileSerializer(I_Info, received_info);
        received_info.Set<double>("elapsed_time", rec_info.Get<double>("elapsed_time"));
        received_info.Set<std::size_t>("memory_usage_ipc", rec_info.Get<std::size_t>("memory_usage_ipc"));
        return received_info;
    }

    return Communication::ImportInfoImpl(I_Info);

    CO_SIM_IO_CATCH
}

Info FileCommunication::ExportInfoImpl(const Info& I_Info)
{
    CO_SIM_IO_TRY

    if (mUseFileSerializer) {
        return GenericSendWithFileSerializer(I_Info, I_Info);
    }

    return Communication::ExportInfoImpl(I_Info);

    CO_SIM_IO_CATCH
}

Info FileCommunication::ImportDataImpl(
    const Info& I_Info,
    Internals::DataContainer<double>& rData)
{
    CO_SIM_IO_TRY

    if (mUseFileSerializer && GetAlwaysUseSerializer()) {
        return GenericReceiveWithFileSerializer(I_Info, rData);
    }

    return Communication::ImportDataImpl(I_Info, rData);

    CO_SIM_IO_CATCH
}

Info FileCommunication::ExportDataImpl(
    const Info& I_Info,
    const Internals::DataContainer<double>& rData)
{
    CO_SIM_IO_TRY

    if (mUseFileSerializer && GetAlwaysUseSerializer()) {
        return GenericSendWithFileSerializer(I_Info, rData);
    }

    return Communication::ExportDataImpl(I_Info, rData);

    CO_SIM_IO_CATCH
}

Info FileCommunication::ImportMeshImpl(
    const Info& I_Info,
    ModelPart& O_ModelPart)
{
    CO_SIM_IO_TRY

    if (mUseFileSerializer) {
        return GenericReceiveWithFileSerializer(I_Info, O_ModelPart);
    }

    return Communication::ImportMeshImpl(I_Info, O_ModelPart);

    CO_SIM_IO_CATCH
}

Info FileCommunication::ExportMeshImpl(
    const Info& I_Info,
    const ModelPart& I_ModelPart)
{
    CO_SIM_IO_TRY

    if (mUseFileSerializer) {
        return GenericSendWithFileSerializer(I_Info, I_ModelPart);
    }

    return Communication::ExportMeshImpl(I_Info, I_ModelPart);

    CO_SIM_IO_CATCH
}

template<class TObjectType>
Info FileCommunication::GenericSendWithFileSerializer(
    const Info& I_Info,
    const TObjectType& rObj)
{
    CO_SIM_IO_TRY

    Info info;

    const std::string identifier = I_Info.Get<std::string>("identifier");

    const fs::path file_name(GetFileName("CoSimIO_data_" + GetConnectionName() + "_" + identifier + "_" + std::to_string(GetDataCommunicator().Rank()), "dat"));

    WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

    const auto start_time(std::chrono::steady_clock::now());
    SerializeToFile(GetTmpFileName(file_name, mUseAuxFileForFileAvailability), rObj, GetSerializerTraceType());

    info.Set<std::size_t>("memory_usage_ipc", fs::file_size(GetTmpFileName(file_name, mUseAuxFileForFileAvailability)));

    MakeFileVisible(file_name, mUseAuxFileForFileAvailability);

    const double elapsed_time = Utilities::ElapsedSeconds(start_time);
    info.Set<double>("elapsed_time", elapsed_time);
    return info;

    CO_SIM_IO_CATCH
}

template<class TObjectType>
Info FileCommunication::GenericReceiveWithFileSerializer(
    const Info& I_Info,
    TObjectType& rObj)
{
    CO_SIM_IO_TRY

    Info info;

    const std::string identifier = I_Info.Get<std::string>("identifier");

    const fs::path file_name(GetFileName("CoSimIO_data_" + GetConnectionName() + "_" + identifier + "_" + std::to_string(GetDataCommunicator().Rank()), "dat"));

    WaitForPath(file_name, mUseAuxFileForFileAvailability);

    info.Set<std::size_t>("memory_usage_ipc", fs::file_size(file_name));

    const auto start_time(std::chrono::steady_clock::now());
    SerializeFromFile(file_name, rObj, GetSerializerTraceType());

    RemovePath(file_name);

    const double elapsed_time = Utilities::ElapsedSeconds(start_time);
    info.Set<double>("elapsed_time", elapsed_time);
    return info;

    CO_SIM_IO_CATCH
}

template<typename T>
double FileCommunication::GenericSend(
    const Info& I_Info,
    const T& rData,
    const int SizeOfData)
{
    CO_SIM_IO_TRY

    const std::string identifier = I_Info.Get<std::string>("identifier");

    const fs::path file_name(GetFileName("CoSimIO_data_" + GetConnectionName() + "_" + identifier + "_" + std::to_string(GetDataCommunicator().Rank()), "dat"));

    WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

    const std::size_t size = rData.size();

    const auto start_time(std::chrono::steady_clock::now());

    std::ofstream output_file(GetTmpFileName(file_name, mUseAuxFileForFileAvailability), std::ios::out|std::ios::binary);
    Utilities::CheckStream(output_file, file_name);

    output_file.write(reinterpret_cast<const char *>(&size), sizeof(std::size_t));

    output_file.write(reinterpret_cast<const char *>(&rData[0]), rData.size()*SizeOfData);

    output_file.close();
    MakeFileVisible(file_name, mUseAuxFileForFileAvailability);

    return Utilities::ElapsedSeconds(start_time);

    CO_SIM_IO_CATCH
}

template<typename T>
double FileCommunication::GenericReceive(
    const Info& I_Info,
    T& rData,
    const int SizeOfData)
{
    CO_SIM_IO_TRY

    const std::string identifier = I_Info.Get<std::string>("identifier");

    const fs::path file_name(GetFileName("CoSimIO_data_" + GetConnectionName() + "_" + identifier + "_" + std::to_string(GetDataCommunicator().Rank()), "dat"));

    WaitForPath(file_name, mUseAuxFileForFileAvailability);

    const auto start_time(std::chrono::steady_clock::now());

    std::ifstream input_file(file_name, std::ios::binary|std::ios::in);
    Utilities::CheckStream(input_file, file_name);

    std::size_t size_read;
    input_file.read((char*)&size_read, sizeof(std::size_t));

    rData.resize(size_read);
    input_file.read((char*)&rData[0], size_read*SizeOfData);

    input_file.close();
    RemovePath(file_name);

    return Utilities::ElapsedSeconds(start_time);

    CO_SIM_IO_CATCH
}

double FileCommunication::SendString(
    const Info& I_Info,
    const std::string& rData)
{
    return GenericSend(I_Info, rData, 1);
}

double FileCommunication::ReceiveString(
    const Info& I_Info,
    std::string& rData)
{
    return GenericReceive(I_Info, rData, 1);
}

double FileCommunication::SendDataContainer(
    const Info& I_Info,
    const Internals::DataContainer<double>& rData)
{
    return GenericSend(I_Info, rData, sizeof(double));
}

double FileCommunication::ReceiveDataContainer(
    const Info& I_Info,
    Internals::DataContainer<double>& rData)
{
    return GenericReceive(I_Info, rData, sizeof(double));
}

} // namespace Internals
} // namespace CoSimIO
