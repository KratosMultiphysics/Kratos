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
#include <mutex>

// Project includes
#include "includes/info.hpp"

namespace CoSimIO {

namespace {
    std::mutex register_mutex;
}

void Info::Print(
    std::ostream& rOStream,
    const std::string& rPrefixString) const
{
    rOStream << "CoSimIO-Info; containing " << Size() << " entries\n";

    for (const auto& r_pair: mOptions) {
        rOStream << rPrefixString << "  name: " << r_pair.first << " | ";
        r_pair.second->Print(rOStream, rPrefixString + "  ");
    }
}

void Info::save(CoSimIO::Internals::Serializer& rSerializer) const
{
    CO_SIM_IO_TRY

    RegisterTypesInSerializer();

    rSerializer.save("mOptions", mOptions);

    CO_SIM_IO_CATCH
}
void Info::load(CoSimIO::Internals::Serializer& rSerializer)
{
    CO_SIM_IO_TRY

    RegisterTypesInSerializer();
    rSerializer.load("mOptions", mOptions);

    CO_SIM_IO_CATCH
}

void Info::RegisterTypesInSerializer()
{
    CO_SIM_IO_TRY

    if (!mpSerializerTypesRegistered) {
        const std::lock_guard<std::mutex> scope_lock(register_mutex);
        if (!mpSerializerTypesRegistered) {
            static bool types_are_registered; // no need to initialize, all that matters is that mpSerializerTypesRegistered is no longer a null_ptr
            mpSerializerTypesRegistered = &types_are_registered;

            static Internals::InfoData<int> info_data_int(1);
            static Internals::InfoData<double> info_data_double(1);
            static Internals::InfoData<bool> info_data_bool(1);
            static Internals::InfoData<std::string> info_data_string("");
            static Internals::InfoData<Info> info_data_info(Info{});

            CoSimIO::Internals::Serializer::Register("info_data_int",    info_data_int);
            CoSimIO::Internals::Serializer::Register("info_data_double", info_data_double);
            CoSimIO::Internals::Serializer::Register("info_data_bool",   info_data_bool);
            CoSimIO::Internals::Serializer::Register("info_data_string", info_data_string);
            CoSimIO::Internals::Serializer::Register("info_data_info",   info_data_info);
        }
    }

    CO_SIM_IO_CATCH
}

bool* Info::mpSerializerTypesRegistered = nullptr;

} // namespace CoSimIO
