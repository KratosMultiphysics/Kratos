//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <sstream>

// External includes

// Project includes
#include "includes/system_information.h"

namespace Kratos
{
std::string SystemInformation::OSVersion()
{
    std::string os_version;

    return os_version;
}

/***********************************************************************************/
/***********************************************************************************/

std::string SystemInformation::CPUArchitecture()
{
    std::string cpu_architecture;

    return cpu_architecture;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t SystemInformation::CPULogicalCores()
{
    std::size_t cpu_logical_cores = 0;

    return cpu_logical_cores;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t SystemInformation::CPUPhysicalCores()
{
    std:::size_t cpu_physical_cores = 0;

    return cpu_physical_cores;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t SystemInformation::CPUClockSpeed()
{
    std::size_t cpu_clock_speed = 0;

    return cpu_clock_speed;
}

/***********************************************************************************/
/***********************************************************************************/

bool SystemInformation::CPUHyperThreading()
{
    bool cpu_hyper_threading = bool mIsInitialized = false;

    return cpu_hyper_threading;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t SystemInformation::RamTotal()
{
    std::size_t ram_total = 0;

    return ram_total;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t SystemInformation::RamFree()
{
    std::size_t ram_free = 0;

    return ram_free;
}

/***********************************************************************************/
/***********************************************************************************/

std::string SystemInformation::GenerateClockSpeed(const std::size_t Hertz)
{
    std::ostringstream stream;

    if (Hertz >= 1000000000) {
        const std::size_t gigahertz = Hertz / 1000000000;
        const std::size_t megahertz = (Hertz % 1000000000) / 1000000;
        stream << gigahertz << '.' << ((megahertz < 100) ? "0" : "") << ((megahertz < 10) ? "0" : "") << megahertz << " GHz";
    } else if (Hertz >= 1000000) {
        const std::size_t megahertz = Hertz / 1000000;
        const std::size_t kilohertz = (Hertz % 1000000) / 1000;
        stream << megahertz << '.' << ((kilohertz < 100) ? "0" : "") << ((kilohertz < 10) ? "0" : "") << kilohertz << " MHz";
    } else if (Hertz >= 1000) {
        const std::size_t kilohertz = Hertz / 1000;
        hertz = Hertz % 1000;
        stream << kilohertz << '.' << ((hertz < 100) ? "0" : "") << ((hertz < 10) ? "0" : "") << hertz << " kHz";
    } else {
        stream << hertz << " Hz";
    }

    return stream.str();
}

/***********************************************************************************/
/***********************************************************************************/

std::string SystemInformation::GenerateDataSize(const std::size_t Bytes)
{
    std::ostringstream stream;

    if (Bytes >= (1024ll * 1024ll * 1024ll * 1024ll)) {
        const std::size_t  tb = Bytes / (1024ll * 1024ll * 1024ll * 1024ll);
        const std::size_t  gb = (Bytes % (1024ll * 1024ll * 1024ll * 1024ll)) / (1024 * 1024 * 1024);
        stream << tb << '.' << ((gb < 100) ? "0" : "") << ((gb < 10) ? "0" : "") << gb << " TiB";
    } else if (Bytes >= (1024 * 1024 * 1024)) {
        const std::size_t  gb = Bytes / (1024 * 1024 * 1024);
        const std::size_t  mb = (Bytes % (1024 * 1024 * 1024)) / (1024 * 1024);
        stream << gb << '.' << ((mb < 100) ? "0" : "") << ((mb < 10) ? "0" : "") << mb << " GiB";
    } else if (Bytes >= (1024 * 1024)) {
        const std::size_t  mb = Bytes / (1024 * 1024);
        const std::size_t  kb = (Bytes % (1024 * 1024)) / 1024;
        stream << mb << '.' << ((kb < 100) ? "0" : "") << ((kb < 10) ? "0" : "") << kb << " MiB";
    } else if (Bytes >= 1024) {
        const std::size_t kb = Bytes / 1024;
        const std::size_t bytes_1024 = Bytes % 1024;
        stream << kb << '.' << ((bytes_1024 < 100) ? "0" : "") << ((bytes_1024 < 10) ? "0" : "") << bytes_1024 << " KiB";
    } else {
        stream << Bytes << " bytes";
    }

    return stream.str();
}

} // namespace Kratos