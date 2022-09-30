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

#ifndef CO_SIM_IO_UTILITIES_INCLUDED
#define CO_SIM_IO_UTILITIES_INCLUDED

// System includes
#include <string>
#include <chrono>
#include <set>

// Project includes
#include "define.hpp"
#include "filesystem_inc.hpp"

namespace CoSimIO {
namespace Utilities {

// Create the name for the connection
// In a function bcs maybe in the future this will
// need to be more elaborate
std::string CO_SIM_IO_API CreateConnectionName(
    const std::string& rName1,
    const std::string& rName2);

void CO_SIM_IO_API CheckEntry(const std::string& rEntry, const std::string& rKey);

template <typename TStream, typename TPath>
static void CheckStream(const TStream& rStream, const TPath& rPath)
{
    CO_SIM_IO_ERROR_IF_NOT(rStream.is_open()) << rPath << " could not be opened!" << std::endl;
}

double CO_SIM_IO_API ElapsedSeconds(const std::chrono::steady_clock::time_point& rStartTime);

// returns if the current system is big endian
bool CO_SIM_IO_API IsBigEndian();

int CO_SIM_IO_API GetNumberOfNodesForElementType(const ElementType I_ElementType);

std::string CO_SIM_IO_API GetElementName(const ElementType I_ElementType);

void CO_SIM_IO_API WaitUntilPathExists(const fs::path& rPath);

std::set<std::size_t> CO_SIM_IO_API ComputePartnerRanksAsImporter(
    const std::size_t MyRank,
    const std::size_t MySize,
    const std::size_t PartnerSize);

std::set<std::size_t> CO_SIM_IO_API ComputePartnerRanksAsExporter(
    const std::size_t MyRank,
    const std::size_t MySize,
    const std::size_t PartnerSize);

} // namespace Utilities
} // namespace CoSimIO

#endif // CO_SIM_IO_UTILITIES_INCLUDED
