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
#include <string>
#include <map>
#include <chrono>
#include <thread>
#include <cmath>
#include <system_error>

// Project includes
#include "includes/utilities.hpp"

namespace CoSimIO {
namespace Utilities {

// Create the name for the connection
// In a function bcs maybe in the future this will
// need to be more elaborate
std::string CreateConnectionName(
    const std::string& rName1,
    const std::string& rName2)
{
    if (rName1 < rName2) {
        return rName1 + "_" + rName2;
    } else {
        return rName2 + "_" + rName1;
    }
}

void CheckEntry(const std::string& rEntry, const std::string& rKey)
{
    // the entries are used e.g. for the folder creation, file names and other things
    // hence they are a little bit restricted to avoid unfortunate failures in rare cases

    const std::size_t max_allowed_size = 1000;
    CO_SIM_IO_ERROR_IF(rEntry.empty()) << "Using an empty entry for \"" << rKey << "\" is not allowed!" << std::endl;
    CO_SIM_IO_ERROR_IF(rEntry.length() > max_allowed_size) << "Entry for \"" << rKey << "\" is too long! Maximum allowed length: " << max_allowed_size << " characters!" << std::endl;

    const char disallowed_chars[] = {'.', ',', ':', ';', '>', '<', '/', '\'', '|', '*', '!', '"', ' '};
    for (const auto ch : disallowed_chars) {
        CO_SIM_IO_ERROR_IF_NOT(rEntry.find(ch) == std::string::npos) << "Entry for \"" << rKey << "\" contains a character that is not allowed: \"" << std::string(1,ch) << "\"!" << std::endl;
    }
}

double ElapsedSeconds(const std::chrono::steady_clock::time_point& rStartTime)
{
    using namespace std::chrono;
    return duration_cast<duration<double>>(steady_clock::now() - rStartTime).count();
}

bool IsBigEndian()
{
    // from: https://stackoverflow.com/a/1001373
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}

int GetNumberOfNodesForElementType(const ElementType I_ElementType)
{
    // using switch over map as the compiler warns if some enum values are not handled in the switch
    switch(I_ElementType) {
        case ElementType::Hexahedra3D20:
            return 20;
        case ElementType::Hexahedra3D27:
            return 27;
        case ElementType::Hexahedra3D8:
            return 8;
        case ElementType::Prism3D15:
            return 15;
        case ElementType::Prism3D6:
            return 6;
        case ElementType::Pyramid3D13:
            return 13;
        case ElementType::Pyramid3D5:
            return 5;
        case ElementType::Quadrilateral2D4:
            return 4;
        case ElementType::Quadrilateral2D8:
            return 8;
        case ElementType::Quadrilateral2D9:
            return 9;
        case ElementType::Quadrilateral3D4:
            return 4;
        case ElementType::Quadrilateral3D8:
            return 8;
        case ElementType::Quadrilateral3D9:
            return 9;
        case ElementType::Tetrahedra3D10:
            return 10;
        case ElementType::Tetrahedra3D4:
            return 4;
        case ElementType::Triangle2D3:
            return 3;
        case ElementType::Triangle2D6:
            return 6;
        case ElementType::Triangle3D3:
            return 3;
        case ElementType::Triangle3D6:
            return 6;
        case ElementType::Line2D2:
            return 2;
        case ElementType::Line2D3:
            return 3;
        case ElementType::Line3D2:
            return 2;
        case ElementType::Line3D3:
            return 3;
        case ElementType::Point2D:
            return 1;
        case ElementType::Point3D:
            return 1;
    };

    CO_SIM_IO_ERROR << "Unknown Element type!" << std::endl;
}

std::string GetElementName(const ElementType I_ElementType)
{
    // using switch over map as the compiler warns if some enum values are not handled in the switch
    switch(I_ElementType) {
        case ElementType::Hexahedra3D20:
            return "Hexahedra3D20";
        case ElementType::Hexahedra3D27:
            return "Hexahedra3D27";
        case ElementType::Hexahedra3D8:
            return "Hexahedra3D8";
        case ElementType::Prism3D15:
            return "Prism3D15";
        case ElementType::Prism3D6:
            return "Prism3D6";
        case ElementType::Pyramid3D13:
            return "Pyramid3D13";
        case ElementType::Pyramid3D5:
            return "Pyramid3D5";
        case ElementType::Quadrilateral2D4:
            return "Quadrilateral2D4";
        case ElementType::Quadrilateral2D8:
            return "Quadrilateral2D8";
        case ElementType::Quadrilateral2D9:
            return "Quadrilateral2D9";
        case ElementType::Quadrilateral3D4:
            return "Quadrilateral3D4";
        case ElementType::Quadrilateral3D8:
            return "Quadrilateral3D8";
        case ElementType::Quadrilateral3D9:
            return "Quadrilateral3D9";
        case ElementType::Tetrahedra3D10:
            return "Tetrahedra3D10";
        case ElementType::Tetrahedra3D4:
            return "Tetrahedra3D4";
        case ElementType::Triangle2D3:
            return "Triangle3D3";
        case ElementType::Triangle2D6:
            return "Triangle2D6";
        case ElementType::Triangle3D3:
            return "Triangle3D3";
        case ElementType::Triangle3D6:
            return "Triangle3D6";
        case ElementType::Line2D2:
            return "Line2D2";
        case ElementType::Line2D3:
            return "Line2D3";
        case ElementType::Line3D2:
            return "Line3D2";
        case ElementType::Line3D3:
            return "Line3D3";
        case ElementType::Point2D:
            return "Point2D";
        case ElementType::Point3D:
            return "Point3D";
    };

    CO_SIM_IO_ERROR << "Unknown Element type!" << std::endl;
}

void WaitUntilPathExists(const fs::path& rPath)
{
    std::error_code ec;
    while(!fs::exists(rPath, ec)) {std::this_thread::sleep_for(std::chrono::microseconds(10));} // wait 0.005s before next check
}

std::set<std::size_t> ComputePartnerRanksAsImporter(
    const std::size_t MyRank,
    const std::size_t MySize,
    const std::size_t PartnerSize)
{
    // validate input
    CO_SIM_IO_ERROR_IF(MySize==0) << "MySize cannot be zero!" << std::endl;
    CO_SIM_IO_ERROR_IF(PartnerSize==0) << "PartnerSize cannot be zero!" << std::endl;
    CO_SIM_IO_ERROR_IF_NOT(MyRank<MySize) << "MyRank must be smaller MySize!" << std::endl;

    if (MySize == 1) {
        // I am serial, communicate with all partner ranks (doesn't matter if partner is distributed or not)
        std::set<std::size_t> partner_ranks;
        for (std::size_t i=0; i<PartnerSize; ++i) {partner_ranks.insert(i);}
        return partner_ranks;
    } else if (MySize == PartnerSize) {
        // special case when both run with the same size
        return {MyRank};
    } else if (MySize > PartnerSize) {
        // partner is serial, only my rank 0 communicates with this rank
        if (MyRank < PartnerSize) {
            return {MyRank};
        } else {
            return {};
        }
    } else {
        // several of partner ranks communicate with one rank of me
        const std::size_t num_ranks_per_partner_rank = static_cast<std::size_t>(std::ceil(PartnerSize / static_cast<double>(MySize)));
        std::set<std::size_t> partner_ranks;
        const std::size_t lower_end = MyRank*num_ranks_per_partner_rank;
        const std::size_t upper_end = (MyRank+1)*num_ranks_per_partner_rank;

        for (std::size_t i=0; i<PartnerSize; ++i) {
            if (i >= lower_end && i < upper_end) {
                partner_ranks.insert(i);
            }
        }
        return partner_ranks;
    }
}

std::set<std::size_t> ComputePartnerRanksAsExporter(
    const std::size_t MyRank,
    const std::size_t MySize,
    const std::size_t PartnerSize)
{
    // validate input
    CO_SIM_IO_ERROR_IF(MySize==0) << "MySize cannot be zero!" << std::endl;
    CO_SIM_IO_ERROR_IF(PartnerSize==0) << "PartnerSize cannot be zero!" << std::endl;
    CO_SIM_IO_ERROR_IF_NOT(MyRank<MySize) << "MyRank must be smaller MySize!" << std::endl;

    if (MySize == 1) {
        return {0};
    } else if (PartnerSize == 1) {
        // partner is serial, all of my ranks communicate with one rank (rank 0)
        return {0};
    } else if (MySize == PartnerSize) {
        // special case when both run with the same size
        return {MyRank};
    } else {
        if (MySize > PartnerSize) {
            // several of my ranks communicate with one rank of partner
            const std::size_t num_ranks_per_partner_rank = static_cast<std::size_t>(std::ceil(MySize / static_cast<double>(PartnerSize)));
            return {MyRank/num_ranks_per_partner_rank};
        } else {
            // partner has more ranks, we only communicate with one
            return {MyRank};
        }
    }
}

} // namespace Utilities
} // namespace CoSimIO
