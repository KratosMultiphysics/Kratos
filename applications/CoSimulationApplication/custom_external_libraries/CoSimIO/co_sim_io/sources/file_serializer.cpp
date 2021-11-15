//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <iomanip>

// Project includes
#include "includes/file_serializer.hpp"

namespace CoSimIO {
namespace Internals {

FileSerializer::FileSerializer(std::string const& rFilename, Serializer::TraceType const& rTrace)
    : Serializer(nullptr, rTrace)
{
    std::fstream* p_file = new std::fstream(rFilename, std::ios::binary|std::ios::in|std::ios::out);

    if (!(*p_file)) {
        delete p_file;
        p_file = new std::fstream(rFilename, std::ios::binary|std::ios::out);
    }

    SetBuffer( p_file );

    CO_SIM_IO_ERROR_IF_NOT(*pGetBuffer()) << "Error opening input file: " << rFilename << std::endl;

    // increase precision when using ascii
    if (rTrace != SERIALIZER_NO_TRACE) {
        *pGetBuffer() << std::setprecision(14);
    }
}

} // namespace Internals
} // namespace CoSimIO
