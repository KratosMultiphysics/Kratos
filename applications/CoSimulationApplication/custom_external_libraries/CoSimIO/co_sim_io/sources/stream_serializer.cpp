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
#include "includes/stream_serializer.hpp"

namespace CoSimIO {
namespace Internals {

StreamSerializer::StreamSerializer(TraceType const& rTrace)
    : Serializer(new std::stringstream(std::ios::binary|std::ios::in|std::ios::out),rTrace)
{
    // increase precision when using ascii
    if (rTrace != SERIALIZER_NO_TRACE) {
        *pGetBuffer() << std::setprecision(14);
    }
}

StreamSerializer::StreamSerializer(const std::string& rData,TraceType const& rTrace)
    : StreamSerializer(rTrace)
{
    *(this->pGetBuffer()) << rData << std::endl;
}

} // namespace Internals
} // namespace CoSimIO
