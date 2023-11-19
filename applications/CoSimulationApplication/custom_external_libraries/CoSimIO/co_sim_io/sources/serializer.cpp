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

// Project includes
#include "includes/serializer.hpp"

namespace CoSimIO {
namespace Internals {

std::string Serializer::TraceTypeToString(Serializer::TraceType Trace)
{
    if (Trace == Serializer::TraceType::SERIALIZER_NO_TRACE) {return "no_trace";}
    else if (Trace == Serializer::TraceType::SERIALIZER_TRACE_ERROR) {return "trace_error";}
    else if (Trace == Serializer::TraceType::SERIALIZER_TRACE_ALL) {return "trace_all";}
    else if (Trace == Serializer::TraceType::SERIALIZER_ASCII) {return "ascii";}
    else {CO_SIM_IO_ERROR << "Invalid trace_type! Valid options are: no_trace, trace_error, trace_all, ascii" << std::endl;}
}

Serializer::TraceType Serializer::StringToTraceType(const std::string& Trace)
{
    if (Trace == "no_trace") {return Serializer::TraceType::SERIALIZER_NO_TRACE;}
    else if (Trace == "trace_error") {return Serializer::TraceType::SERIALIZER_TRACE_ERROR;}
    else if (Trace == "trace_all") {return Serializer::TraceType::SERIALIZER_TRACE_ALL;}
    else if (Trace == "ascii") {return Serializer::TraceType::SERIALIZER_ASCII;}
    else {CO_SIM_IO_ERROR << "Invalid trace_type! Valid options are: no_trace, trace_error, trace_all, ascii" << std::endl;}
}

Serializer::RegisteredObjectsContainerType Serializer::msRegisteredObjects;
Serializer::RegisteredObjectsNameContainerType Serializer::msRegisteredObjectsName;

} // namespace Internals
}
