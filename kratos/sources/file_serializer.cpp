//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "includes/file_serializer.h"

namespace Kratos
{
FileSerializer::FileSerializer(std::string const& Filename, Serializer::TraceType const& rTrace, const bool DataOnly)
    : Serializer(nullptr, rTrace, DataOnly)
{
    std::fstream* p_file = new std::fstream(std::string(Filename+".rest").c_str(), std::ios::binary|std::ios::in|std::ios::out);
    if(!(*p_file)) {
        delete p_file;
        p_file = new std::fstream(std::string(Filename+".rest").c_str(), std::ios::binary|std::ios::out);
    }
    SetBuffer( p_file );
    KRATOS_ERROR_IF_NOT(*pGetBuffer()) << "Error opening input file : "
                                << std::string(Filename+".rest") << std::endl;
}
}