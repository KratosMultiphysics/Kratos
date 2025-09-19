//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <cstdint>

// Project includes
#include "containers/nd_data.h"
#include "includes/define.h"
#include "xml_element.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) XmlUtilities
{
public:
    ///@name Public static operations
    ///@{

    template<class TDataType>
    void static AddDataArrayAttributes(
        XmlElement& rXmlElement,
        const NDData<TDataType>& rNDData);


    template<class TDataType>
    void static AppendData(
        std::vector<char>& rOutput,
        TDataType const * pData,
        const std::uint64_t DataSize,
        std::uint64_t const * pHeader,
        const std::uint64_t HeaderSize)
    {
        auto data_type_size = sizeof(TDataType);

        // calculate number of bytes in the data
        const std::uint64_t total_data_length = data_type_size * DataSize;
        const std::uint64_t total_header_length = sizeof(std::uint64_t) * HeaderSize;

        auto p_casted_header = reinterpret_cast<char const *>(pHeader);
        auto p_casted_data = reinterpret_cast<char const *>(pData);

        // now we reserve space for the header to write the size of pData, and
        // for the pData
        const auto current_size = rOutput.size();
        rOutput.resize(rOutput.size() + total_header_length + total_data_length);
        std::copy(p_casted_header, p_casted_header + total_header_length, rOutput.begin() + current_size);
        std::copy(p_casted_data, p_casted_data + total_data_length, rOutput.begin() + current_size + total_header_length);
    }

    ///@}
};

///@}

} // namespace Kratos