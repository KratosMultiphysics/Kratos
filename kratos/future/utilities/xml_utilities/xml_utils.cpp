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

// System includes
#include <numeric>

// Project includes

// Include base h
#include "xml_utils.h"

namespace Kratos::Future {

template<class TDataType>
void XmlUtilities::AddDataArrayAttributes(
    XmlElement& rXmlElement,
    const NDData<TDataType>& rNDData)
{
    KRATOS_TRY

    // add the data type information
    if constexpr(std::is_same_v<TDataType, unsigned char>) {
        rXmlElement.AddAttribute("type", "UInt" + std::to_string(sizeof(TDataType) * 8));
    } else if constexpr(std::is_same_v<TDataType, bool>) {
        rXmlElement.AddAttribute("type", "UInt" + std::to_string(sizeof(TDataType) * 8));
    } else if constexpr(std::is_same_v<TDataType, int>) {
        rXmlElement.AddAttribute("type", "Int" + std::to_string(sizeof(TDataType) * 8));
    } else if constexpr(std::is_same_v<TDataType, double>) {
        rXmlElement.AddAttribute("type", "Float" + std::to_string(sizeof(TDataType) * 8));
    } else {
        KRATOS_ERROR << "Unsupported data type.";
    }

    const auto& shape = rNDData.Shape();

    KRATOS_ERROR_IF(shape.size() == 0)
        << "NDData must have at least one dimension representing number of entities.\n";

    rXmlElement.AddAttribute("NumberOfComponents", std::to_string(std::accumulate(shape.begin() + 1, shape.end(), 1u, std::multiplies<unsigned int>{})));

    KRATOS_CATCH("");
}

// template instantiations
template void XmlUtilities::AddDataArrayAttributes(XmlElement&, const NDData<unsigned char>&);
template void XmlUtilities::AddDataArrayAttributes(XmlElement&, const NDData<bool>&);
template void XmlUtilities::AddDataArrayAttributes(XmlElement&, const NDData<int>&);
template void XmlUtilities::AddDataArrayAttributes(XmlElement&, const NDData<double>&);

} // namespace Kratos::Future