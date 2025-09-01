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
#include <iomanip>
#include <type_traits>

// Project includes

// Include base h
#include "xml_nd_data_element.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TDataType>
XmlNDDataElement<TDataType>::XmlNDDataElement(
    const std::string& rDataArrayName,
    typename NDData<TDataType>::Pointer pNDData)
    : BaseType("DataArray"),
      mpNDData(pNDData)
{
    // add the data array name
    mAttributes["Name"] = rDataArrayName;

    // add the data type information
    if constexpr(std::is_same_v<TDataType, unsigned char>) {
        mAttributes["type"] = "UInt" + std::to_string(sizeof(TDataType) * 8);
    } else if constexpr(std::is_same_v<TDataType, bool>) {
        mAttributes["type"] = "UInt" + std::to_string(sizeof(TDataType) * 8);
    } else if constexpr(std::is_same_v<TDataType, int>) {
        mAttributes["type"] = "Int" + std::to_string(sizeof(TDataType) * 8);
    } else if constexpr(std::is_same_v<TDataType, double>) {
        mAttributes["type"] = "Float" + std::to_string(sizeof(TDataType) * 8);
    } else {
        KRATOS_ERROR << "Unsupported data type.";
    }

    const auto& shape = pNDData->Shape();
    mAttributes["NumberOfComponents"] = std::to_string(std::accumulate(shape.begin() + 1, shape.end(), 1u, std::multiplies<unsigned int>{}));
}

// template instantiations
template class XmlNDDataElement<unsigned char>;
template class XmlNDDataElement<bool>;
template class XmlNDDataElement<int>;
template class XmlNDDataElement<double>;

} // namespace Kratos
