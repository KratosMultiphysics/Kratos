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
#include <type_traits>

// Project includes

// Include base h
#include "xml_ascii_nd_data_element.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TDataType>
XmlAsciiNDDataElement<TDataType>::XmlAsciiNDDataElement(
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

    // add the shape information
    const auto shape = mpNDData->Shape();
    mAttributes["NumberOfComponents"] = std::to_string(std::accumulate(shape.begin() + 1, shape.end(), 1u, std::multiplies<unsigned int>{}));

    // add the format information
    mAttributes["format"] = "ascii";
}

template<class TDataType>
void XmlAsciiNDDataElement<TDataType>::Write(
    std::ostream& rOStream,
    const IndexType Level) const
{
    const auto span = mpNDData->ViewData();

    if (span.size() == 0) {
        WriteEmptyElementTag(rOStream, Level);
    } else {
        WriteElementTagStart(rOStream, Level);

        // write the data
        const std::string tabbing(Level * 3, ' ');

        rOStream << tabbing;
        for (auto itr = span.begin(); itr != span.end(); ++itr) {
            if constexpr(std::is_same_v<TDataType, unsigned char>) {
                rOStream << "  " << static_cast<int>(*itr);
            } else {
                rOStream << "  " << *itr;
            }
        }

        WriteElementTagEnd(rOStream, Level);
    }
}

// template instantiations
template class XmlAsciiNDDataElement<unsigned char>;
template class XmlAsciiNDDataElement<bool>;
template class XmlAsciiNDDataElement<int>;
template class XmlAsciiNDDataElement<double>;

} // namespace Kratos
