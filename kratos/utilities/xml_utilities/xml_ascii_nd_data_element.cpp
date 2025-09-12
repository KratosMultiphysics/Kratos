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
#include "xml_ascii_nd_data_element.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TDataType>
XmlAsciiNDDataElement<TDataType>::XmlAsciiNDDataElement(
    const std::string& rDataArrayName,
    typename NDData<TDataType>::Pointer pNDData,
    const IndexType Precision)
    : BaseType(rDataArrayName, pNDData),
      mPrecision(Precision)
{
    // add the format information
    this->AddAttribute("format", "ascii");
}

template<class TDataType>
void XmlAsciiNDDataElement<TDataType>::Write(
    std::ostream& rOStream,
    const IndexType Level) const
{
    KRATOS_TRY

    rOStream << std::scientific << std::setprecision(mPrecision);

    const auto span = this->mpNDData->ViewData();

    if (span.size() == 0) {
        this->WriteEmptyElementTag(rOStream, Level);
    } else {
        this->WriteElementTagStart(rOStream, Level);

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

        rOStream << "\n";

        this->WriteElementTagEnd(rOStream, Level);
    }

    KRATOS_CATCH("");
}

// template instantiations
template class XmlAsciiNDDataElement<unsigned char>;
template class XmlAsciiNDDataElement<bool>;
template class XmlAsciiNDDataElement<int>;
template class XmlAsciiNDDataElement<double>;

} // namespace Kratos
