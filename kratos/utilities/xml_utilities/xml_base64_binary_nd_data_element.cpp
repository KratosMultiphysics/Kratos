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
#include <cstdint>

// Project includes
#include "input_output/base_64_encoded_output.h"

// Include base h
#include "xml_base64_binary_nd_data_element.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TDataType>
XmlBase64BinaryNDDataElement<TDataType>::XmlBase64BinaryNDDataElement(
    const std::string& rDataArrayName,
    typename NDData<TDataType>::Pointer pNDData)
    : BaseType(rDataArrayName, pNDData)
{
    // modify the format information
    this->AddAttribute("format", "binary");
}

template<class TDataType>
void XmlBase64BinaryNDDataElement<TDataType>::Write(
    std::ostream& rOStream,
    const IndexType Level) const
{
    KRATOS_TRY

    const auto span = this->mpNDData->ViewData();

    if (span.size() == 0) {
        this->WriteEmptyElementTag(rOStream, Level);
    } else {
        this->WriteElementTagStart(rOStream, Level);

        // write the data
        const std::string tabbing(Level * 3, ' ');

        rOStream << tabbing << "  ";

        const std::uint64_t total_data_size = sizeof(TDataType) * span.size();

        {
            // Base 64 encoded output should be
            // on a separate scope, because destructor will write the padding
            Base64EncodedOutput base64_encoder(rOStream);
            base64_encoder.WriteData(&total_data_size, 1);
            base64_encoder.WriteData(span.begin(), span.size());
        }

        rOStream << "\n";

        this->WriteElementTagEnd(rOStream, Level);
    }

    KRATOS_CATCH("");
}

// template instantiations
template class XmlBase64BinaryNDDataElement<unsigned char>;
template class XmlBase64BinaryNDDataElement<bool>;
template class XmlBase64BinaryNDDataElement<int>;
template class XmlBase64BinaryNDDataElement<double>;

} // namespace Kratos
