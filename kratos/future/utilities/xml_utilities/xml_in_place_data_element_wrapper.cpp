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

// External includes

// Project includes
#include "utilities/data_type_traits.h"
#include "xml_ascii_nd_data_element.h"
#include "xml_base64_binary_nd_data_element.h"
#include "xml_utils.h"

// Include base h
#include "xml_in_place_data_element_wrapper.h"

namespace Kratos::Future {

XmlInPlaceDataElementWrapper::XmlInPlaceDataElementWrapper(
    const XmlOutputType OutputType,
    const IndexType Precision)
    : BaseType("InPlace"),
      mOutputType(OutputType),
      mPrecision(Precision)
{
}

XmlElement::Pointer XmlInPlaceDataElementWrapper::Get(
    const std::string& rDataArrayName,
    NDDataPointerType pNDData)
{
    KRATOS_TRY

    return std::visit([this, &rDataArrayName](auto p_nd_data) -> XmlElement::Pointer {
        const auto& shape = p_nd_data->Shape();

        KRATOS_ERROR_IF(std::accumulate(shape.begin() + 1, shape.end(), 1u, std::multiplies<unsigned int>{}) == 0)
            << "Writing data arrays with zero components is prohibited. [ data array name = \""
            << rDataArrayName << "\", shape = " << shape << ", nd data = " << *p_nd_data << " ].\n";

        using nd_data_data_type = typename BareType<decltype(*p_nd_data)>::DataType;

        switch (this->mOutputType) {
            case ASCII:
                return Kratos::make_shared<XmlAsciiNDDataElement<nd_data_data_type>>(rDataArrayName, p_nd_data, this->mPrecision);
            case BINARY:
                return Kratos::make_shared<XmlBase64BinaryNDDataElement<nd_data_data_type>>(rDataArrayName, p_nd_data);
            default:
                return nullptr;
        }
    }, pNDData);

    KRATOS_CATCH("");
}

void XmlInPlaceDataElementWrapper::Write(
    std::ostream& rOStream,
    const IndexType Level) const
{
    KRATOS_TRY

    KRATOS_ERROR << "This is dummy wrapper element only.";

    KRATOS_CATCH("");
}

} // namespace Kratos::Future