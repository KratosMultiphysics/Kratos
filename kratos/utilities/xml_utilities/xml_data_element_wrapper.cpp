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

// Project includes
#include "xml_elements_array.h"
#include "xml_utils.h"

// Include base h
#include "xml_data_element_wrapper.h"

namespace Kratos {

XmlDataElementWrapper::XmlDataElementWrapper(const std::string& rTagName)
    : XmlElement(rTagName)
{
}

XmlElement::Pointer XmlDataElementWrapper::Get(
    const std::string& rDataArrayName,
    NDDataPointerType pNDData)
{
    KRATOS_TRY

    auto p_element = Kratos::make_shared<XmlElementsArray>("DataArray");
    p_element->AddAttribute("Name", rDataArrayName);
    std::visit([&p_element](auto p_nd_data) {
        XmlUtilities::AddDataArrayAttributes(*p_element, *p_nd_data);
    }, pNDData);

    return p_element;

    KRATOS_CATCH("");
}

} // namespace Kratos