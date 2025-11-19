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

// Include base h
#include "xml_elements_array.h"

namespace Kratos::Future {

XmlElementsArray::XmlElementsArray(const std::string& rTagName)
    : BaseType(rTagName)
{
}

void XmlElementsArray::AddElement(XmlElement::Pointer pElement)
{
    mElementsArray.push_back(pElement);
}

std::vector<XmlElement::Pointer> XmlElementsArray::GetElements() const
{
    return mElementsArray;
}

void XmlElementsArray::ClearElements()
{
    mElementsArray.clear();
}

void XmlElementsArray:: Write(
    std::ostream& rOStream,
    const IndexType Level) const
{
    KRATOS_TRY

    if (mElementsArray.empty()) {
        WriteEmptyElementTag(rOStream, Level);
    } else {
        WriteElementTagStart(rOStream, Level);

        for (const auto& p_xml_element : mElementsArray) {
            p_xml_element->Write(rOStream, Level + 1);
        }

        WriteElementTagEnd(rOStream, Level);
    }

    KRATOS_CATCH("");
}

} // namespace Kratos::Future