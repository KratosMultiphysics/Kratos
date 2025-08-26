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

namespace Kratos {

XmlElementsArray::XmlElementsArray(const std::string& rTagName)
    : BaseType(rTagName)
{
}

void XmlElementsArray::AddElement(XmlElement::Pointer pElement)
{
    mElementsArray.push_back(pElement);
}

void XmlElementsArray::AddAttribute(
    const std::string& rName,
    const std::string& rValue)
{
    mAttributes[rName] = rValue;
}

void XmlElementsArray:: Write(
    std::ostream& rOStream,
    const IndexType Level) const
{
    if (mElementsArray.empty()) {
        WriteEmptyElementTag(rOStream, Level);
    } else {
        WriteElementTagStart(rOStream, Level);

        for (const auto& p_xml_element : mElementsArray) {
            p_xml_element->Write(rOStream, Level + 1);
        }

        WriteElementTagEnd(rOStream, Level);
    }
}

} // namespace Kratos