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
#include "xml_element.h"

namespace Kratos {

XmlElement::XmlElement(const std::string& rTagName)
    : mTagName(rTagName)
{
}

void XmlElement::WriteElementTagStart(
    std::ostream& rOStream,
    const IndexType Level) const
{
    const std::string tabbing(Level * 3, ' ');

    rOStream << tabbing << "<" << mTagName;

    for (const auto& r_pair : mAttributes) {
        rOStream << " " << r_pair.first << "=\"" << r_pair.second << "\"";
    }

    rOStream << ">\n";
}

void XmlElement::WriteElementTagEnd(
    std::ostream& rOStream,
    const IndexType Level) const
{
    const std::string tabbing(Level * 3, ' ');
    rOStream << tabbing << "</" << mTagName << ">\n";
}
void XmlElement::WriteEmptyElementTag(
    std::ostream& rOStream,
    const IndexType Level) const
{
    const std::string tabbing(Level * 3, ' ');

    rOStream << tabbing << "<" << mTagName;

    for (const auto& r_pair : mAttributes) {
        rOStream << " " << r_pair.first << "=\"" << r_pair.second << "\"";
    }

    rOStream << "/>\n";
}

} // namespace Kratos