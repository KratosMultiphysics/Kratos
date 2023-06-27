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
#include <string>
#include <vector>
#include <tuple>

// Project includes
#include "expression/literal_flat_expression.h"

// Include base h
#include "xml_ostream_writer.h"

namespace Kratos {

XmlOStreamWriter::XmlOStreamWriter(std::ostream& rOStream)
    : mrOStream(rOStream)
{
}

void XmlOStreamWriter::WriteElement(
    const XmlExpressionElement& rElement,
    const IndexType Level)
{
    const std::string tabbing(Level * 3, ' ');

    mrOStream << tabbing << "<" << rElement.GetTagName();
    const auto& r_attributes = rElement.GetAttributes();
    if (r_attributes.size() > 0) {
        for (const auto& r_pair : r_attributes) {
            mrOStream << " " << r_pair.first << "=\"" << r_pair.second << "\"";
        }
    }

    const auto& r_sub_elements = rElement.GetElements();
    const auto& r_expressions = rElement.GetExpressions();

    if (r_sub_elements.size() > 0) {
        // write sub elements
        mrOStream << ">\n";
        for (const auto& p_sub_element : r_sub_elements) {
            this->WriteElement(*p_sub_element, Level + 1);
        }
        mrOStream << tabbing << "</" << rElement.GetTagName() << ">\n";
    } else if (r_expressions.size() > 0) {
        // write a data expression element
        this->WriteExpressions(r_expressions, tabbing);
        // close the element
        mrOStream << "\n" << tabbing << "</" << rElement.GetTagName() << ">\n";
    } else {
        // then it is an empty element.
        mrOStream << "/>\n";
    }
}

} // namespace Kratos