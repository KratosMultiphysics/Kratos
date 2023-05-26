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
#include <iomanip>

// Project includes
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"

// Include base h
#include "xml_ostream_ascii_writer.h"

namespace Kratos {

namespace XmlOStreamAsciiWriterHelperUtilities
{
template <class TExpressionType>
void WriteExpression(
    std::ostream& rOStream,
    const std::vector<Expression::Pointer>& rExpressions)
{
    std::vector<const TExpressionType*> transformed_expressions(rExpressions.size());
    std::transform(rExpressions.begin(), rExpressions.end(),
                   transformed_expressions.begin(), [](auto& pExpression) {
                       return dynamic_cast<const TExpressionType*>(&*(pExpression));
                   });

    for (const auto& p_expression : transformed_expressions) {
        auto itr = p_expression->cbegin();
        auto data_end = p_expression->cend();
        for (; itr != data_end; ++itr) {
            if constexpr (std::is_same_v<std::remove_pointer_t<decltype(itr)>, const char>) {
                rOStream << "  " << static_cast<int>(*(itr));
            }
            else {
                rOStream << "  " << *itr;
            }
        }
    }
}
} // namespace XmlOStreamAsciiWriterHelperUtilities

XmlOStreamAsciiWriter::XmlOStreamAsciiWriter(
    std::ostream& rOStream,
    const IndexType Precision)
    : mrOStream(rOStream)
{
    mrOStream << std::scientific << std::setprecision(Precision);
}

void XmlOStreamAsciiWriter::WriteElement(
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
        mrOStream << " format=\"ascii\">\n" << tabbing;

        if (std::all_of(r_expressions.begin(), r_expressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<char>*>(&*pExpression); })) {
            XmlOStreamAsciiWriterHelperUtilities::WriteExpression<LiteralFlatExpression<char>>(mrOStream, r_expressions);
        } else if (std::all_of(r_expressions.begin(), r_expressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<int>*>(&*pExpression); })) {
            XmlOStreamAsciiWriterHelperUtilities::WriteExpression<LiteralFlatExpression<int>>(mrOStream, r_expressions);
        } else if (std::all_of(r_expressions.begin(), r_expressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<double>*>(&*pExpression); })) {
            XmlOStreamAsciiWriterHelperUtilities::WriteExpression<LiteralFlatExpression<double>>(mrOStream, r_expressions);
        } else {
            XmlOStreamAsciiWriterHelperUtilities::WriteExpression<Expression>(mrOStream, r_expressions);
        }
        mrOStream << "\n" << tabbing << "</" << rElement.GetTagName() << ">\n";
    } else {
        // then it is an empty element.
        mrOStream << "/>\n";
    }
}

} // namespace Kratos