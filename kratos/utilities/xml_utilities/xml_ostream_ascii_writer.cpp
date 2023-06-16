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
#include "expression/literal_flat_expression.h"

// Include base h
#include "xml_ostream_ascii_writer.h"

namespace Kratos {

namespace XmlOStreamAsciiWriterHelperUtilities
{
template <class TExpressionType>
void WriteExpression(
    std::ostream& rOStream,
    const std::vector<Expression::ConstPointer>& rExpressions)
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
    : XmlOStreamWriter(rOStream)
{
    mrOStream << std::scientific << std::setprecision(Precision);
}

void XmlOStreamAsciiWriter::WriteExpressions(
    const std::vector<Expression::ConstPointer>& rExpressions,
    const std::string& rTabbing)
{
    mrOStream << " format=\"ascii\">\n" << rTabbing;

    if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<char>*>(&*pExpression); })) {
        XmlOStreamAsciiWriterHelperUtilities::WriteExpression<LiteralFlatExpression<char>>(mrOStream, rExpressions);
    } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<int>*>(&*pExpression); })) {
        XmlOStreamAsciiWriterHelperUtilities::WriteExpression<LiteralFlatExpression<int>>(mrOStream, rExpressions);
    } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<double>*>(&*pExpression); })) {
        XmlOStreamAsciiWriterHelperUtilities::WriteExpression<LiteralFlatExpression<double>>(mrOStream, rExpressions);
    } else {
        XmlOStreamAsciiWriterHelperUtilities::WriteExpression<Expression>(mrOStream, rExpressions);
    }
}

} // namespace Kratos