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
#include <numeric>

// Project includes
#include "input_output/base_64_encoded_output.h"
#include "expression/literal_flat_expression.h"

// Include base h
#include "xml_ostream_base64_binary_writer.h"

namespace Kratos {

namespace XmlOStreamBase64BinaryWriterHelperUtilities
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

    constexpr std::size_t data_type_size = sizeof(decltype(*(transformed_expressions[0]->cbegin())));
    const std::size_t total_entities = std::accumulate(rExpressions.begin(), rExpressions.end(), 0U, [](const std::size_t LHS, const auto& pExpression) { return LHS + pExpression->NumberOfEntities();});
    const std::size_t flattened_shape_size = rExpressions[0]->GetItemComponentCount();
    const std::size_t total_number_of_values = total_entities * flattened_shape_size;
    const unsigned int total_data_size = total_number_of_values * data_type_size;

    {
        Base64EncodedOutput base64_encoder(rOStream);

        base64_encoder.WriteData(&total_data_size, 1);
        for (const auto p_expression : transformed_expressions) {
            base64_encoder.WriteData(p_expression->cbegin(), p_expression->NumberOfEntities() * flattened_shape_size);
        }
    }
}
} // namespace XmlOStreamBase64BinaryWriterHelperUtilities

XmlOStreamBase64BinaryWriter::XmlOStreamBase64BinaryWriter(std::ostream& rOStream)
    : XmlOStreamWriter(rOStream)
{
}

void XmlOStreamBase64BinaryWriter::WriteExpressions(
    const std::vector<Expression::ConstPointer>& rExpressions,
    const std::string& rTabbing)
{
    mrOStream << " format=\"binary\">\n" << rTabbing << "  ";

    if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<char>*>(&*pExpression); })) {
        XmlOStreamBase64BinaryWriterHelperUtilities::WriteExpression<LiteralFlatExpression<char>>(mrOStream, rExpressions);
    } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<int>*>(&*pExpression); })) {
        XmlOStreamBase64BinaryWriterHelperUtilities::WriteExpression<LiteralFlatExpression<int>>(mrOStream, rExpressions);
    } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<double>*>(&*pExpression); })) {
        XmlOStreamBase64BinaryWriterHelperUtilities::WriteExpression<LiteralFlatExpression<double>>(mrOStream, rExpressions);
    } else {
        XmlOStreamBase64BinaryWriterHelperUtilities::WriteExpression<Expression>(mrOStream, rExpressions);
    }
}

} // namespace Kratos