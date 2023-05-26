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
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"

// Include base h
#include "xml_ostream_base64_binary_writer.h"

namespace Kratos {

namespace XmlOStreamBase64BinaryWriterHelperUtilities
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
    : mrOStream(rOStream)
{
}

void XmlOStreamBase64BinaryWriter::WriteElement(
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
        mrOStream << " format=\"binary\">\n" << tabbing << "  ";

        if (std::all_of(r_expressions.begin(), r_expressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<char>*>(&*pExpression); })) {
            XmlOStreamBase64BinaryWriterHelperUtilities::WriteExpression<LiteralFlatExpression<char>>(mrOStream, r_expressions);
        } else if (std::all_of(r_expressions.begin(), r_expressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<int>*>(&*pExpression); })) {
            XmlOStreamBase64BinaryWriterHelperUtilities::WriteExpression<LiteralFlatExpression<int>>(mrOStream, r_expressions);
        } else if (std::all_of(r_expressions.begin(), r_expressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<double>*>(&*pExpression); })) {
            XmlOStreamBase64BinaryWriterHelperUtilities::WriteExpression<LiteralFlatExpression<double>>(mrOStream, r_expressions);
        } else {
            XmlOStreamBase64BinaryWriterHelperUtilities::WriteExpression<Expression>(mrOStream, r_expressions);
        }
        mrOStream << "\n" << tabbing << "</" << rElement.GetTagName() << ">\n";
    } else {
        // then it is an empty element.
        mrOStream << "/>\n";
    }
}

} // namespace Kratos