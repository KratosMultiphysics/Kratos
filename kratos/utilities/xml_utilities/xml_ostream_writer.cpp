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
#include <numeric>

// Project includes
#include "input_output/base_64_encoded_output.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"

// Include base h
#include "xml_ostream_writer.h"

namespace Kratos {

XmlOStreamWriter::XmlOStreamWriter(
    std::ostream& rOStream,
    const WriterFormat OutputFormat,
    const IndexType Precision)
    : mrOStream(rOStream),
      mOutputFormat(OutputFormat)
{
    mrOStream << std::scientific << std::setprecision(Precision);
}

void XmlOStreamWriter::WriteElement(
    const std::string& rTagName,
    const std::vector<std::pair<std::string, std::string>>& rAttributes,
    const IndexType Level,
    const bool IsEmptyElement)
{
    WriteAttributes(rTagName, rAttributes, Level);

    if (IsEmptyElement) {
        mrOStream << "/>\n";
    } else {
        mrOStream << ">\n";
    }

}

void XmlOStreamWriter::CloseElement(
    const std::string& rTagName,
    const IndexType Level)
{
    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
    mrOStream << tabbing << "</" << rTagName << ">\n";
}

void XmlOStreamWriter::WriteDataElement(
    const std::string& rTagName,
    const std::vector<std::pair<std::string, std::string>>& rAttributes,
    const std::vector<Expression::Pointer>& rExpressions,
    const IndexType Level)
{

    switch (mOutputFormat){
        case WriterFormat::ASCII:
            if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<char>*>(&*pExpression); })) {
                WriteDataElementAscii<LiteralFlatExpression<char>>(rTagName, rAttributes, Level, rExpressions);
            } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<int>*>(&*pExpression); })) {
                WriteDataElementAscii<LiteralFlatExpression<int>>(rTagName, rAttributes, Level, rExpressions);
            } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<double>*>(&*pExpression); })) {
                WriteDataElementAscii<LiteralFlatExpression<double>>(rTagName, rAttributes, Level, rExpressions);
            } else {
                WriteDataElementAscii<Expression>(rTagName, rAttributes, Level, rExpressions);
            }
            break;
        case WriterFormat::BINARY:
            if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<char>*>(&*pExpression); })) {
                WriteDataElementBinary<LiteralFlatExpression<char>>(rTagName, rAttributes, Level, rExpressions);
            } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<int>*>(&*pExpression); })) {
                WriteDataElementBinary<LiteralFlatExpression<int>>(rTagName, rAttributes, Level, rExpressions);
            } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<const LiteralFlatExpression<double>*>(&*pExpression); })) {
                WriteDataElementBinary<LiteralFlatExpression<double>>(rTagName, rAttributes, Level, rExpressions);
            } else {
                WriteDataElementBinary<Expression>(rTagName, rAttributes, Level, rExpressions);
            }
            break;
    }
}


void XmlOStreamWriter::WriteAttributes(
    const std::string& rTagName,
    const std::vector<std::pair<std::string, std::string>>& rAttributes,
    const IndexType Level)
{
    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
    mrOStream << tabbing << "<" << rTagName;
    if (rAttributes.size() > 0) {
        for (const auto& r_pair : rAttributes) {
            mrOStream << " " << r_pair.first << "=\"" << r_pair.second << "\"";
        }
    }
}

template<class TExpressionType>
void XmlOStreamWriter::WriteDataElementAscii(
    const std::string& rTagName,
    const std::vector<std::pair<std::string, std::string>>& rAttributes,
    const IndexType Level,
    const std::vector<Expression::Pointer>& rExpressions)
{
    WriteAttributes(rTagName, rAttributes, Level);
    // add format
    const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
    mrOStream << " format=\"ascii\">\n" << tabbing;

    std::vector<const TExpressionType*> transformed_expressions(rExpressions.size());
    std::transform(rExpressions.begin(), rExpressions.end(),
                transformed_expressions.begin(), [](auto& pExpression) {
                    return dynamic_cast<const TExpressionType*>(&*(pExpression));
                });

    for (const auto& p_expression : transformed_expressions) {
        auto itr = p_expression->cbegin();
        auto data_end = p_expression->cend();
        for (; itr != data_end; ++itr) {
            if constexpr(std::is_same_v<std::remove_pointer_t<decltype(itr)>, const char>) {
                mrOStream << "  " << static_cast<int>(*(itr));
            } else {
                mrOStream << "  " << *itr;
            }
        }
    }

    mrOStream << "\n";
    CloseElement(rTagName, Level);
}

template<class TExpressionType>
void XmlOStreamWriter::WriteDataElementBinary(
    const std::string& rTagName,
    const std::vector<std::pair<std::string, std::string>>& rAttributes,
    const IndexType Level,
    const std::vector<Expression::Pointer>& rExpressions)
{
    WriteAttributes(rTagName, rAttributes, Level);

    std::vector<const TExpressionType*> transformed_expressions(rExpressions.size());
    std::transform(rExpressions.begin(), rExpressions.end(),
                transformed_expressions.begin(), [](auto& pExpression) {
                    return dynamic_cast<const TExpressionType*>(&*(pExpression));
                });

    if (rExpressions.size() == 0) {
        mrOStream << " format=\"binary\"/>\n";
        return;
    } else {
        const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
        mrOStream << " format=\"binary\">\n" << tabbing << "  ";
    }

    constexpr IndexType data_type_size = sizeof(decltype(*(transformed_expressions[0]->cbegin())));
    const IndexType total_entities = std::accumulate(rExpressions.begin(), rExpressions.end(), 0U, [](const IndexType LHS, const auto& pExpression) { return LHS + pExpression->NumberOfEntities();});
    const IndexType flattened_shape_size = rExpressions[0]->GetItemComponentCount();
    const IndexType total_number_of_values = total_entities * flattened_shape_size;
    const unsigned int total_data_size = total_number_of_values * data_type_size;

    {
        Base64EncodedOutput base64_encoder(mrOStream);

        base64_encoder.WriteData(&total_data_size, 1);
        for (const auto p_expression : transformed_expressions) {
            base64_encoder.WriteData(p_expression->cbegin(), p_expression->NumberOfEntities() * flattened_shape_size);
        }
    }

    mrOStream << "\n";
    CloseElement(rTagName, Level);
}

const std::string XmlOStreamWriter::GetTabbing(const IndexType Level)
{
    std::stringstream ss_tabbing;
    for (IndexType i = 0; i < Level; ++i) {
        ss_tabbing << "   ";
    }
    return ss_tabbing.str();
}

}