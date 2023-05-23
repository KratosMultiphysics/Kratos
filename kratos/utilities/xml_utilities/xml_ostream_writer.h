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

#pragma once

// System includes
#include <string>
#include <vector>

// Project includes
#include "includes/define.h"
#include "containers/container_expression/expressions/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) XmlOStreamWriter
{
public:
    ///@name Life cycle
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Public enums
    ///@{

    enum WriterFormat
    {
        ASCII,
        BINARY
    };

    ///@}
    ///@name Life cycle
    ///@{

    XmlOStreamWriter(
        std::ostream& rOStream,
        const WriterFormat OutputFormat,
        const IndexType Precision);

    ///@}
    ///@name Public operations
    ///@{

    void WriteElement(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const IndexType Level,
        const bool IsEmptyElement);

    void CloseElement(
        const std::string& rTagName,
        const IndexType Level);

    void WriteDataElement(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const std::vector<Expression::Pointer>& rExpressions,
        const IndexType Level);

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::ostream& mrOStream;

    const WriterFormat mOutputFormat;

    ///@}
    ///@name Private operations
    ///@{

    void WriteAttributes(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const IndexType Level);

    template <class TExpressionType>
    void WriteDataElementAscii(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const IndexType Level,
        const std::vector<Expression::Pointer>& rExpressions);

    template <class TExpressionType>
    void WriteDataElementBinary(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const IndexType Level,
        const std::vector<Expression::Pointer>& rExpressions);

    static const std::string GetTabbing(const IndexType Level);

    ///@}
};

} // namespace Kratos