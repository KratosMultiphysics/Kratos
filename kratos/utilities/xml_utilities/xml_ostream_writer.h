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

/* @class XmlOStreamWriter
 * @ingroup KratosCore
 * @brief Output stream writer for XML format.
 * @author Suneth Warnakulasuriya
 */
class KRATOS_API(KRATOS_CORE) XmlOStreamWriter
{
public:
    ///@name Life cycle
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Public enums
    ///@{

    /// Enumerations for the output writer format.
    enum WriterFormat
    {
        ASCII,  /// ASCII format.
        BINARY  /// Binary format.
    };

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Constructor.
     * @param rOStream The output stream to write to.
     * @param OutputFormat The format of the writer.
     * @param Precision The precision for floating-point numbers.
     */
    XmlOStreamWriter(
        std::ostream& rOStream,
        const WriterFormat OutputFormat,
        const IndexType Precision);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Writes an XML element with attributes.
     * @param rTagName The tag name of the element.
     * @param rAttributes The attributes of the element.
     * @param Level The indentation level.
     * @param IsEmptyElement Flag indicating if the element is empty.
     */
    void WriteElement(
        const std::string& rTagName,
        const std::vector<std::pair<std::string, std::string>>& rAttributes,
        const IndexType Level,
        const bool IsEmptyElement);

    /**
     * @brief Closes an XML element.
     * @param rTagName The tag name of the element to close.
     * @param Level The indentation level.
     */
    void CloseElement(
        const std::string& rTagName,
        const IndexType Level);

    /**
     * @brief Writes an XML data element with attributes and expressions.
     * @param rTagName The tag name of the data element.
     * @param rAttributes The attributes of the data element.
     * @param rExpressions The expressions to write as data.
     * @param Level The indentation level.
     */
    void WriteDataElement(
        const std::string& rTagName,
        const std::vector<std::pair<std::string, std::string>>& rAttributes,
        const std::vector<Expression::Pointer>& rExpressions,
        const IndexType Level);

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::ostream& mrOStream; /// The output stream

    const WriterFormat mOutputFormat; /// The output format.

    ///@}
    ///@name Private operations
    ///@{

    /**
     * @brief Writes the attributes of an XML element.
     * @param rTagName The tag name of the element.
     * @param rAttributes The attributes of the element.
     * @param Level The indentation level.
     */
    void WriteAttributes(
        const std::string& rTagName,
        const std::vector<std::pair<std::string, std::string>>& rAttributes,
        const IndexType Level);

    /**
     * @brief Writes an ASCII data element.
     * @tparam TExpressionType The type of expression to write.
     * @param rTagName The tag name of the data element.
     * @param rAttributes The attributes of the data element.
     * @param Level The indentation level.
     * @param rExpressions The expressions to write as data.
     */
    template <class TExpressionType>
    void WriteDataElementAscii(
        const std::string& rTagName,
        const std::vector<std::pair<std::string, std::string>>& rAttributes,
        const IndexType Level,
        const std::vector<Expression::Pointer>& rExpressions);

    /**
     * @brief Writes a binary data element.
     * @tparam TExpressionType The type of expression to write.
     * @param rTagName The tag name of the data element.
     * @param rAttributes The attributes of the data element.
     * @param Level The indentation level.
     * @param rExpressions The expressions to write as data.
     */
    template <class TExpressionType>
    void WriteDataElementBinary(
        const std::string& rTagName,
        const std::vector<std::pair<std::string, std::string>>& rAttributes,
        const IndexType Level,
        const std::vector<Expression::Pointer>& rExpressions);

    /**
     * @brief Gets the indentation string based on the level.
     * @param Level The indentation level.
     * @return The indentation string.
     */
    static const std::string GetTabbing(const IndexType Level);

    ///@}
};

} // namespace Kratos