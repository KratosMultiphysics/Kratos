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
#include "expression/literal_flat_expression.h"
#include "utilities/xml_utilities/xml_expression_element.h"

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
    ///@name Life cycle
    ///@{

    /**
     * @brief Constructor.
     * @param rOStream The output stream to write to.
     */
    XmlOStreamWriter(std::ostream& rOStream);

    virtual ~XmlOStreamWriter() = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Writes an XML expression element.
     * @param XmlExpressionElement Expression xml element to be written.
     * @param Level The indentation level.
     */
    void WriteElement(
        const XmlExpressionElement& rElement,
        const IndexType Level = 0);

    ///@}

protected:
    ///@name Protected member variables
    ///@{

    std::ostream& mrOStream; /// The output stream

    ///@}
    ///@name Protected operations
    ///@{

    /**
     * @brief Writes generic lazy type expressions
     *
     * @param rExpressions      Expressions list to write.
     * @param rTabbing          Tabbing used for expression writing.
     */
    virtual void WriteExpressions(
        const std::vector<Expression::ConstPointer>& rExpressions,
        const std::string& rTabbing) = 0;

    ///@}
};

} // namespace Kratos