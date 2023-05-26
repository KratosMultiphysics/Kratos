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
#include <ostream>

// Project includes
#include "includes/define.h"
#include "xml_expression_element.h"

namespace Kratos {

///@name Kratos Classes
///@{

/* @class XmlOStreamAsciiWriter
 * @ingroup KratosCore
 * @brief Output stream ascii writer for XML format.
 * @author Suneth Warnakulasuriya
 */
class KRATOS_API(KRATOS_CORE) XmlOStreamAsciiWriter
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
     * @param Precision The precision for floating-point numbers.
     */
    XmlOStreamAsciiWriter(
        std::ostream& rOStream,
        const IndexType Precision);

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

private:
    ///@name Private member variables
    ///@{

    std::ostream& mrOStream; /// The output stream

    ///@}
};

} // namespace Kratos