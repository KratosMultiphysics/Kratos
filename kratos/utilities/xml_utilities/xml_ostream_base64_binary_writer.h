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
#include "xml_ostream_writer.h"

namespace Kratos {

///@name Kratos Classes
///@{

/* @class XmlOStreamBase64BinaryWriter
 * @ingroup KratosCore
 * @brief Output stream ascii writer for XML format.
 * @author Suneth Warnakulasuriya
 */
class KRATOS_API(KRATOS_CORE) XmlOStreamBase64BinaryWriter: public XmlOStreamWriter
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
    XmlOStreamBase64BinaryWriter(std::ostream& rOStream);

    ///@}

protected:
    ///@name Protected operations
    ///@{

    void WriteExpressions(
        const std::vector<Expression::ConstPointer>& rExpressions,
        const std::string& rTabbing) override;

    ///@}
};

} // namespace Kratos