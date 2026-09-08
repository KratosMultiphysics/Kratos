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
#include <map>

// Project includes
#include "includes/define.h"

namespace Kratos {

///@name Kratos Classes
///@{

/* @class XmlElement
 * @ingroup KratosCore
 * @brief Output stream writer for XML format.
 * @author Suneth Warnakulasuriya
 */
class KRATOS_API(KRATOS_CORE) XmlElement
{
public:
    ///@name Life cycle
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(XmlElement);

    ///@}
    ///@name Life cycle
    ///@{

    XmlElement(const std::string& rTagName);

    virtual ~XmlElement() = default;

    ///@}
    ///@name Public operations
    ///@{

    void AddAttribute(
        const std::string& rName,
        const std::string& rValue);

    std::string GetTagName() const;

    std::map<std::string, std::string> GetAttributes() const;

    void ClearAttributes();

    /**
     * @brief Writes an XML expression element.
     * @param XmlExpressionElement Expression xml element to be written.
     * @param Level The indentation level.
     */
    virtual void Write(
        std::ostream& rOStream,
        const IndexType Level = 0) const = 0;

    ///@}

protected:
    ///@name Protected operations
    ///@{

    void WriteElementTagStart(
        std::ostream& rOStream,
        const IndexType Level) const;

    void WriteElementTagEnd(
        std::ostream& rOStream,
        const IndexType Level) const;

    void WriteEmptyElementTag(
        std::ostream& rOStream,
        const IndexType Level) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::map<std::string, std::string> mAttributes;

    const std::string mTagName;

    ///@}
};

} // namespace Kratos