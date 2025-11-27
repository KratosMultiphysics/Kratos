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
#include "xml_element.h"

namespace Kratos::Future {

///@name Kratos Classes
///@{

/* @class XmlElementsArray
 * @ingroup KratosCore
 * @brief Output stream writer for XML format.
 * @author Suneth Warnakulasuriya
 */
class KRATOS_API(KRATOS_CORE) XmlElementsArray : public XmlElement
{
public:
    ///@name Life cycle
    ///@{

    using IndexType = std::size_t;

    using BaseType = XmlElement;

    KRATOS_CLASS_POINTER_DEFINITION(XmlElementsArray);

    ///@}
    ///@name Life cycle
    ///@{

    XmlElementsArray(const std::string& rTagName);

    ///@}
    ///@name Public operations
    ///@{

    void AddElement(XmlElement::Pointer pElement);

    std::vector<XmlElement::Pointer> GetElements() const;

    void ClearElements();

    /**
     * @brief Writes an XML expression element.
     * @param XmlExpressionElement Expression xml element to be written.
     * @param Level The indentation level.
     */
    void Write(
        std::ostream& rOStream,
        const IndexType Level = 0) const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::vector<XmlElement::Pointer> mElementsArray;

    ///@}
};

} // namespace Kratos::Future