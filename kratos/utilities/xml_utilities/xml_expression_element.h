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
#include "expression/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) XmlExpressionElement {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(XmlExpressionElement);

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Constructor.
     * @param rTagName The tag name of the XML element.
     */
    XmlExpressionElement(const std::string& rTagName);

    /**
     * @brief Constructor.
     * @param rDataName The name of the data element.
     * @param rExpressions The expressions to write as data.
     */
    XmlExpressionElement(
        const std::string& rDataName,
        const std::vector<Expression::ConstPointer>& rExpressions);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Get the tag name of the XML element.
     * @return The tag name.
     */
    const std::string GetTagName() const;

    /**
     * @brief Add an attribute to the XML element.
     * @param rName The name of the attribute.
     * @param rValue The value of the attribute.
     */
    void AddAttribute(
        const std::string& rName,
        const std::string& rValue);

    /**
     * @brief Get the attributes of the XML element.
     * @return The attributes.
     */
    const std::vector<std::pair<std::string, std::string>>& GetAttributes() const;

    /**
     * @brief Clear the attributes of the XML element.
     */
    void ClearAttributes();

    /**
     * @brief Add a sub-element to the XML element.
     * @param pXmlElement The sub-element to add.
     */
    void AddElement(const XmlExpressionElement::Pointer pXmlElement);

    /**
     * @brief Get sub-elements with a specific tag name.
     * @param rTagName The tag name of the sub-elements.
     * @return The vector of sub-elements.
     */
    std::vector<XmlExpressionElement::Pointer> GetElements(const std::string& rTagName) const;

    /**
     * @brief Get all sub-elements of the XML element.
     * @return The vector of sub-elements.
     */
    const std::vector<XmlExpressionElement::Pointer>& GetElements() const;

    /**
     * @brief Clear all sub-elements of the XML element.
     */
    void ClearElements();

    const std::vector<Expression::ConstPointer> GetExpressions() const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const std::string mTagName; /// The tag name of the XML element.

    std::vector<std::pair<std::string, std::string>> mAttributes; /// The attributes of the XML element.

    std::vector<XmlExpressionElement::Pointer> mXmlElements; /// The sub-elements of the XML element.

    const std::vector<Expression::ConstPointer> mExpressions; /// The expressions to write as data.

    ///@}
};

///@}

} // namespace Kratos
