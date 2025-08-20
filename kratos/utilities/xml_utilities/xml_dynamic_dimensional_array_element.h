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
#include <variant>

// Project includes
#include "includes/define.h"
#include "containers/dynamic_dimensional_array.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) XmlDynamicDimensionalArrayElement {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using ArrayPointerType = std::variant<
                                    DynamicDimensionalArray<unsigned char>::Pointer,
                                    DynamicDimensionalArray<bool>::Pointer,
                                    DynamicDimensionalArray<int>::Pointer,
                                    DynamicDimensionalArray<double>::Pointer
                                >;

    KRATOS_CLASS_POINTER_DEFINITION(XmlDynamicDimensionalArrayElement);

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Constructor.
     * @param rTagName The tag name of the XML element.
     */
    XmlDynamicDimensionalArrayElement(const std::string& rTagName);

    /**
     * @brief Constructor.
     * @param rDataName The name of the data element.
     * @param rDynamicDimensionalArrays The expressions to write as data.
     */
    XmlDynamicDimensionalArrayElement(
        const std::string& rDataName,
        const std::vector<ArrayPointerType>& rDynamicDimensionalArrays);

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
    void AddElement(const XmlDynamicDimensionalArrayElement::Pointer pXmlElement);

    /**
     * @brief Get sub-elements with a specific tag name.
     * @param rTagName The tag name of the sub-elements.
     * @return The vector of sub-elements.
     */
    std::vector<XmlDynamicDimensionalArrayElement::Pointer> GetElements(const std::string& rTagName) const;

    /**
     * @brief Get all sub-elements of the XML element.
     * @return The vector of sub-elements.
     */
    const std::vector<XmlDynamicDimensionalArrayElement::Pointer>& GetElements() const;

    /**
     * @brief Clear all sub-elements of the XML element.
     */
    void ClearElements();

    const std::vector<ArrayPointerType> GetDynamicDimensionalArrays() const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const std::string mTagName; /// The tag name of the XML element.

    std::vector<std::pair<std::string, std::string>> mAttributes; /// The attributes of the XML element.

    std::vector<XmlDynamicDimensionalArrayElement::Pointer> mXmlElements; /// The sub-elements of the XML element.

    const std::vector<ArrayPointerType> mDynamicDimensionalArrays; /// The arrays to write as data.

    ///@}
};

///@}

} // namespace Kratos
