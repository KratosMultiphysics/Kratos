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
#include <numeric>

// Project includes
#include "utilities/data_type_traits.h"

// Include base h
#include "xml_dynamic_dimensional_array_element.h"

namespace Kratos {

XmlDynamicDimensionalArrayElement::XmlDynamicDimensionalArrayElement(const std::string& rTagName)
    : mTagName(rTagName)
{
}

XmlDynamicDimensionalArrayElement::XmlDynamicDimensionalArrayElement(
    const std::string& rDataName,
    const std::vector<ArrayPointerType>& rDynamicDimensionalArrays)
    : mTagName("DataArray"),
    mDynamicDimensionalArrays(rDynamicDimensionalArrays)
{
    KRATOS_ERROR_IF(mDynamicDimensionalArrays.size() == 0)
        << "Empty dynamic dimensional array lists are not allowed.";

    IndexType number_of_components = 0;
    std::string type = "";
    for (auto& p_variant_array : mDynamicDimensionalArrays) {
        std::visit([&number_of_components, &type](auto pArray) {
            if (number_of_components == 0) {
                number_of_components = pArray->Size() / pArray->Shape()[0];
            }

            KRATOS_ERROR_IF(pArray->Shape()[0] == 0 || number_of_components != pArray->Size() / pArray->Shape()[0])
                << "Found dynamic dimensional array with zero elements or with mismatching shapes. "
                << " [ dynamic dimensional array shape = " << pArray->Shape() << ", number of components = "
                << number_of_components << " ].\n";

            if constexpr(std::is_same_v<typename BareType<decltype(*pArray)>::DataType, unsigned char>) {
                if (type == "") type = "UInt" + std::to_string(sizeof(char) * 8);
                else if (type != "UInt") type = "Float" + std::to_string(sizeof(double) * 8);
            } else if constexpr(std::is_same_v<typename BareType<decltype(*pArray)>::DataType, bool>) {
                if (type == "") type = "UInt" + std::to_string(sizeof(char) * 8);
                else if (type != "UInt") type = "Float" + std::to_string(sizeof(double) * 8);
            } else if constexpr(std::is_same_v<typename BareType<decltype(*pArray)>::DataType, int>) {
                if (type == "") type = "Int" + std::to_string(sizeof(int) * 8);
                else if (type != "Int") type = "Float" + std::to_string(sizeof(double) * 8);
            } else if constexpr(std::is_same_v<typename BareType<decltype(*pArray)>::DataType, double>) {
                if (type == "") type = "Float" + std::to_string(sizeof(double) * 8);
                else if (type != "Float") type = "Float" + std::to_string(sizeof(double) * 8);
            } else {
                KRATOS_ERROR << "Unsupported data type in the dynamic dimensional array.";
            }
        }, p_variant_array);
    }

    AddAttribute("type", type);
    AddAttribute("Name", rDataName);

    if (number_of_components > 0) {
        AddAttribute("NumberOfComponents", std::to_string(number_of_components));
    }
}

    ///@}
    ///@name Public operations
    ///@{

const std::string XmlDynamicDimensionalArrayElement::GetTagName() const
{
    return mTagName;
};

void XmlDynamicDimensionalArrayElement::AddAttribute(
    const std::string& rName,
    const std::string& rValue)
{
    for (const auto& r_attribute : mAttributes) {
        KRATOS_ERROR_IF(r_attribute.first == rName)
            << "There exists an attribute named \"" << rName
            << "\" in the xml element with value = \""
            << r_attribute.second << "\" [ given new value = \""
            << rValue << "\" ].\n";
    }
    mAttributes.push_back(std::make_pair(rName, rValue));
}

const std::vector<std::pair<std::string, std::string>>& XmlDynamicDimensionalArrayElement::GetAttributes() const
{
    return mAttributes;
}

void XmlDynamicDimensionalArrayElement::ClearAttributes()
{
    mAttributes.clear();
}

void XmlDynamicDimensionalArrayElement::AddElement(const XmlDynamicDimensionalArrayElement::Pointer pXmlElement)
{
    if (mDynamicDimensionalArrays.size() == 0) {
        for (const auto& p_element : mXmlElements) {
            KRATOS_ERROR_IF(&*(p_element) == &*(pXmlElement))
                << "The xml element is already aded.";
        }
        mXmlElements.push_back(pXmlElement);
    } else {
        KRATOS_ERROR << "Cannot add element to an Xml element which has "
                        "data [ current xml tag = \""
                    << GetTagName() << "\", new element tag name = \""
                    << pXmlElement->GetTagName() << "\" ].\n";
    }
}

std::vector<XmlDynamicDimensionalArrayElement::Pointer> XmlDynamicDimensionalArrayElement::GetElements(const std::string& rTagName) const
{
    std::vector<XmlDynamicDimensionalArrayElement::Pointer> results;
    for (const auto& p_element : mXmlElements) {
        if (p_element->GetTagName() == rTagName) {
            results.push_back(p_element);
        }
    }
    return results;
}

const std::vector<XmlDynamicDimensionalArrayElement::Pointer>& XmlDynamicDimensionalArrayElement::GetElements() const
{
    return mXmlElements;
}

void XmlDynamicDimensionalArrayElement::ClearElements()
{
    mXmlElements.clear();
}

const std::vector<XmlDynamicDimensionalArrayElement::ArrayPointerType> XmlDynamicDimensionalArrayElement::GetDynamicDimensionalArrays() const
{
    return mDynamicDimensionalArrays;
}

} // namespace Kratos
