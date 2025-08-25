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
#include "xml_nd_data_element.h"

namespace Kratos {

XmlNDDataElement::XmlNDDataElement(const std::string& rTagName)
    : mTagName(rTagName)
{
}

XmlNDDataElement::XmlNDDataElement(
    const std::string& rDataName,
    const std::vector<ArrayPointerType>& rListOfNDDataPointers)
    : mTagName("DataArray"),
      mListOfNDData(rListOfNDDataPointers)
{
    KRATOS_ERROR_IF(mListOfNDData.size() == 0)
        << "Empty list of N-Dimensional data is not allowed.";

    IndexType number_of_components = 0;
    std::string type = "";
    for (auto& p_variant_data : mListOfNDData) {
        std::visit([&number_of_components, &type](auto pData) {
            if (number_of_components == 0) {
                number_of_components = pData->Size() / pData->Shape()[0];
            }

            KRATOS_ERROR_IF(pData->Shape()[0] == 0 || number_of_components != pData->Size() / pData->Shape()[0])
                << "Found dynamic dimensional array with zero elements or with mismatching shapes. "
                << " [ dynamic dimensional array shape = " << pData->Shape() << ", number of components = "
                << number_of_components << " ].\n";

            if constexpr(std::is_same_v<typename BareType<decltype(*pData)>::DataType, unsigned char>) {
                if (type == "") type = "UInt" + std::to_string(sizeof(char) * 8);
                else if (type != "UInt") type = "Float" + std::to_string(sizeof(double) * 8);
            } else if constexpr(std::is_same_v<typename BareType<decltype(*pData)>::DataType, bool>) {
                if (type == "") type = "UInt" + std::to_string(sizeof(char) * 8);
                else if (type != "UInt") type = "Float" + std::to_string(sizeof(double) * 8);
            } else if constexpr(std::is_same_v<typename BareType<decltype(*pData)>::DataType, int>) {
                if (type == "") type = "Int" + std::to_string(sizeof(int) * 8);
                else if (type != "Int") type = "Float" + std::to_string(sizeof(double) * 8);
            } else if constexpr(std::is_same_v<typename BareType<decltype(*pData)>::DataType, double>) {
                if (type == "") type = "Float" + std::to_string(sizeof(double) * 8);
                else if (type != "Float") type = "Float" + std::to_string(sizeof(double) * 8);
            } else {
                KRATOS_ERROR << "Unsupported data type in the dynamic dimensional array.";
            }
        }, p_variant_data);
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

const std::string XmlNDDataElement::GetTagName() const
{
    return mTagName;
};

void XmlNDDataElement::AddAttribute(
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

const std::vector<std::pair<std::string, std::string>>& XmlNDDataElement::GetAttributes() const
{
    return mAttributes;
}

void XmlNDDataElement::ClearAttributes()
{
    mAttributes.clear();
}

void XmlNDDataElement::AddElement(const XmlNDDataElement::Pointer pXmlElement)
{
    if (mListOfNDData.size() == 0) {
        for (const auto& p_element : mXmlElements) {
            KRATOS_ERROR_IF(&*(p_element) == &*(pXmlElement))
                << "The xml element is already added.";
        }
        mXmlElements.push_back(pXmlElement);
    } else {
        KRATOS_ERROR << "Cannot add element to an Xml element which has "
                        "data [ current xml tag = \""
                    << GetTagName() << "\", new element tag name = \""
                    << pXmlElement->GetTagName() << "\" ].\n";
    }
}

std::vector<XmlNDDataElement::Pointer> XmlNDDataElement::GetElements(const std::string& rTagName) const
{
    std::vector<XmlNDDataElement::Pointer> results;
    for (const auto& p_element : mXmlElements) {
        if (p_element->GetTagName() == rTagName) {
            results.push_back(p_element);
        }
    }
    return results;
}

const std::vector<XmlNDDataElement::Pointer>& XmlNDDataElement::GetElements() const
{
    return mXmlElements;
}

void XmlNDDataElement::ClearElements()
{
    mXmlElements.clear();
}

const std::vector<XmlNDDataElement::ArrayPointerType> XmlNDDataElement::GetListOfNDData() const
{
    return mListOfNDData;
}

} // namespace Kratos
