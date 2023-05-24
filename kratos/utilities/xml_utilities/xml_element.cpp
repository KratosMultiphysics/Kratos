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

// Project includes
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"

// Include base h
#include "xml_element.h"

namespace Kratos {

XmlElement::XmlElement(const std::string& rTagName)
    : mTagName(rTagName)
{
}

XmlElement::XmlElement(
    const std::string& rDataName,
    const std::vector<Expression::Pointer>& rExpressions)
    : mTagName("DataArray"),
    mExpressions(rExpressions)
{
    KRATOS_ERROR_IF(rExpressions.size() == 0)
        << "Empty expression lists are not allowed.";

    IndexType number_of_components = 0;

    for (const auto& p_expression : rExpressions) {
        if (number_of_components == 0) {
            number_of_components = p_expression->GetItemComponentCount();
        }
        KRATOS_ERROR_IF(number_of_components != p_expression->GetItemComponentCount())
            << "Found expressions with mismatching shapes.";
    }

    if (std::all_of(mExpressions.begin(), mExpressions.end(), [](const auto& pExpression) {
            return dynamic_cast<const LiteralFlatExpression<char>*>(&*pExpression);
        })) {
        AddAttribute("type", "UInt8");
    } else if (std::all_of(mExpressions.begin(), mExpressions.end(), [](const auto& pExpression) {
            return dynamic_cast<const LiteralFlatExpression<int>*>(&*pExpression);
        })) {
        AddAttribute("type", "Int32");
    } else {
        AddAttribute("type", "Float64");
    }

    AddAttribute("Name", rDataName);

    if (number_of_components > 0) {
        AddAttribute("NumberOfComponents", std::to_string(number_of_components));
    }
}

    ///@}
    ///@name Public operations
    ///@{

const std::string XmlElement::GetTagName() const
{
    return mTagName;
};

void XmlElement::AddAttribute(
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

const std::vector<std::pair<std::string, std::string>>& XmlElement::GetAttributes() const
{
    return mAttributes;
}

void XmlElement::ClearAttributes()
{
    mAttributes.clear();
}

void XmlElement::AddElement(const XmlElement::Pointer pXmlElement)
{
    if (mExpressions.size() == 0) {
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

std::vector<XmlElement::Pointer> XmlElement::GetElements(const std::string& rTagName) const
{
    std::vector<XmlElement::Pointer> results;
    for (const auto& p_element : mXmlElements) {
        if (p_element->GetTagName() == rTagName) {
            results.push_back(p_element);
        }
    }
    return results;
}

const std::vector<XmlElement::Pointer>& XmlElement::GetElements() const
{
    return mXmlElements;
}

void XmlElement::ClearElements()
{
    mXmlElements.clear();
}

void XmlElement::Write(
    XmlOStreamWriter& rWriter,
    const IndexType Level) const
{
    if (mXmlElements.size() > 0) {
        rWriter.WriteElement(GetTagName(), GetAttributes(), Level, false);
        for (const auto& p_element : mXmlElements) {
            p_element->Write(rWriter, Level + 1);
        }
        rWriter.CloseElement(GetTagName(), Level);
    } else if (mExpressions.size() > 0) {
        rWriter.WriteDataElement(GetTagName(), GetAttributes(), mExpressions, Level);
    } else {
        rWriter.WriteElement(GetTagName(), GetAttributes(), Level, true);
    }
}

} // namespace Kratos
