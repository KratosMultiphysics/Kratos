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
#include "expression/literal_flat_expression.h"

// Include base h
#include "xml_expression_element.h"

namespace Kratos {

XmlExpressionElement::XmlExpressionElement(const std::string& rTagName)
    : mTagName(rTagName)
{
}

XmlExpressionElement::XmlExpressionElement(
    const std::string& rDataName,
    const std::vector<Expression::ConstPointer>& rExpressions)
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
        AddAttribute("type", "UInt" + std::to_string(sizeof(char) * 8));
    } else if (std::all_of(mExpressions.begin(), mExpressions.end(), [](const auto& pExpression) {
            return dynamic_cast<const LiteralFlatExpression<int>*>(&*pExpression);
        })) {
        AddAttribute("type", "Int" + std::to_string(sizeof(int) * 8));
    } else {
        AddAttribute("type", "Float" + std::to_string(sizeof(double) * 8));
    }

    AddAttribute("Name", rDataName);

    if (number_of_components > 0) {
        AddAttribute("NumberOfComponents", std::to_string(number_of_components));
    }
}

    ///@}
    ///@name Public operations
    ///@{

const std::string XmlExpressionElement::GetTagName() const
{
    return mTagName;
};

void XmlExpressionElement::AddAttribute(
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

const std::vector<std::pair<std::string, std::string>>& XmlExpressionElement::GetAttributes() const
{
    return mAttributes;
}

void XmlExpressionElement::ClearAttributes()
{
    mAttributes.clear();
}

void XmlExpressionElement::AddElement(const XmlExpressionElement::Pointer pXmlElement)
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

std::vector<XmlExpressionElement::Pointer> XmlExpressionElement::GetElements(const std::string& rTagName) const
{
    std::vector<XmlExpressionElement::Pointer> results;
    for (const auto& p_element : mXmlElements) {
        if (p_element->GetTagName() == rTagName) {
            results.push_back(p_element);
        }
    }
    return results;
}

const std::vector<XmlExpressionElement::Pointer>& XmlExpressionElement::GetElements() const
{
    return mXmlElements;
}

void XmlExpressionElement::ClearElements()
{
    mXmlElements.clear();
}

const std::vector<Expression::ConstPointer> XmlExpressionElement::GetExpressions() const
{
    return mExpressions;
}

} // namespace Kratos
