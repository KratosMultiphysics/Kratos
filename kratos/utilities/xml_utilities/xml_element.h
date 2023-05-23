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
#include "containers/container_expression/expressions/expression.h"
#include "xml_ostream_writer.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) XmlElement {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(XmlElement);

    ///@}
    ///@name Life cycle
    ///@{

    XmlElement(const std::string& rTagName);

    XmlElement(
        const std::string& rDataName,
        const std::vector<Expression::Pointer>& rExpressions);

    ///@}
    ///@name Public operations
    ///@{

    const std::string GetTagName() const;

    void AddAttribute(
        const std::string& rName,
        const std::string& rValue);

    const std::vector<std::pair<const std::string, const std::string>>& GetAttributes() const;

    void ClearAttributes();

    void AddElement(const XmlElement::Pointer pXmlElement);

    std::vector<XmlElement::Pointer> GetElements(const std::string& rTagName) const;

    const std::vector<XmlElement::Pointer>& GetElements() const;

    void ClearElements();

    void Write(
        XmlOStreamWriter& rWriter,
        const IndexType Level = 0) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const std::string mTagName;

    std::vector<std::pair<const std::string, const std::string>> mAttributes;

    std::vector<XmlElement::Pointer> mXmlElements;

    const std::vector<Expression::Pointer> mExpressions;

    ///@}
};

///@}

} // namespace Kratos
