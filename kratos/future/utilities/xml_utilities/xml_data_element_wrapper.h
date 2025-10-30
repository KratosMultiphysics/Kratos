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
#include <variant>

// Project includes
#include "includes/define.h"
#include "containers/nd_data.h"
#include "xml_element.h"

namespace Kratos::Future {

///@name Kratos Classes
///@{

/* @class XmlDataElementWrapper
 * @ingroup KratosCore
 * @brief Output stream writer for XML format.
 * @author Suneth Warnakulasuriya
 */
class KRATOS_API(KRATOS_CORE) XmlDataElementWrapper : public XmlElement
{
public:
    ///@name Life cycle
    ///@{

    using NDDataPointerType = std::variant<
                                        NDData<unsigned char>::Pointer,
                                        NDData<bool>::Pointer,
                                        NDData<int>::Pointer,
                                        NDData<double>::Pointer
                                    >;

    KRATOS_CLASS_POINTER_DEFINITION(XmlDataElementWrapper);

    ///@}
    ///@name Life cycle
    ///@{

    XmlDataElementWrapper(const std::string& rTagName);

    ~XmlDataElementWrapper() override = default;

    ///@}
    ///@name Public operations
    ///@{

    virtual XmlElement::Pointer Get(
        const std::string& rDataArrayName,
        NDDataPointerType pNDData);

    ///@}
};

} // namespace Kratos::Future