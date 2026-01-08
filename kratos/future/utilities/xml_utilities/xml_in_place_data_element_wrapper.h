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

// Project includes
#include "includes/define.h"
#include "containers/nd_data.h"
#include "xml_data_element_wrapper.h"

namespace Kratos::Future {

///@name Kratos Classes
///@{

/* @class XmlInPlaceDataElementWrapper
 * @ingroup KratosCore
 * @brief Output stream writer for XML format.
 * @author Suneth Warnakulasuriya
 */
class KRATOS_API(KRATOS_CORE) XmlInPlaceDataElementWrapper : public XmlDataElementWrapper
{
public:
    ///@name Type definitions
    ///@{

    using BaseType = XmlDataElementWrapper;

    KRATOS_CLASS_POINTER_DEFINITION(XmlInPlaceDataElementWrapper);

    ///@}
    ///@name Enums
    ///@{

    enum XmlOutputType
    {
        ASCII = 0,
        BINARY = 1
    };

    ///@}
    ///@name Life cycle
    ///@{

    XmlInPlaceDataElementWrapper(
        const XmlOutputType OutputType,
        const IndexType Precision);

    ~XmlInPlaceDataElementWrapper() override = default;

    ///@}
    ///@name Public operations
    ///@{

    XmlElement::Pointer Get(
        const std::string& rDataArrayName,
        NDDataPointerType pNDData) override;

    void Write(
        std::ostream& rOStream,
        const IndexType Level = 0) const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const XmlOutputType mOutputType;

    const IndexType mPrecision;

    ///@}
};

} // namespace Kratos::Future