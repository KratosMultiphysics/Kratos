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

// Project includes
#include "includes/define.h"
#include "containers/nd_data.h"
#include "xml_ascii_nd_data_element.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TDataType>
class KRATOS_API(KRATOS_CORE) XmlBase64BinaryNDDataElement : public XmlAsciiNDDataElement<TDataType>
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using BaseType = XmlAsciiNDDataElement<TDataType>;

    KRATOS_CLASS_POINTER_DEFINITION(XmlBase64BinaryNDDataElement);

    ///@}
    ///@name Life cycle
    ///@{

    XmlBase64BinaryNDDataElement(
        const std::string& rDataArrayName,
        typename NDData<TDataType>::Pointer pNDData);

    ///@}
    ///@name Public operations
    ///@{

    void Write(
        std::ostream& rOStream,
        const IndexType Level = 0) const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    typename NDData<TDataType>::Pointer mpNDData;

    ///@}
};

///@}

} // namespace Kratos
