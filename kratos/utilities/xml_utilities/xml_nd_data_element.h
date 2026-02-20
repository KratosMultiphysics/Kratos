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
#include "xml_element.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TDataType>
class KRATOS_API(KRATOS_CORE) XmlNDDataElement : public XmlElement {
public:
    ///@name Type definitions
    ///@{

    using BaseType = XmlElement;

    KRATOS_CLASS_POINTER_DEFINITION(XmlNDDataElement);

    ///@}
    ///@name Life cycle
    ///@{

    XmlNDDataElement(
        const std::string& rDataArrayName,
        typename NDData<TDataType>::Pointer pNDData);

    ///@}

protected:
    ///@name Private member variables
    ///@{

    typename NDData<TDataType>::Pointer mpNDData;

    ///@}
};

///@}

} // namespace Kratos
