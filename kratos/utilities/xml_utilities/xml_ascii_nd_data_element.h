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
#include "xml_nd_data_element.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TDataType>
class KRATOS_API(KRATOS_CORE) XmlAsciiNDDataElement : public XmlNDDataElement<TDataType> {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using BaseType = XmlNDDataElement<TDataType>;

    KRATOS_CLASS_POINTER_DEFINITION(XmlAsciiNDDataElement);

    ///@}
    ///@name Life cycle
    ///@{

    XmlAsciiNDDataElement(
        const std::string& rDataArrayName,
        typename NDData<TDataType>::Pointer pNDData,
        const IndexType Precision);

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

    const IndexType mPrecision;

    ///@}
};

///@}

} // namespace Kratos
