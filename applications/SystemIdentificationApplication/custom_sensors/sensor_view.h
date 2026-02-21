//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "tensor_adaptors/tensor_adaptor.h"
#include "includes/define.h"

// Application includes
#include "sensor.h"

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SensorView
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensorView);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SensorView(
        Sensor::Pointer pSensor,
        const std::string& rTensorAdaptorName);

    ///@}
    ///@name Public operations
    ///@{

    Sensor::Pointer GetSensor() const;

    TensorAdaptor<double>::Pointer GetTensorAdaptor() const;

    std::string GetTensorAdaptorName() const;

    void  AddAuxiliaryTensorAdaptor(
        const std::string& rSuffix,
        TensorAdaptor<double>::Pointer pTensorAdaptor);

    TensorAdaptor<double>::Pointer GetAuxiliaryTensorAdaptor(const std::string& rSuffix) const;

    std::vector<std::string> GetAuxiliarySuffixes() const;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const;

    void PrintInfo(std::ostream& rOStream) const;

    void PrintData(std::ostream& rOStream) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    Sensor::Pointer mpSensor;

    const std::string mTensorAdaptorName;

    TensorAdaptor<double>::Pointer mpTensorAdaptor;

    ///@}
};

/// output stream functions
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const SensorView& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@} // Kratos Classes

} /* namespace Kratos.*/