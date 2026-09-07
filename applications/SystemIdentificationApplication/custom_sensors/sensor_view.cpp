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

// System includes
#include <type_traits>
#include <sstream>

// External includes

// Project includes
#include "includes/model_part.h"

// Application includes

// Include base h
#include "sensor_view.h"

namespace Kratos {

SensorView::SensorView(
    Sensor::Pointer pSensor,
    const std::string& rTensorAdaptorName)
    : mpSensor(pSensor),
      mTensorAdaptorName(rTensorAdaptorName),
      mpTensorAdaptor(pSensor->GetTensorAdaptor(rTensorAdaptorName))
{
}

Sensor::Pointer SensorView::GetSensor() const
{
    return mpSensor;
}

TensorAdaptor<double>::Pointer SensorView::GetTensorAdaptor() const
{
    return mpTensorAdaptor;
}

std::string SensorView::GetTensorAdaptorName() const
{
    return mTensorAdaptorName;
}

void  SensorView::AddAuxiliaryTensorAdaptor(
    const std::string& rSuffix,
    TensorAdaptor<double>::Pointer pTensorAdaptor)
{
    std::stringstream name;
    name << this->mTensorAdaptorName << "_" << rSuffix;
    mpSensor->AddTensorAdaptor(name.str(), pTensorAdaptor);
}

TensorAdaptor<double>::Pointer SensorView::GetAuxiliaryTensorAdaptor(const std::string& rSuffix) const
{
    std::stringstream name;
    name << this->mTensorAdaptorName << "_" << rSuffix;
    return mpSensor->GetTensorAdaptor(name.str());
}

std::vector<std::string> SensorView::GetAuxiliarySuffixes() const
{
    KRATOS_TRY

    std::vector<std::string> suffixes;
    for (const auto& r_pair : mpSensor->GetTensorAdaptorsMap()) {
        const auto& r_name = r_pair.first;
        if (r_name.rfind(mTensorAdaptorName + "_", 0) == 0) {
            suffixes.push_back(r_name.substr(mTensorAdaptorName.size() + 1));
        }
    }

    return suffixes;

    KRATOS_CATCH("");
}

std::string SensorView::Info() const
{
    return mpSensor->Info();
}

void SensorView::PrintInfo(std::ostream& rOStream) const
{
    mpSensor->PrintInfo(rOStream);
}

void SensorView::PrintData(std::ostream& rOStream) const
{
    mpSensor->PrintData(rOStream);
}

} /* namespace Kratos.*/