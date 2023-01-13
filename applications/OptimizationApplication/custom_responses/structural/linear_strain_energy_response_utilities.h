//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Reza Najian Asl, https://github.com/RezaNajian
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) LinearStrainEnergyResponseUtilities
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using GeometryType = ModelPart::ElementType::GeometryType;

    ///@}
    ///@name Static operations
    ///@{

    static double CalculateStrainEnergy(ModelPart& rModelPart);

    static void CalculateStrainEnergyShapeSensitivity(
        ModelPart& rModelPart,
        const double Delta,
        const Variable<array_1d<double, 3>>& rOutputSensitivityVariable);

    static void CalculateStrainEnergyElementPropertiesSensitivity(
        ModelPart& rModelPart,
        const double Delta,
        const Variable<double>& rPrimalVariable,
        const Variable<double>& rOutputSensitivityVariable);

    ///@}
private:
    ///@name Private operations
    ///@{

    ///@}
};

///@}
}