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

class KRATOS_API(OPTIMIZATION_APPLICATION) MassResponseUtilities
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using GeometryType = ModelPart::ElementType::GeometryType;

    ///@}
    ///@name Static operations
    ///@{

    static double CalculateMass(const ModelPart& rModelPart);

    static void CalculateMassShapeSensitivity(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rOutputSensitivityVariable);

    static void CalculateMassDensitySensitivity(
        ModelPart& rModelPart,
        const Variable<double>& rOutputSensitivityVariable);

    static void CalculateMassThicknessSensitivity(
        ModelPart& rModelPart,
        const Variable<double>& rOutputSensitivityVariable);

    static void CalculateMassCrossAreaSensitivity(
        ModelPart& rModelPart,
        const Variable<double>& rOutputSensitivityVariable);

    ///@}
private:
    ///@name Private operations
    ///@{

    static bool HasVariableInProperties(
        const ModelPart& rModelPart,
        const Variable<double>& rVariable);

    static void CalculateMassGeometricalPropertySensitivity(
        ModelPart& rModelPart,
        const Variable<double>& rGeometricalPropertySensitivityVariable,
        const Variable<double>& rGeometricalCoflictingPropertySensitivityVariable,
        const Variable<double>& rOutputSensitivityVariable);

    ///@}
};

///@}
}