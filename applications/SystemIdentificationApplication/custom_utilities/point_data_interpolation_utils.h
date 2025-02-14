//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <tuple>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

template<class TEntity>
class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) PointDataInterpolationUtils {
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PointDataInterpolationUtils);

    ///@}
    ///@name Life cycle
    ///@{

    PointDataInterpolationUtils(ModelPart& rModelPart);

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    void CalculateInterpolatedNodalValues(
        std::vector<TDataType>& rOutput,
        const Variable<TDataType>& rVariable) const;

    template<class TDataType>
    void CalculateInterpolatedEntityValues(
        std::vector<TDataType>& rOutput,
        const Variable<TDataType>& rVariable) const;

    void UpdatePoints(const std::vector<Point>& rCoordinates);

    ///@}


private:
    ///@name Private member variables
    ///@{

    ModelPart * const mpModelPart;

    std::vector<Point> mCoordinates;

    std::vector<std::tuple<typename TEntity::Pointer, Vector>> mCoordinateEntityData;

    ///@}
};

///@}

} // namespace Kratos