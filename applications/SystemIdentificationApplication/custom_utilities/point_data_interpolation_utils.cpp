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

// System includes

// External includes

// Project includes

// Application includes
#include "utilities/brute_force_point_locator.h"

// Include base h
#include "point_data_interpolation_utils.h"

namespace Kratos {

template<class TEntity>
PointDataInterpolationUtils<TEntity>::PointDataInterpolationUtils(ModelPart& rModelPart)
    : mpModelPart(&rModelPart)
{
}

template<class TEntity>
template<class TDataType>
void PointDataInterpolationUtils<TEntity>::CalculateInterpolatedNodalValues(
    std::vector<TDataType>& rOutput,
    const Variable<TDataType>& rVariable) const
{
    KRATOS_TRY

    if (rOutput.size() != mCoordinates.size()) {
        rOutput.resize(mCoordinates.size());
    }

    IndexPartition<IndexType>(mCoordinates.size()).for_each([&rOutput, &rVariable, this](const auto Index) {
        const auto& r_geometry = std::get<0>(this->mCoordinateEntityData[Index])->GetGeometry();
        const auto& Ns = std::get<1>(this->mCoordinateEntityData[Index]);
        auto value = rVariable.Zero();
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            value += r_geometry[i].GetValue(rVariable) * Ns[i];
        }

        rOutput[Index] = value;
    });

    KRATOS_CATCH("");
}

template<class TEntity>
template<class TDataType>
void PointDataInterpolationUtils<TEntity>::CalculateInterpolatedEntityValues(
    std::vector<TDataType>& rOutput,
    const Variable<TDataType>& rVariable) const
{
    KRATOS_TRY

    if (rOutput.size() != mCoordinates.size()) {
        rOutput.resize(mCoordinates.size());
    }

    IndexPartition<IndexType>(mCoordinates.size()).for_each([&rOutput, &rVariable, this](const auto Index) {
        rOutput[Index] = std::get<0>(this->mCoordinateEntityData[Index])->GetValue(rVariable);
    });

    KRATOS_CATCH("");
}

template<class TEntity>
void PointDataInterpolationUtils<TEntity>::UpdatePoints(const std::vector<Point>& rCoordinates)
{
    KRATOS_TRY

    mCoordinates = rCoordinates;

    if (mCoordinateEntityData.size() != mCoordinates.size()) {
        mCoordinateEntityData.resize(mCoordinates.size());
    }

    BruteForcePointLocator point_locator(*mpModelPart);

    IndexPartition<IndexType>(rCoordinates.size()).for_each([&point_locator, this](const auto Index) {
        if constexpr(std::is_same_v<TEntity, Element>) {
            const auto entity_id = point_locator.FindElement(this->mCoordinates[Index], std::get<1>(this->mCoordinateEntityData[Index]));
            std::get<0>(this->mCoordinateEntityData[Index]) = this->mpModelPart->pGetElement(entity_id);
        } else if constexpr(std::is_same_v<TEntity, Condition>) {
            const auto entity_id = point_locator.FindCondition(this->mCoordinates[Index], std::get<1>(this->mCoordinateEntityData[Index]));
            std::get<0>(this->mCoordinateEntityData[Index]) = this->mpModelPart->pGetCondition(entity_id);
        }
    });


    KRATOS_CATCH("");
}

// template instantiations
template class PointDataInterpolationUtils<ModelPart::ElementType>;
template void PointDataInterpolationUtils<ModelPart::ElementType>::CalculateInterpolatedNodalValues(std::vector<double>&, const Variable<double>&) const;
template void PointDataInterpolationUtils<ModelPart::ElementType>::CalculateInterpolatedNodalValues(std::vector<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&) const;
template void PointDataInterpolationUtils<ModelPart::ElementType>::CalculateInterpolatedNodalValues(std::vector<array_1d<double, 4>>&, const Variable<array_1d<double, 4>>&) const;
template void PointDataInterpolationUtils<ModelPart::ElementType>::CalculateInterpolatedNodalValues(std::vector<array_1d<double, 6>>&, const Variable<array_1d<double, 6>>&) const;
template void PointDataInterpolationUtils<ModelPart::ElementType>::CalculateInterpolatedNodalValues(std::vector<array_1d<double, 9>>&, const Variable<array_1d<double, 9>>&) const;
template void PointDataInterpolationUtils<ModelPart::ElementType>::CalculateInterpolatedEntityValues(std::vector<double>&, const Variable<double>&) const;
template void PointDataInterpolationUtils<ModelPart::ElementType>::CalculateInterpolatedEntityValues(std::vector<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&) const;
template void PointDataInterpolationUtils<ModelPart::ElementType>::CalculateInterpolatedEntityValues(std::vector<array_1d<double, 4>>&, const Variable<array_1d<double, 4>>&) const;
template void PointDataInterpolationUtils<ModelPart::ElementType>::CalculateInterpolatedEntityValues(std::vector<array_1d<double, 6>>&, const Variable<array_1d<double, 6>>&) const;
template void PointDataInterpolationUtils<ModelPart::ElementType>::CalculateInterpolatedEntityValues(std::vector<array_1d<double, 9>>&, const Variable<array_1d<double, 9>>&) const;

template class PointDataInterpolationUtils<ModelPart::ConditionType>;
template void PointDataInterpolationUtils<ModelPart::ConditionType>::CalculateInterpolatedNodalValues(std::vector<double>&, const Variable<double>&) const;
template void PointDataInterpolationUtils<ModelPart::ConditionType>::CalculateInterpolatedNodalValues(std::vector<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&) const;
template void PointDataInterpolationUtils<ModelPart::ConditionType>::CalculateInterpolatedNodalValues(std::vector<array_1d<double, 4>>&, const Variable<array_1d<double, 4>>&) const;
template void PointDataInterpolationUtils<ModelPart::ConditionType>::CalculateInterpolatedNodalValues(std::vector<array_1d<double, 6>>&, const Variable<array_1d<double, 6>>&) const;
template void PointDataInterpolationUtils<ModelPart::ConditionType>::CalculateInterpolatedNodalValues(std::vector<array_1d<double, 9>>&, const Variable<array_1d<double, 9>>&) const;
template void PointDataInterpolationUtils<ModelPart::ConditionType>::CalculateInterpolatedEntityValues(std::vector<double>&, const Variable<double>&) const;
template void PointDataInterpolationUtils<ModelPart::ConditionType>::CalculateInterpolatedEntityValues(std::vector<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&) const;
template void PointDataInterpolationUtils<ModelPart::ConditionType>::CalculateInterpolatedEntityValues(std::vector<array_1d<double, 4>>&, const Variable<array_1d<double, 4>>&) const;
template void PointDataInterpolationUtils<ModelPart::ConditionType>::CalculateInterpolatedEntityValues(std::vector<array_1d<double, 6>>&, const Variable<array_1d<double, 6>>&) const;
template void PointDataInterpolationUtils<ModelPart::ConditionType>::CalculateInterpolatedEntityValues(std::vector<array_1d<double, 9>>&, const Variable<array_1d<double, 9>>&) const;

}