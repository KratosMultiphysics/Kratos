

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "closest_points.h"

namespace Kratos {

PointWithId::PointWithId(const IndexType NewId, const CoordinatesArrayType& rCoords, const double Distance)
    : IndexedObject(NewId), Point(rCoords), mDistance(Distance)
{
    KRATOS_ERROR_IF(mDistance<0.0) << "Distance cannot be negative!" << std::endl;
}

bool PointWithId::operator<(const PointWithId& rOther) const
{
    return mDistance < rOther.mDistance;
}

void PointWithId::save(Serializer &rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject);
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Point);
    rSerializer.save("distance", mDistance);
}

void PointWithId::load(Serializer &rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject);
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Point);
    rSerializer.load("distance", mDistance);
}


ClosestPointsContainer::ClosestPointsContainer(const std::size_t MaxSize) : mMaxSize(MaxSize) { }

ClosestPointsContainer::ClosestPointsContainer(const std::size_t MaxSize, const double MaxDistance)
    : mMaxSize(MaxSize), mMaxDistance(MaxDistance) { }


void ClosestPointsContainer::Add(const PointWithId& rPoint)
{
    if (rPoint.GetDistance() > mMaxDistance) return;
    if (mClosestPoints.size() >= mMaxSize && rPoint.GetDistance() > mClosestPoints.rbegin()->GetDistance()) return;

    mClosestPoints.insert(rPoint);

    LimitToMaxSize();
}


void ClosestPointsContainer::Merge(const ClosestPointsContainer& rOther)
{
    KRATOS_ERROR_IF(std::abs(mMaxDistance - rOther.mMaxDistance) > 1E-12) << "Maximum allowed distanced don't match!" << std::endl;

    mClosestPoints.insert(rOther.GetPoints().begin(), rOther.GetPoints().end());

    LimitToMaxSize();
}


void ClosestPointsContainer::LimitToMaxSize()
{
    // if maximum size is exceeded, then remove points with largest distance
    if (mClosestPoints.size() > mMaxSize) {
        auto new_end = mClosestPoints.begin();
        std::advance(new_end, mMaxSize);
        mClosestPoints.erase(new_end, mClosestPoints.end());
    }
}

void ClosestPointsContainer::save(Serializer &rSerializer) const
{
    // rSerializer.save("closest_points", mClosestPoints);
    rSerializer.save("max_size", mMaxSize);
    rSerializer.save("max_distance", mMaxDistance);
}

void ClosestPointsContainer::load(Serializer &rSerializer)
{
    // rSerializer.load("closest_points", mClosestPoints);
    rSerializer.load("max_size", mMaxSize);
    rSerializer.load("max_distance", mMaxDistance);
}

}  // namespace Kratos.
