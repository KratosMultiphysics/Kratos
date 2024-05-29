

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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
    KRATOS_ERROR_IF(mDistance < 0.0) << "Distance cannot be negative!" << std::endl;
}

PointWithId::PointWithId(const PointWithId& rOther)
    : IndexedObject(rOther), Point(rOther), mDistance(rOther.mDistance) { }

bool PointWithId::operator<(const PointWithId& rOther) const
{
    return (!Point::operator==(rOther)) && (mDistance < rOther.mDistance);
}

bool PointWithId::operator==(const PointWithId& rOther) const
{
    return ( Point::operator==(rOther));
}

bool PointWithId::operator!=(const PointWithId& rOther) const
{
    return (!Point::operator==(rOther));
}

bool PointWithId::operator==(const PointWithId& rOther) const
{
    return Point::operator==(rOther);
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

ClosestPointsContainer::ClosestPointsContainer(const ClosestPointsContainer& rOther)
    : mClosestPoints(rOther.mClosestPoints), mMaxSize(rOther.mMaxSize), mMaxDistance(rOther.mMaxDistance) { }

bool ClosestPointsContainer::operator==(const ClosestPointsContainer& rOther) const
{
    // basic checks
    if (this->mClosestPoints.size() != rOther.mClosestPoints.size()) return false;
    if (this->mMaxSize != rOther.mMaxSize) return false;
    if (this->mMaxDistance != rOther.mMaxDistance) return false;

    // check points
    auto it_point = this->mClosestPoints.begin();
    auto it_ref_point = rOther.mClosestPoints.begin();

    for (std::size_t i=0; i<mClosestPoints.size(); ++i) {
        const PointWithId& r_point = *it_point;
        const PointWithId& r_ref_point = *it_ref_point;
        if (!(r_point == r_ref_point)) return false; // calls Point::operator==
        if (r_point.GetId() != r_ref_point.GetId()) return false;
        if (std::abs(r_point.GetDistance()-r_ref_point.GetDistance()) > 1E-12) return false;

        std::advance(it_point, 1);
        std::advance(it_ref_point, 1);
    }

    return true;
}

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
    rSerializer.save("closest_points", mClosestPoints);
    rSerializer.save("max_size", mMaxSize);
    rSerializer.save("max_distance", mMaxDistance);
}

void ClosestPointsContainer::load(Serializer &rSerializer)
{
    rSerializer.load("closest_points", mClosestPoints);
    rSerializer.load("max_size", mMaxSize);
    rSerializer.load("max_distance", mMaxDistance);
}

}  // namespace Kratos.
