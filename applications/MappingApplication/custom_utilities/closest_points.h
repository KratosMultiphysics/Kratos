
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

#pragma once

// System includes
#include <set>
#include <limits>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/indexed_object.h"
#include "includes/serializer.h"
#include "geometries/point.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(MAPPING_APPLICATION) PointWithId : public IndexedObject, public Point
{
public:
    using IndexedObject::IndexType;

    PointWithId(const IndexType NewId, const CoordinatesArrayType& rCoords, const double Distance);

    // Copy constructor
    PointWithId(const PointWithId& rOther);

    PointWithId& operator=(const PointWithId& rOther) = delete;

    bool operator<(const PointWithId& rOther) const;

    bool operator!=(const PointWithId& rOther) const;

    bool operator==(const PointWithId& rOther) const;

    double GetDistance() const { return mDistance; }

private:
    ///@name Member Variables
    ///@{

    double mDistance;

    ///@}

    ///@name Serialization
    ///@{

    // default CTor for Serializer
    PointWithId() = default;

    friend class Serializer;

    void save(Serializer &rSerializer) const override;

    void load(Serializer &rSerializer) override;

    ///@}
};


class KRATOS_API(MAPPING_APPLICATION) ClosestPointsContainer
{
public:
    using ContainerType = std::set<PointWithId>;

    explicit ClosestPointsContainer(const std::size_t MaxSize);
    ClosestPointsContainer(const std::size_t MaxSize, const double MaxDistance);

    // Copy constructor
    explicit ClosestPointsContainer(const ClosestPointsContainer& rOther);

    ClosestPointsContainer& operator=(const ClosestPointsContainer& rOther) = delete;

    bool operator==(const ClosestPointsContainer& rOther) const;

    void Add(const PointWithId& rPoint);

    void Merge(const ClosestPointsContainer& rOther);

    ContainerType& GetPoints() { return mClosestPoints; }

    const ContainerType& GetPoints() const { return mClosestPoints; }

private:
    ///@name Member Variables
    ///@{

    ContainerType mClosestPoints;
    std::size_t mMaxSize;
    double mMaxDistance = std::numeric_limits<double>::max();

    ///@}

    void LimitToMaxSize();

    ///@name Serialization
    ///@{

    // default CTor for Serializer
    ClosestPointsContainer() = default;

    friend class Serializer;

    void save(Serializer &rSerializer) const;

    void load(Serializer &rSerializer);

    ///@}
};


///@} addtogroup block
}  // namespace Kratos.