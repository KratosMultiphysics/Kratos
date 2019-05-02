#pragma once

#include <array>
#include <assert.h>

namespace ANurbs {

template <typename TScalar, int TDimension>
class Point
{
public:
    using PointType = Point<TScalar, TDimension>;
    using ScalarType = TScalar;

private:
    std::array<double, TDimension> m_coordinates;

public:
    Point()
        : m_coordinates({})
    {
    }

    Point(
        const Point& other)
        : m_coordinates(other.m_coordinates)
    {
    }

    Point(
        std::initializer_list<ScalarType> coordinates)
    {
        assert(coordinates.size() == TDimension && "Wrong dimension");

        auto it = coordinates.begin();
        for (int i = 0; i < Dimension(); i++) {
            m_coordinates[i] = *it++;
        }
    }

    static constexpr int
    Dimension()
    {
        return TDimension;
    }

    const ScalarType
    operator[](
        const int index) const
    {
        return m_coordinates[index];
    }

    ScalarType&
    operator[](
        const int index)
    {
        assert(index < Dimension() && "Index out of range");

        return m_coordinates[index];
    }

    const ScalarType
    X() const
    {
        static_assert(Dimension() >= 1, "Index out of range");
        return m_coordinates[0];
    }

    void
    SetX(
        const ScalarType value)
    {
        static_assert(Dimension() >= 1, "Index out of range");
        m_coordinates[0] = value;
    }

    const ScalarType
    Y() const
    {
        static_assert(Dimension() >= 2, "Index out of range");
        return m_coordinates[1];
    }

    void
    SetY(
        const ScalarType value)
    {
        static_assert(Dimension() >= 2, "Index out of range");
        m_coordinates[1] = value;
    }

    const ScalarType
    Z() const
    {
        static_assert(Dimension() >= 3, "Index out of range");
        return m_coordinates[2];
    }

    void
    SetZ(
        const ScalarType value)
    {
        static_assert(Dimension() >= 3, "Index out of range");
        m_coordinates[2] = value;
    }

    friend PointType
    operator+(
        PointType lhs,
        const PointType& rhs)
    {
        for (int i = 0; i < Dimension(); i++) {
            lhs[i] += rhs[i];
        }

        return lhs;
    }

    friend PointType
    operator+=(
        PointType& lhs,
        const PointType& rhs)
    {
        for (int i = 0; i < Dimension(); i++) {
            lhs[i] += rhs[i];
        }

        return lhs;
    }

    friend PointType
    operator-(
        PointType lhs,
        const PointType& rhs)
    {
        for (int i = 0; i < Dimension(); i++) {
            lhs[i] -= rhs[i];
        }

        return lhs;
    }

    friend PointType
    operator*(
        PointType lhs,
        const ScalarType rhs)
    {
        for (int i = 0; i < Dimension(); i++) {
            lhs[i] *= rhs;
        }

        return lhs;
    }
};

using Point1D = Point<double, 1>;
using Point2D = Point<double, 2>;
using Point3D = Point<double, 3>;

} // namespace ANurbs
