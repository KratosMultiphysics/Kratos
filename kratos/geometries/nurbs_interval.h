//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

#if !defined(KRATOS_NURBS_INTERVAL_H_INCLUDED )
#define  KRATOS_NURBS_INTERVAL_H_INCLUDED

#include <algorithm>
#include <limits>
#include <sstream>

namespace Kratos {

class Interval
{
private:    // variables
    double m_t0;
    double m_t1;

public:     // constructor
    Interval() : m_t0(std::numeric_limits<double>::quiet_NaN()),
        m_t1(std::numeric_limits<double>::quiet_NaN())
    {
    }

    Interval(double t0, double t1) : m_t0(t0), m_t1(t1)
    {
    }

    Interval(std::pair<double, double> interval) : m_t0(interval.first),
        m_t1(interval.second)
    {
    }

public:     // methods
    double t0() const
    {
        return m_t0;
    }

    void set_t0(double value)
    {
        m_t0 = value;
    }

    double t1() const
    {
        return m_t1;
    }

    void set_t1(double value)
    {
        m_t1 = value;
    }

    double min() const
    {
        return std::min(m_t0, m_t1);
    }

    double max() const
    {
        return std::max(m_t0, m_t1);
    }

    double delta() const
    {
        return m_t1 - m_t0;
    }

    double length() const
    {
        return std::abs(delta());
    }

    double normalized_at(const double t) const
    {
        return (t - m_t0) / length();
    }

    double parameter_at_normalized(const double t) const
    {
        return m_t0 + delta() * t;
    }

    static double parameter_at_normalized(const double a, const double b,
        const double t)
    {
        return a + (b - a) * t;
    }

    Interval normalized_interval(const double t0, const double t1) const
    {
        double t0Normalized = normalized_at(t0);
        double t1Normalized = normalized_at(t1);

        return Interval(t0Normalized, t1Normalized);
    }

    Interval normalized_interval(const Interval interval) const
    {
        return normalized_interval(interval.m_t0, interval.m_t1);
    }

    bool contains(const double t) const
    {
        return (min() <= t) && (t <= max());
    }

    double clamp(const double value) const
    {
        return clamp(value, min(), max());
    }

public:     // static methods
    static double clamp(const double value, const double min, const double max)
    {
        if (value < min) {
            return min;
        }

        if (value > max) {
            return max;
        }

        return value;
    }

    static double clamp_01(const double value)
    {
        return clamp(value, 0, 1);
    }
}; // Interval

} // namespace Kratos

#endif // KRATOS_NURBS_INTERVAL_H_INCLUDED defined