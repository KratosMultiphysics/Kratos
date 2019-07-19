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

    Interval(const double T0, const double T1) : m_t0(T0), m_t1(T1)
    {
    }

    Interval(const std::pair<double, double> Bounds) : m_t0(Bounds.first),
        m_t1(Bounds.second)
    {
    }

public:     // methods
    double GetT0() const
    {
        return m_t0;
    }

    void SetT0(const double Value)
    {
        m_t0 = Value;
    }

    double GetT1() const
    {
        return m_t1;
    }

    void SetT1(const double Value)
    {
        m_t1 = Value;
    }

    double GetMin() const
    {
        return std::min(m_t0, m_t1);
    }

    double GetMax() const
    {
        return std::max(m_t0, m_t1);
    }

    double GetDelta() const
    {
        return m_t1 - m_t0;
    }

    double GetLength() const
    {
        return std::abs(GetDelta());
    }

    double GetNormalizedAt(const double Parameter) const
    {
        return (Parameter - m_t0) / GetLength();
    }

    double GetParameterAtNormalized(const double Parameter) const
    {
        return m_t0 + GetDelta() * Parameter;
    }

    static double GetParameterAtNormalized(const double A, const double B,
        const double Parameter)
    {
        return A + (B - A) * Parameter;
    }

    Interval GetNormalizedInterval(const double T0, const double T1) const
    {
        double t0Normalized = GetNormalizedAt(T0);
        double t1Normalized = GetNormalizedAt(T1);

        return Interval(t0Normalized, t1Normalized);
    }

    Interval GetNormalizedInterval(const Interval Bounds) const
    {
        return GetNormalizedInterval(Bounds.m_t0, Bounds.m_t1);
    }

    bool Contains(const double Parameter) const
    {
        return (GetMin() <= Parameter) && (Parameter <= GetMax());
    }

    double GetClamp(const double Value) const
    {
        return GetClamp(Value, GetMin(), GetMax());
    }

public:     // static methods
    static double GetClamp(const double Value, const double Min,
        const double Max)
    {
        if (Value < Min) {
            return Min;
        }

        if (Value > Max) {
            return Max;
        }

        return Value;
    }

    static double GetClamp01(const double Value)
    {
        return GetClamp(Value, 0, 1);
    }
}; // Interval

} // namespace Kratos

#endif // KRATOS_NURBS_INTERVAL_H_INCLUDED defined