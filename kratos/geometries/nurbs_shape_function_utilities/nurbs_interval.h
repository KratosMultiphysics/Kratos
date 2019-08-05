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

namespace Kratos
{

class Interval
{
public:
    ///@name Life Cycle
    ///@{
    Interval()
        : m_t0(std::numeric_limits<double>::quiet_NaN())
        , m_t1(std::numeric_limits<double>::quiet_NaN())
    {
    }

    Interval(const double T0, const double T1)
        : m_t0(T0)
        , m_t1(T1)
    {
    }

    Interval(const std::pair<double, double> Bounds)
        : m_t0(Bounds.first)
        , m_t1(Bounds.second)
    {
    }

    ///@}
    ///@name Get/ Set Operators
    ///@{

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

    double MinParameter() const
    {
        return std::min(m_t0, m_t1);
    }

    double MaxParameter() const
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

    ///@}
    ///@name Normalized Operators
    ///@{

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

    ///@}
    ///@name Operators
    ///@{

    bool IsInside(double& ParameterT) const
    {
        const double min_parameter = MinParameter();
        if (ParameterT < min_parameter) {
            ParameterT = min_parameter;
            return false;
        }

        const double max_parameter = MaxParameter();
        if (ParameterT > max_parameter) {
            ParameterT = max_parameter;
            return false;
        }

        return true;
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    double m_t0;
    double m_t1;

    ///@}
}; // Interval

} // namespace Kratos

#endif // KRATOS_NURBS_INTERVAL_H_INCLUDED defined 