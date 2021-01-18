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

// System includes
#include <algorithm>
#include <limits>
#include <sstream>

// External includes

// Project includes

namespace Kratos
{

/// Class for optimized use of intervals
/*
* Provides universal geometrical utiltity functions for the computation of
* curve and surface NURBS/ B-Spline shape functions.
*/
class NurbsInterval
{
public:
    ///@name Life Cycle
    ///@{
    NurbsInterval()
    {
    }

    NurbsInterval(const double T0, const double T1)
        : mT0(T0)
        , mT1(T1)
    {
    }

    NurbsInterval(const std::pair<double, double> Bounds)
        : mT0(Bounds.first)
        , mT1(Bounds.second)
    {
    }

    ///@}
    ///@name Get/ Set Operators
    ///@{

    double GetT0() const
    {
        return mT0;
    }

    void SetT0(const double Value)
    {
        mT0 = Value;
    }

    double GetT1() const
    {
        return mT1;
    }

    void SetT1(const double Value)
    {
        mT1 = Value;
    }

    double MinParameter() const
    {
        return std::min(mT0, mT1);
    }

    double MaxParameter() const
    {
        return std::max(mT0, mT1);
    }

    double GetDelta() const
    {
        return mT1 - mT0;
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
        return (Parameter - mT0) / GetLength();
    }

    double GetParameterAtNormalized(const double Parameter) const
    {
        return mT0 + GetDelta() * Parameter;
    }

    static double GetParameterAtNormalized(const double A, const double B,
        const double Parameter)
    {
        return A + (B - A) * Parameter;
    }

    NurbsInterval GetNormalizedInterval(const double T0, const double T1) const
    {
        double t0Normalized = GetNormalizedAt(T0);
        double t1Normalized = GetNormalizedAt(T1);

        return NurbsInterval(t0Normalized, t1Normalized);
    }

    NurbsInterval GetNormalizedInterval(const NurbsInterval Bounds) const
    {
        return GetNormalizedInterval(Bounds.mT0, Bounds.mT1);
    }

    ///@}
    ///@name Operators
    ///@{

    /*
    * @brief Detects if ParameterT lays within the limits of this interval.
    *        Keeps the ParameterT or sets it to the value of the insulted
    *        bounding interval limit.
    * @return true if ParameterT is inside.
    */
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

    /// @brief Lower bound of the NurbsInterval
    double mT0;

    /// @brief Upper bound of the NurbsInterval
    double mT1;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        rSerializer.save("T0", mT0);
        rSerializer.save("T1", mT1);
    }

    void load(Serializer& rSerializer)
    {
        rSerializer.load("T0", mT0);
        rSerializer.load("T1", mT1);
    }

    ///@}
}; // NurbsInterval

} // namespace Kratos

#endif // KRATOS_NURBS_INTERVAL_H_INCLUDED defined 