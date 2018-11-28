//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined(KRATOS_TIME_DISCRETIZATION_H_INCLUDED )
#define  KRATOS_TIME_DISCRETIZATION_H_INCLUDED

// System includes
#include <array>

// External includes

// Project includes

namespace Kratos {
namespace TimeDiscretization {

class BDF1
{
public:
    std::array<double, 2> ComputeBDFCoefficients(double DeltaTime) const;
};

class BDF2
{
public:
    std::array<double, 3> ComputeBDFCoefficients(double DeltaTime) const;

    std::array<double, 3> ComputeBDFCoefficients(double DeltaTime,
                                                 double PreviousDeltaTime) const;
};

class BDF3
{
public:
    std::array<double, 4> ComputeBDFCoefficients(double DeltaTime) const;
};

class BDF4
{
public:
    std::array<double, 5> ComputeBDFCoefficients(double DeltaTime) const;
};

class BDF5
{
public:
    std::array<double, 6> ComputeBDFCoefficients(double DeltaTime) const;
};

class BDF6
{
public:
    std::array<double, 7> ComputeBDFCoefficients(double DeltaTime) const;
};

class Newmark
{
public:
    double GetBeta()  const { return 0.25; }
    double GetGamma() const { return 0.5; }
};

class Bossak
{
public:
    Bossak(const double AlphaM)
        : mAlphaM(AlphaM) {}
    double GetAlphaM() const { return mAlphaM; }
    double GetBeta()   const { return 0.5 + mAlphaM; }
    double GetGamma()  const { return 0.25 * (1+mAlphaM) * (1+mAlphaM); }
private:
    double mAlphaM;
};

class GeneralizedAlpha
{
public:
    GeneralizedAlpha(const double AlphaM, const double AlphaF)
        : mAlphaM(AlphaM), mAlphaF(AlphaF) {}
    double GetAlphaM() const { return mAlphaM; }
    double GetAlphaF() const { return mAlphaF; }
    double GetBeta()   const { return 0.5 + mAlphaM - mAlphaF; }
    double GetGamma()  const { return 0.25 * (1+mAlphaM-mAlphaF) * (1+mAlphaM-mAlphaF); }
private:
    double mAlphaM;
    double mAlphaF;
};

std::size_t GetMinimumBufferSize(const BDF1& rTimeDiscr) { return 2;}
std::size_t GetMinimumBufferSize(const BDF2& rTimeDiscr) { return 3;}
std::size_t GetMinimumBufferSize(const BDF3& rTimeDiscr) { return 4;}
std::size_t GetMinimumBufferSize(const BDF4& rTimeDiscr) { return 5;}
std::size_t GetMinimumBufferSize(const BDF5& rTimeDiscr) { return 6;}
std::size_t GetMinimumBufferSize(const BDF6& rTimeDiscr) { return 7;}

std::size_t GetMinimumBufferSize(const Newmark& rTimeDiscr)          { return 2;}
std::size_t GetMinimumBufferSize(const Bossak& rTimeDiscr)           { return 2;}
std::size_t GetMinimumBufferSize(const GeneralizedAlpha& rTimeDiscr) { return 2;}

} // namespace TimeDiscretization.
}  // namespace Kratos.

#endif // KRATOS_TIME_DISCRETIZATION_H_INCLUDED  defined
