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

#if !defined( KRATOS_TIME_DISCRETIZATION_H_INCLUDED )
#define  KRATOS_TIME_DISCRETIZATION_H_INCLUDED

// System includes
#include <array>

// External includes

// Project includes

namespace Kratos {
namespace TimeDiscretization {

class KRATOS_API(KRATOS_CORE) BDF1
{
public:
    std::array<double, 2> ComputeBDFCoefficients(double DeltaTime) const;
};

class KRATOS_API(KRATOS_CORE) BDF2
{
public:
    std::array<double, 3> ComputeBDFCoefficients(double DeltaTime,
                                                 double PreviousDeltaTime) const;
};

class KRATOS_API(KRATOS_CORE) BDF3
{
public:
    std::array<double, 4> ComputeBDFCoefficients(double DeltaTime) const;
};

class KRATOS_API(KRATOS_CORE) BDF4
{
public:
    std::array<double, 5> ComputeBDFCoefficients(double DeltaTime) const;
};

class KRATOS_API(KRATOS_CORE) BDF5
{
public:
    std::array<double, 6> ComputeBDFCoefficients(double DeltaTime) const;
};

class KRATOS_API(KRATOS_CORE) BDF6
{
public:
    std::array<double, 7> ComputeBDFCoefficients(double DeltaTime) const;
};

class Newmark
{
public:
    // Attention, when changing the defaults do this also in the python-interface!
    // pybind cannot detect them automatically!
    Newmark(const double NewmarkBeta=0.25)
        : mNewmarkBeta(NewmarkBeta) {}
    double GetBeta()  const { return mNewmarkBeta; }
    double GetGamma() const { return 0.5; }
private:
    double mNewmarkBeta;
};

class Bossak
{
public:
    // Attention, when changing the defaults do this also in the python-interface!
    // pybind cannot detect them automatically!
    Bossak(const double AlphaM=-0.3,
           const double NewmarkBeta=0.25)
        : mAlphaM(AlphaM),
          mNewmarkBeta(NewmarkBeta) {}
    double GetAlphaM() const { return mAlphaM; }
    double GetBeta()   const { return 0.5 + mAlphaM; }
    double GetGamma()  const { return mNewmarkBeta * (1+mAlphaM) * (1+mAlphaM); }
private:
    double mAlphaM;
    double mNewmarkBeta;
};

class GeneralizedAlpha
{
public:
    // Attention, when changing the defaults do this also in the python-interface!
    // pybind cannot detect them automatically!
    GeneralizedAlpha(const double AlphaM=-0.3,
                     const double AlphaF=0.0,
                     const double NewmarkBeta=0.25)
        : mAlphaM(AlphaM),
          mAlphaF(AlphaF),
          mNewmarkBeta(NewmarkBeta) {}
    double GetAlphaM() const { return mAlphaM; }
    double GetAlphaF() const { return mAlphaF; }
    double GetBeta()   const { return 0.5 + mAlphaM - mAlphaF; }
    double GetGamma()  const { return mNewmarkBeta * (1+mAlphaM-mAlphaF) * (1+mAlphaM-mAlphaF); }
private:
    double mAlphaM;
    double mAlphaF;
    double mNewmarkBeta;
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
