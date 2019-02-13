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
class ProcessInfo; // forward-declaring to not having to include it here

namespace TimeDiscretization {

class KRATOS_API(KRATOS_CORE) BDF1
{
public:
    std::array<double, 2> ComputeBDFCoefficients(double DeltaTime) const;
    std::array<double, 2> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const;
};

class KRATOS_API(KRATOS_CORE) BDF2
{
public:
    std::array<double, 3> ComputeBDFCoefficients(double DeltaTime,
                                                 double PreviousDeltaTime) const;
    std::array<double, 3> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const;
};

class KRATOS_API(KRATOS_CORE) BDF3
{
public:
    std::array<double, 4> ComputeBDFCoefficients(double DeltaTime) const;
    std::array<double, 4> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const;
};

class KRATOS_API(KRATOS_CORE) BDF4
{
public:
    std::array<double, 5> ComputeBDFCoefficients(double DeltaTime) const;
    std::array<double, 5> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const;
};

class KRATOS_API(KRATOS_CORE) BDF5
{
public:
    std::array<double, 6> ComputeBDFCoefficients(double DeltaTime) const;
    std::array<double, 6> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const;
};

class KRATOS_API(KRATOS_CORE) BDF6
{
public:
    std::array<double, 7> ComputeBDFCoefficients(double DeltaTime) const;
    std::array<double, 7> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const;
};

class Newmark
{
public:
    Newmark() = default;

    explicit Newmark(const double NewmarkBeta,
                     const double NewmarkGamma)
        : mNewmarkBeta(NewmarkBeta),
          mNewmarkGamma(NewmarkGamma) {}

    double GetBeta()  const { return mNewmarkBeta; }
    double GetGamma() const { return mNewmarkGamma; }

private:
    double mNewmarkBeta=0.25;
    double mNewmarkGamma=0.5;
};

class Bossak
{
public:
    Bossak() = default;

    explicit Bossak(const double AlphaM)
        : mAlphaM(AlphaM) {}

    explicit Bossak(const double AlphaM,
                    const double NewmarkBeta,
                    const double NewmarkGamma)
        : mAlphaM(AlphaM),
          mNewmarkBeta(NewmarkBeta),
          mNewmarkGamma(NewmarkGamma) {}

    double GetAlphaM() const { return mAlphaM; }
    double GetBeta()   const { return mNewmarkBeta * (1-mAlphaM) * (1-mAlphaM); }
    double GetGamma()  const { return mNewmarkGamma - mAlphaM; }

private:
    double mAlphaM=-0.3;
    double mNewmarkBeta=0.25;
    double mNewmarkGamma=0.5;
};

class GeneralizedAlpha
{
public:
    GeneralizedAlpha() = default;

    explicit GeneralizedAlpha(const double AlphaM,
                              const double AlphaF)
        : mAlphaM(AlphaM),
          mAlphaF(AlphaF) {}

    explicit GeneralizedAlpha(const double AlphaM,
                              const double AlphaF,
                              const double NewmarkBeta,
                              const double NewmarkGamma)
        : mAlphaM(AlphaM),
          mAlphaF(AlphaF),
          mNewmarkBeta(NewmarkBeta),
          mNewmarkGamma(NewmarkGamma) {}

    double GetAlphaM() const { return mAlphaM; }
    double GetAlphaF() const { return mAlphaF; }
    double GetBeta()   const { return mNewmarkBeta * (1-mAlphaM+mAlphaF) * (1-mAlphaM+mAlphaF); }
    double GetGamma()  const { return mNewmarkGamma - mAlphaM + mAlphaF; }

private:
    double mAlphaM=-0.3;
    double mAlphaF=0.0;
    double mNewmarkBeta=0.25;
    double mNewmarkGamma=0.5;
};

inline std::size_t GetMinimumBufferSize(const BDF1& rTimeDiscr) { return 2;}
inline std::size_t GetMinimumBufferSize(const BDF2& rTimeDiscr) { return 3;}
inline std::size_t GetMinimumBufferSize(const BDF3& rTimeDiscr) { return 4;}
inline std::size_t GetMinimumBufferSize(const BDF4& rTimeDiscr) { return 5;}
inline std::size_t GetMinimumBufferSize(const BDF5& rTimeDiscr) { return 6;}
inline std::size_t GetMinimumBufferSize(const BDF6& rTimeDiscr) { return 7;}

inline std::size_t GetMinimumBufferSize(const Newmark& rTimeDiscr)          { return 2;}
inline std::size_t GetMinimumBufferSize(const Bossak& rTimeDiscr)           { return 2;}
inline std::size_t GetMinimumBufferSize(const GeneralizedAlpha& rTimeDiscr) { return 2;}

} // namespace TimeDiscretization.
}  // namespace Kratos.

#endif // KRATOS_TIME_DISCRETIZATION_H_INCLUDED  defined
