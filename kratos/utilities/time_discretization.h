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
#include <vector>

// External includes

// Project includes

namespace Kratos {
namespace TimeDiscretization {

class BDF1
{
public:
    void ComputeBDFCoefficients(const double DeltaTime,
                                std::vector<double>& rCoefficients) const;
};

class BDF2
{
public:
    void ComputeBDFCoefficients(const double DeltaTime,
                                std::vector<double>& rCoefficients) const;
};

class BDF3
{
public:
    void ComputeBDFCoefficients(const double DeltaTime,
                                std::vector<double>& rCoefficients) const;
};

class BDF4
{
public:
    void ComputeBDFCoefficients(const double DeltaTime,
                                std::vector<double>& rCoefficients) const;
};

class BDF5
{
public:
    void ComputeBDFCoefficients(const double DeltaTime,
                                std::vector<double>& rCoefficients) const;
};

class BDF6
{
public:
    void ComputeBDFCoefficients(const double DeltaTime,
                                std::vector<double>& rCoefficients) const;
};

class Newmark
{
public:
    double GetBeta() const { return 0.25; }
    double GetGamma() const { return 0.5; }
};

class Bossak
{
public:
    Bossak(const double AlphaM)
        : mAlphaM(AlphaM) {}
    double GetBeta() const { return 0.5 + mAlphaM; }
    double GetGamma() const { return 0.25 * (1+mAlphaM) * (1+mAlphaM); }
private:
    double mAlphaM;
};

class GeneralizedAlpha
{
public:
    GeneralizedAlpha(const double AlphaM, const double AlphaF)
        : mAlphaM(AlphaM), mAlphaF(AlphaF) {}
    double GetBeta() const { return 0.5 + mAlphaM - mAlphaF; }
    double GetGamma() const { return 0.25 * (1+mAlphaM-mAlphaF) * (1+mAlphaM-mAlphaF); }
private:
    double mAlphaM;
    double mAlphaF;
};

std::size_t GetMinimumBufferSize(const BDF1& rTimeDisc) { return 2;}
std::size_t GetMinimumBufferSize(const BDF2& rTimeDisc) { return 3;}
std::size_t GetMinimumBufferSize(const BDF3& rTimeDisc) { return 4;}
std::size_t GetMinimumBufferSize(const BDF4& rTimeDisc) { return 5;}
std::size_t GetMinimumBufferSize(const BDF5& rTimeDisc) { return 6;}
std::size_t GetMinimumBufferSize(const BDF6& rTimeDisc) { return 7;}

std::size_t GetMinimumBufferSize(const Newmark& rTimeDisc)          { return 2;}
std::size_t GetMinimumBufferSize(const Bossak& rTimeDisc)           { return 2;}
std::size_t GetMinimumBufferSize(const GeneralizedAlpha& rTimeDisc) { return 2;}


} // namespace TimeDiscretization.
}  // namespace Kratos.

#endif // KRATOS_TIME_DISCRETIZATION_H_INCLUDED  defined
