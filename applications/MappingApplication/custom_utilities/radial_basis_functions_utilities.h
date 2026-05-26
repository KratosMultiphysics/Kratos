//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti

#pragma once

// System includes

// External includes

// Project includes
#include "custom_utilities/mapper_utilities.h"

namespace Kratos
{
namespace RadialBasisFunctionsUtilities
{
using SizeType = std::size_t;
using IndexType = std::size_t;

enum class RBFType {
    InverseMultiquadric,
    Multiquadric,
    Gaussian,
    ThinPlateSpline,
    WendlandC2
};

/// Inverse Multiquadric RBF
struct InverseMultiquadric{
    double h;                      
    double operator()(double r) const {
        const double q = h * r;
        return 1.0 / std::sqrt(1.0 + q * q);
    }
};

/// Multiquadric RBF
struct Multiquadric{
    double h;                      
    double operator()(double r) const {
        const double q = r / h;
        return std::sqrt(1.0 + q * q);
    }
};

/// Gaussian RBF
struct Gaussian{
    double h;                      
    double operator()(double r) const {
        const double q = r / h;
        return std::exp(-0.5 * q * q);
    }
};

/// Thin Plate Spline RBF
struct ThinPlateSpline {
    double operator()(double r) const {
        if (r < 1.0e-12)
            return 0.0; 
        const double r2 = r * r;
        return r2 * std::log(r2);
    }
};

/// Wendland C2 RBF
struct WendlandC2 {
    double h;                      
    double operator()(double r) const {
        const double q = r / h;
        if (q >= 1.0)
            return 0.0;
        return std::pow(1.0 - q, 4) * (4.0 * q + 1.0); // (1-q)^4 * (4q+1)
    }
};

double KRATOS_API(MAPPING_APPLICATION) CalculateWendlandC2SupportRadius(const Matrix& rPoints, const double k = 2.5);

}  // namespace RadialBasisFunctionsUtilities.

}  // namespace Kratos.