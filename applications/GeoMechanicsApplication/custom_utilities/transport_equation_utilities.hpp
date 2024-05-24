// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#pragma once

// Project includes

// Application includes
#include "custom_retention/retention_law.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class GeoTransportEquationUtilities
{
public:
    template <unsigned int TDim, unsigned int TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TNumNodes> CalculatePermeabilityMatrix(
        const Matrix&                            rGradNpT,
        double                                   DynamicViscosityInverse,
        const BoundedMatrix<double, TDim, TDim>& rMaterialPermeabilityMatrix,
        double                                   RelativePermeability,
        double                                   PermeabilityUpdateFactor,
        double                                   IntegrationCoefficient)
    {
        return CalculatePermeabilityMatrix(rGradNpT, DynamicViscosityInverse,
                                           rMaterialPermeabilityMatrix, RelativePermeability,
                                           PermeabilityUpdateFactor, IntegrationCoefficient);
    }

    static inline Matrix CalculatePermeabilityMatrix(const Matrix& rGradNpT,
                                                     double        DynamicViscosityInverse,
                                                     const Matrix& rMaterialPermeabilityMatrix,
                                                     double        RelativePermeability,
                                                     double        PermeabilityUpdateFactor,
                                                     double        IntegrationCoefficient)
    {
        return CalculatePermeabilityMatrix(rGradNpT, DynamicViscosityInverse, rMaterialPermeabilityMatrix,
                                           RelativePermeability * PermeabilityUpdateFactor,
                                           IntegrationCoefficient);
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TNumNodes> CalculatePermeabilityMatrix(
        const Matrix&                            rGradNpT,
        double                                   DynamicViscosityInverse,
        const BoundedMatrix<double, TDim, TDim>& rMaterialPermeabilityMatrix,
        double                                   RelativePermeability,
        double                                   IntegrationCoefficient)
    {
        return CalculatePermeabilityMatrix(rGradNpT, DynamicViscosityInverse, rMaterialPermeabilityMatrix,
                                           RelativePermeability, IntegrationCoefficient);
    }

    static inline Matrix CalculatePermeabilityMatrix(const Matrix& rGradNpT,
                                                     double        DynamicViscosityInverse,
                                                     const Matrix& rMaterialPermeabilityMatrix,
                                                     double        RelativePermeability,
                                                     double        IntegrationCoefficient)
    {
        return -PORE_PRESSURE_SIGN_FACTOR * DynamicViscosityInverse *
               prod(rGradNpT, Matrix(prod(rMaterialPermeabilityMatrix, trans(rGradNpT)))) *
               RelativePermeability * IntegrationCoefficient;
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    static inline BoundedMatrix<double, TNumNodes * TDim, TNumNodes> CalculateCouplingMatrix(
        const Matrix& rB, const Vector& rVoigtVector, const Vector& rNp, double BiotCoefficient, double BishopCoefficient, double IntegrationCoefficient)
    {
        return CalculateCouplingMatrix(rB, rVoigtVector, rNp, BiotCoefficient, BishopCoefficient,
                                       IntegrationCoefficient);
    }

    static inline Matrix CalculateCouplingMatrix(const Matrix& rB,
                                                 const Vector& rVoigtVector,
                                                 const Vector& rNp,
                                                 double        BiotCoefficient,
                                                 double        BishopCoefficient,
                                                 double        IntegrationCoefficient)
    {
        return PORE_PRESSURE_SIGN_FACTOR * BiotCoefficient * BishopCoefficient *
               outer_prod(Vector(prod(trans(rB), rVoigtVector)), rNp) * IntegrationCoefficient;
    }

    template <unsigned int TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TNumNodes> CalculateCompressibilityMatrix(
        const Vector& rNp, double BiotModulusInverse, double IntegrationCoefficient)
    {
        return CalculateCompressibilityMatrix(rNp, BiotModulusInverse, IntegrationCoefficient);
    }

    static inline Matrix CalculateCompressibilityMatrix(const Vector& rNp, double BiotModulusInverse, double IntegrationCoefficient)
    {
        return -PORE_PRESSURE_SIGN_FACTOR * BiotModulusInverse * outer_prod(rNp, rNp) * IntegrationCoefficient;
    }

    static double CalculateSoilDensity(double DegreeOfSaturation, const Properties& rProp)
    {
        return DegreeOfSaturation * rProp[POROSITY] * rProp[DENSITY_WATER] +
               (1.0 - rProp[POROSITY]) * rProp[DENSITY_SOLID];
    }

    static Vector CalculateSoilDensities(const Vector& rDegreesSaturation, const Properties& rProp)
    {
        Vector result(rDegreesSaturation.size());
        std::transform(rDegreesSaturation.cbegin(), rDegreesSaturation.cend(), result.begin(),
                       [&rProp](const auto& degree_saturation) {
            return CalculateSoilDensity(degree_saturation, rProp);
        });
        return result;
    }

    static Vector CalculateDegreesSaturation(const Vector& rPressureSolution,
                                             const Matrix& rNContainer,
                                             const std::vector<RetentionLaw::Pointer>& rRetentionLawVector,
                                             const Properties& rProp)
    {
        RetentionLaw::Parameters retention_parameters(rProp);
        Vector                   result(rRetentionLawVector.size());
        for (unsigned int g_point = 0; g_point < rRetentionLawVector.size(); ++g_point) {
            const double fluid_pressure = CalculateFluidPressure(row(rNContainer, g_point), rPressureSolution);
            retention_parameters.SetFluidPressure(fluid_pressure);
            result(g_point) = rRetentionLawVector[g_point]->CalculateSaturation(retention_parameters);
        }
        return result;
    }

    [[nodiscard]] static double CalculateFluidPressure(const Vector& rN, const Vector& rPressureVector)
    {
        return inner_prod(rN, rPressureVector);
    }

    [[nodiscard]] static double CalculateBiotModulusInverse(double BiotCoefficient,
                                                            double DegreeOfSaturation,
                                                            double DerivativeOfSaturation,
                                                            const Properties& rProperties)
    {
        const auto bulk_modulus_fluid = rProperties[IGNORE_UNDRAINED] ? TINY : rProperties[BULK_MODULUS_FLUID];
        double result = (BiotCoefficient - rProperties[POROSITY]) / rProperties[BULK_MODULUS_SOLID] +
                        rProperties[POROSITY] / bulk_modulus_fluid;

        result *= DegreeOfSaturation;
        result -= DerivativeOfSaturation * rProperties[POROSITY];

        return result;
    }

    [[nodiscard]] static double CalculateBulkModulus(const Matrix& rConstitutiveMatrix)
    {
        KRATOS_ERROR_IF(rConstitutiveMatrix.size1() == 0)
            << "Constitutive matrix is empty, aborting bulk modulus calculation.\n";
        const SizeType index_G = rConstitutiveMatrix.size1() - 1;
        return rConstitutiveMatrix(0, 0) - (4.0 / 3.0) * rConstitutiveMatrix(index_G, index_G);
    }

    [[nodiscard]] static std::vector<double> CalculateBiotCoefficients(const std::vector<Matrix>& rConstitutiveMatrices,
                                                                       const Properties& rProperties)
    {
        std::vector<double> result;
        std::transform(rConstitutiveMatrices.begin(), rConstitutiveMatrices.end(),
                       std::back_inserter(result), [&rProperties](const Matrix& rConstitutiveMatrix) {
            return CalculateBiotCoefficient(rConstitutiveMatrix, rProperties);
        });

        return result;
    }

    [[nodiscard]] static std::vector<double> CalculateFluidPressures(const Matrix& rNContainer, const Vector& rPressureVector)
    {
        auto result = std::vector<double>{};
        for (auto i = std::size_t{0}; i < rNContainer.size1(); ++i) {
            result.emplace_back(CalculateFluidPressure(row(rNContainer, i), rPressureVector));
        }
        return result;
    }

    [[nodiscard]] static double CalculatePermeabilityUpdateFactor(const Vector&     rStrainVector,
                                                                  const Properties& rProperties)
    {
        if (rProperties[PERMEABILITY_CHANGE_INVERSE_FACTOR] > 0.0) {
            const double InverseCK = rProperties[PERMEABILITY_CHANGE_INVERSE_FACTOR];
            const double epsV      = StressStrainUtilities::CalculateTrace(rStrainVector);
            const double ePrevious = rProperties[POROSITY] / (1.0 - rProperties[POROSITY]);
            const double eCurrent  = (1.0 + ePrevious) * std::exp(epsV) - 1.0;
            const double permLog10 = (eCurrent - ePrevious) * InverseCK;
            return std::pow(10.0, permLog10);
        }

        return 1.0;
    }

    [[nodiscard]] static std::vector<double> CalculatePermeabilityUpdateFactors(const std::vector<Vector>& rStrainVectors,
                                                                                const Properties& rProperties)
    {
        auto result = std::vector<double>{};
        std::transform(rStrainVectors.cbegin(), rStrainVectors.cend(), std::back_inserter(result),
                       [&rProperties](const auto& rStrainVector) {
            return CalculatePermeabilityUpdateFactor(rStrainVector, rProperties);
        });
        return result;
    }

private:
    [[nodiscard]] static double CalculateBiotCoefficient(const Matrix&     rConstitutiveMatrix,
                                                         const Properties& rProperties)
    {
        return rProperties.Has(BIOT_COEFFICIENT)
                   ? rProperties[BIOT_COEFFICIENT]
                   : 1.0 - CalculateBulkModulus(rConstitutiveMatrix) / rProperties[BULK_MODULUS_SOLID];
    }
}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
