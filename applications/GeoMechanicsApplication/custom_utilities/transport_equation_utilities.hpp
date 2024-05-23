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
        return -PORE_PRESSURE_SIGN_FACTOR * DynamicViscosityInverse *
               prod(rGradNpT, Matrix(prod(rMaterialPermeabilityMatrix, trans(rGradNpT)))) *
               RelativePermeability * PermeabilityUpdateFactor * IntegrationCoefficient;
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
                                             const Properties&  rProp,
                                             const ProcessInfo& rCurrentProcessInfo)
    {
        RetentionLaw::Parameters retention_parameters(rProp, rCurrentProcessInfo);
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
}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
