// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "transport_equation_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"

namespace Kratos
{

Matrix GeoTransportEquationUtilities::CalculatePermeabilityMatrix(const Matrix& rGradNpT,
                                                                  double DynamicViscosityInverse,
                                                                  const Matrix& rMaterialPermeabilityMatrix,
                                                                  double RelativePermeability,
                                                                  double IntegrationCoefficient)
{
    return -PORE_PRESSURE_SIGN_FACTOR * DynamicViscosityInverse *
           prod(rGradNpT, Matrix(prod(rMaterialPermeabilityMatrix, trans(rGradNpT)))) *
           RelativePermeability * IntegrationCoefficient;
}

std::vector<double> GeoTransportEquationUtilities::CalculateFluidPressures(const Matrix& rNContainer,
                                                                           const Vector& rPressureVector)
{
    const auto          num_rows = rNContainer.size1();
    std::vector<double> result;
    result.reserve(num_rows);

    for (std::size_t i = 0; i < num_rows; ++i) {
        result.push_back(CalculateFluidPressure(row(rNContainer, i), rPressureVector));
    }

    return result;
}

std::vector<double> GeoTransportEquationUtilities::CalculateInverseBiotModuli(
    const std::vector<double>& rBiotCoefficients,
    const std::vector<double>& rDegreesOfSaturation,
    const std::vector<double>& DerivativesOfSaturation,
    const Properties&          rProperties)
{
    std::vector<double> result;
    result.reserve(rBiotCoefficients.size());
    for (std::size_t i = 0; i < rBiotCoefficients.size(); ++i) {
        result.push_back(CalculateInverseBiotModulus(rBiotCoefficients[i], rDegreesOfSaturation[i],
                                                     DerivativesOfSaturation[i], rProperties));
    }
    return result;
}

double GeoTransportEquationUtilities::CalculateBulkModulus(const Matrix& rConstitutiveMatrix)
{
    KRATOS_ERROR_IF(rConstitutiveMatrix.size1() == 0)
        << "Constitutive matrix is empty, aborting bulk modulus calculation.\n";
    const SizeType index_G = rConstitutiveMatrix.size1() - 1;
    return rConstitutiveMatrix(0, 0) - (4.0 / 3.0) * rConstitutiveMatrix(index_G, index_G);
}

std::vector<double> GeoTransportEquationUtilities::CalculateBiotCoefficients(const std::vector<Matrix>& rConstitutiveMatrices,
                                                                             const Properties& rProperties)
{
    std::vector<double> result;
    result.reserve(rConstitutiveMatrices.size());
    std::transform(rConstitutiveMatrices.begin(), rConstitutiveMatrices.end(),
                   std::back_inserter(result), [&rProperties](const Matrix& rConstitutiveMatrix) {
        return CalculateBiotCoefficient(rConstitutiveMatrix, rProperties);
    });

    return result;
}

std::vector<double> GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactors(
    const std::vector<Vector>& rStrainVectors, const Properties& rProperties)
{
    auto result = std::vector<double>{};
    result.reserve(rStrainVectors.size());
    std::transform(rStrainVectors.cbegin(), rStrainVectors.cend(), std::back_inserter(result),
                   [&rProperties](const auto& rStrainVector) {
        return CalculatePermeabilityUpdateFactor(rStrainVector, rProperties);
    });
    return result;
}

double GeoTransportEquationUtilities::CalculateBiotCoefficient(const Matrix& rConstitutiveMatrix,
                                                               const Properties& rProperties)
{
    return rProperties.Has(BIOT_COEFFICIENT)
               ? rProperties[BIOT_COEFFICIENT]
               : 1.0 - CalculateBulkModulus(rConstitutiveMatrix) / rProperties[BULK_MODULUS_SOLID];
}

double GeoTransportEquationUtilities::CalculateInverseBiotModulus(double BiotCoefficient,
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

double GeoTransportEquationUtilities::CalculatePermeabilityUpdateFactor(const Vector& rStrainVector,
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

double GeoTransportEquationUtilities::CalculateParticleDiameter(const Properties& rProperties)
{
    return rProperties.Has(PIPE_MODIFIED_D) && rProperties[PIPE_MODIFIED_D]
               ? 2.08e-4 * std::pow((rProperties[PIPE_D_70] / 2.08e-4), 0.4)
               : rProperties[PIPE_D_70];
}

} // namespace Kratos
