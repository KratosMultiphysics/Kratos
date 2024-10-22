// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#include "compressibility_calculator.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "custom_retention/retention_law_factory.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

CompressibilityCalculator::CompressibilityCalculator(const InputProvider& rInputProvider)
    : mInputProvider(rInputProvider)
{
}

Matrix CompressibilityCalculator::LHSContribution() { return {}; }

Vector CompressibilityCalculator::RHSContribution() { return {}; }

void CompressibilityCalculator::CalculateLeftAndRightHandSide(Matrix& rLeftHandSideMatrix, Vector& rRightHandSideVector)
{
    auto compressibility_matrix = CalculateCompressibilityMatrix(mInputProvider.GetNContainer(),
                                                                 mInputProvider.GetIntegrationCoefficients());
    rLeftHandSideMatrix += compressibility_matrix * mInputProvider.GetDtPressureCoefficient();
    rRightHandSideVector += -prod(compressibility_matrix, mInputProvider.GetNodalValues(DT_WATER_PRESSURE));
}

Matrix CompressibilityCalculator::CalculateCompressibilityMatrix(const Matrix& rNContainer,
                                                                 const Vector& rIntegrationCoefficients)
{
    const auto&              r_properties = mInputProvider.GetElementProperties();
    RetentionLaw::Parameters parameters(r_properties);
    auto                     retention_law = RetentionLawFactory::Clone(r_properties);

    auto result = Matrix{ZeroMatrix{rNContainer.size2(), rNContainer.size2()}};
    for (unsigned int integration_point_index = 0;
         integration_point_index < rIntegrationCoefficients.size(); ++integration_point_index) {
        const auto   N                  = Vector{row(rNContainer, integration_point_index)};
        const double BiotModulusInverse = CalculateBiotModulusInverse(integration_point_index);
        result += GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
            N, BiotModulusInverse, rIntegrationCoefficients[integration_point_index]);
    }
    return result;
}

double CompressibilityCalculator::CalculateBiotModulusInverse(const unsigned int integrationPointIndex) const
{
    const auto&  r_properties     = mInputProvider.GetElementProperties();
    const double biot_coefficient = r_properties[BIOT_COEFFICIENT];

    double bulk_fluid = TINY;
    if (!r_properties[IGNORE_UNDRAINED]) {
        bulk_fluid = r_properties[BULK_MODULUS_FLUID];
    }
    double result = (biot_coefficient - r_properties[POROSITY]) / r_properties[BULK_MODULUS_SOLID] +
                    r_properties[POROSITY] / bulk_fluid;

    RetentionLaw::Parameters RetentionParameters(mInputProvider.GetElementProperties());
    const double             degree_of_saturation =
        mInputProvider.GetRetentionLaws()[integrationPointIndex]->CalculateSaturation(RetentionParameters);
    const double derivative_of_saturation =
        mInputProvider.GetRetentionLaws()[integrationPointIndex]->CalculateDerivativeOfSaturation(RetentionParameters);

    result *= degree_of_saturation;
    result -= derivative_of_saturation * r_properties[POROSITY];
    return result;
}

} // namespace Kratos
