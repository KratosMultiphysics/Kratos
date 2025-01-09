// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   Richard Faasse
//
#include "filter_compressibility_calculator.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

FilterCompressibilityCalculator::FilterCompressibilityCalculator(InputProvider rInputProvider)
    : mInputProvider(std::move(rInputProvider))
{
}

Matrix FilterCompressibilityCalculator::LHSContribution()
{
    return LHSContribution(CalculateCompressibilityMatrix());
}

Vector FilterCompressibilityCalculator::RHSContribution()
{
    return RHSContribution(CalculateCompressibilityMatrix());
}

Vector FilterCompressibilityCalculator::RHSContribution(const Matrix& rCompressibilityMatrix) const
{
    return -prod(rCompressibilityMatrix, mInputProvider.GetNodalValues(DT_WATER_PRESSURE));
}

Matrix FilterCompressibilityCalculator::LHSContribution(const Matrix& rCompressibilityMatrix) const
{
    return mInputProvider.GetMatrixScalarFactor() * rCompressibilityMatrix;
}

std::pair<Matrix, Vector> FilterCompressibilityCalculator::LocalSystemContribution()
{
    const auto compressibility_matrix = CalculateCompressibilityMatrix();
    return {LHSContribution(compressibility_matrix), RHSContribution(compressibility_matrix)};
}

Matrix FilterCompressibilityCalculator::CalculateCompressibilityMatrix() const
{
    const auto& r_N_container            = mInputProvider.GetNContainer();
    const auto& integration_coefficients = mInputProvider.GetIntegrationCoefficients();
    const auto& projected_gravity_on_integration_points =
        mInputProvider.GetProjectedGravityForIntegrationPoints();
    auto result = Matrix{ZeroMatrix{r_N_container.size2(), r_N_container.size2()}};
    for (unsigned int integration_point_index = 0;
         integration_point_index < integration_coefficients.size(); ++integration_point_index) {
        const auto N = Vector{row(r_N_container, integration_point_index)};
        result += GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
            N, CalculateBiotModulusInverse(projected_gravity_on_integration_points[integration_point_index]),
            integration_coefficients[integration_point_index]);
    }
    return result;
}

double FilterCompressibilityCalculator::CalculateBiotModulusInverse(double ProjectedGravity) const
{
    const auto&  r_properties = mInputProvider.GetElementProperties();
    const double result =
        1.0 / (r_properties[DENSITY_WATER] * ProjectedGravity * r_properties[PIPE_ELEMENT_LENGTH]) +
        1.0 / r_properties[BULK_MODULUS_FLUID];
    return result;
}

const Properties& FilterCompressibilityCalculator::InputProvider::GetElementProperties() const
{
    return mGetElementProperties();
}

const std::vector<RetentionLaw::Pointer>& FilterCompressibilityCalculator::InputProvider::GetRetentionLaws() const
{
    return mGetRetentionLaws();
}

const Matrix& FilterCompressibilityCalculator::InputProvider::GetNContainer() const
{
    return mGetNContainer();
}

Vector FilterCompressibilityCalculator::InputProvider::GetIntegrationCoefficients() const
{
    return mGetIntegrationCoefficients();
}

Vector FilterCompressibilityCalculator::InputProvider::GetProjectedGravityForIntegrationPoints() const
{
    return mGetProjectedGravityForIntegrationPoints();
}

double FilterCompressibilityCalculator::InputProvider::GetMatrixScalarFactor() const
{
    return mGetMatrixScalarFactor();
}

Vector FilterCompressibilityCalculator::InputProvider::GetNodalValues(const Variable<double>& rVariable) const
{
    return mGetNodalValues(rVariable);
}

} // namespace Kratos
