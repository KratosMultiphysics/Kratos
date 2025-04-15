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
//                   Anne van de Graaf
//

#include "integration_coefficients_calculator.h"

namespace Kratos
{

IntegrationCoefficientsCalculator::IntegrationCoefficientsCalculator(std::unique_ptr<IntegrationCoefficientModifier> CoefficientModifier)
    : mCoefficientModifier{std::move(CoefficientModifier)}
{
}

void IntegrationCoefficientsCalculator::save(const Serializer&) const
{
    // No data members to be saved (yet)
}

void IntegrationCoefficientsCalculator::load(const Serializer&) const
{
    // No data members to be loaded (yet)
}

std::unique_ptr<IntegrationCoefficientModifier> IntegrationCoefficientsCalculator::CloneModifier() const
{
    return mCoefficientModifier ? mCoefficientModifier->Clone() : nullptr;
};

} // namespace Kratos