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

CalculateIntegrationCoefficients0::CalculateIntegrationCoefficients0(std::unique_ptr<IntegrationCoefficientModifier> CoefficientModifier)
    : mCoefficientModifier{std::move(CoefficientModifier)}
{
}

std::unique_ptr<IntegrationCoefficientModifier> CalculateIntegrationCoefficients0::CloneModifier() const
{
    if (mCoefficientModifier) {
        return mCoefficientModifier->Clone();
    } else {
        return nullptr;
    }
};

} // namespace Kratos