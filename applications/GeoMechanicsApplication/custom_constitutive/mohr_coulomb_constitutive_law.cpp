// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Wijtze Pieter Kikstra,
//                   Mohamed Nabi
//

// Application includes
#include "custom_constitutive/mohr_coulomb_constitutive_law.hpp"
#include "utilities/math_utils.h"

namespace Kratos
{

MohrCoulombConstitutiveLaw::MohrCoulombConstitutiveLaw() = default;

int MohrCoulombConstitutiveLaw::Check(const Properties&   rMaterialProperties,
                                 const GeometryType& rElementGeometry,
                                 const ProcessInfo&  rCurrentProcessInfo) const
{
    // Verify Properties variables
    if (!rMaterialProperties.Has(CRITICAL_DISPLACEMENT) || rMaterialProperties[CRITICAL_DISPLACEMENT] <= 0.0)
        KRATOS_ERROR << "CRITICAL_DISPLACEMENT has Key zero, is not defined or has an invalid "
                        "value for property: "
                     << rMaterialProperties.Id() << std::endl;

    if (!rMaterialProperties.Has(YOUNG_MODULUS) || rMaterialProperties[YOUNG_MODULUS] <= 0.00)
        KRATOS_ERROR
            << "YOUNG_MODULUS has Key zero, is not defined or has an invalid value for property"
            << rMaterialProperties.Id() << std::endl;

    if (!rMaterialProperties.Has(YIELD_STRESS) || rMaterialProperties[YIELD_STRESS] < 0.0)
        KRATOS_ERROR
            << "YIELD_STRESS has Key zero, is not defined or has an invalid value for property"
            << rMaterialProperties.Id() << std::endl;

    if (!rMaterialProperties.Has(FRICTION_COEFFICIENT) || rMaterialProperties[FRICTION_COEFFICIENT] < 0.0)
        KRATOS_ERROR << "FRICTION_COEFFICIENT has Key zero, is not defined or has an invalid value "
                        "for property"
                     << rMaterialProperties.Id() << std::endl;

    if (!rMaterialProperties.Has(DAMAGE_THRESHOLD) ||
        rMaterialProperties[DAMAGE_THRESHOLD] <= 0.0 || rMaterialProperties[DAMAGE_THRESHOLD] > 1.0)
        KRATOS_ERROR
            << "DAMAGE_THRESHOLD has Key zero, is not defined or has an invalid value for property"
            << rMaterialProperties.Id() << std::endl;

    return 0;
}

ConstitutiveLaw::Pointer MohrCoulombConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<MohrCoulombConstitutiveLaw>(*this);
}

double& MohrCoulombConstitutiveLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    if (rThisVariable == DAMAGE_VARIABLE || rThisVariable == STATE_VARIABLE)
        rValue = mStateVariable;
    return rValue;
}

void MohrCoulombConstitutiveLaw::SetValue(const Variable<double>& rThisVariable,
                                     const double&           rValue,
                                     const ProcessInfo&      rCurrentProcessInfo)
{
    if (rThisVariable == STATE_VARIABLE) mStateVariable = rValue;
}

double MohrCoulombConstitutiveLaw::CalculateCoulombYieldFunction(Vector& principalStress)
{
    double result = 0.0;
        //0.25 * (2.0 * principalStress(0) - principalStress(1) - principalStress(2)) * (1.0 - sPhi) - cCosPhi;


    return result;
}
} // Namespace Kratos