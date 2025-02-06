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
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_constitutive/coulomb_yield_function.hpp"
#include "custom_constitutive/tension_cutoff_function.hpp"

namespace Kratos
{

MohrCoulombConstitutiveLaw::MohrCoulombConstitutiveLaw() = default;

int MohrCoulombConstitutiveLaw::Check(const Properties&   rMaterialProperties,
                                 const GeometryType& rElementGeometry,
                                 const ProcessInfo&  rCurrentProcessInfo) const
{
    // Verify Properties variables
    if (!rMaterialProperties.Has(GEO_COHESION) || rMaterialProperties[GEO_COHESION] <= 0.0)
        KRATOS_ERROR << "GEO_COHESION is not defined or has an invalid value for property: "
        << rMaterialProperties.Id() << std::endl;

    if (!rMaterialProperties.Has(GEO_FRICTION_ANGLE) || rMaterialProperties[GEO_FRICTION_ANGLE] <= 0.0)
        KRATOS_ERROR << "GEO_FRICTION_ANGLE is not defined or has an invalid value for property: "
        << rMaterialProperties.Id() << std::endl;


    if (!rMaterialProperties.Has(GEO_DILATION_ANGLE) || rMaterialProperties[GEO_DILATION_ANGLE] <= 0.0)
        KRATOS_ERROR << "GEO_DILATION_ANGLE is not defined or has an invalid value for property: "
                     << rMaterialProperties.Id() << std::endl;

    if (!rMaterialProperties.Has(GEO_TENSION_CUTOFF) || rMaterialProperties[GEO_TENSION_CUTOFF] <= 0.0)
        KRATOS_ERROR << "GEO_TENSION_CUTOFF is not defined or has an invalid value for property: "
                     << rMaterialProperties.Id() << std::endl;
    

    if (!rMaterialProperties.Has(YOUNG_MODULUS) || rMaterialProperties[YOUNG_MODULUS] <= 0.0)
        KRATOS_ERROR << "YOUNG_MODULUS has Key zero, is not defined or has an invalid value for property"
        << rMaterialProperties.Id() << std::endl;

    if (!rMaterialProperties.Has(POISSON_RATIO) || rMaterialProperties[POISSON_RATIO] < 0.0)
        KRATOS_ERROR << "POISSON_RATIO is not defined or has an invalid value for property"
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

// ================================================================================================
void MohrCoulombConstitutiveLaw::CalculateMohrCoulomb(const Properties& rProp, Vector& rCautchyStressVector)
{
    Vector principalStress = StressStrainUtilities::CalculatePrincipalStresses(rCautchyStressVector);

    double friction_angle = rProp[GEO_FRICTION_ANGLE] * Globals::Pi / 180.0;
    double cohesion       = rProp[GEO_COHESION];
    double tension_cutoff = rProp[GEO_TENSION_CUTOFF];

    auto const coulombYieldFunction = CoulombYieldFunction(friction_angle, cohesion);
    auto const tensionCutoffFunction = TensionCutoffFunction(tension_cutoff);

    double fme = coulombYieldFunction(principalStress);
    double fte = tensionCutoffFunction(principalStress);

    int region_index = this->FindRegionIndex(fme, fte);



}

// ================================================================================================
int MohrCoulombConstitutiveLaw::FindRegionIndex(double fme, double fte) 
{
    if (fme <= 0.0 && fte <= 0.0) return 0;
    if (fme < 0.0 && fte <= 0.0) return 1;


}

} // Namespace Kratos