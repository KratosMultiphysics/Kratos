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
    Vector principalStress = this->CalculatePrincipalStresses(rCautchyStressVector);
    double phi             = rProp[GEO_FRICTION_ANGLE];
    double cohesion        = rProp[GEO_COHESION];
    double psi             = rProp[GEO_DILATION_ANGLE];
    double tensionCutOff   = rProp[GEO_TENSION_CUTOFF];
    
    double fme = this->CalculateCoulombYieldFunction(principalStress, phi, cohesion);
    double fte = this->CalculateTensionYieldFunction(principalStress, tensionCutOff);




}

// ================================================================================================
Vector MohrCoulombConstitutiveLaw::CalculatePrincipalStresses(Vector& rCauchyStressVector)
{
    auto stress_tensor = MathUtils<double>::StressVectorToTensor(rCauchyStressVector);
    Matrix PrincipalStressMatrix;
    Matrix EigenVectorsMatrix;
    MathUtils<double>::GaussSeidelEigenSystem(stress_tensor, EigenVectorsMatrix,
                                              PrincipalStressMatrix, 1.0e-16, 20);
    Vector result = ZeroVector(3);
    for (int i = 0; i < 3; ++i) {
        result(i) = PrincipalStressMatrix(i, i);
    }
    return result;
}


// ================================================================================================
double MohrCoulombConstitutiveLaw::CalculateCoulombYieldFunction(Vector& rPrincipalStress, double phi, double cohesion)
{
    double result = 0.5 * (rPrincipalStress(0) - rPrincipalStress(2)) -
                    0.5 * (rPrincipalStress(0) + rPrincipalStress(2)) * std::sin(phi) -
                    cohesion * std::cos(phi);
    return result;
}

// ================================================================================================
double MohrCoulombConstitutiveLaw::CalculateTensionYieldFunction(Vector& principalStress, double tensionCutOff)
{
    double result = tensionCutOff - principalStress(2);
    return result;
}


} // Namespace Kratos