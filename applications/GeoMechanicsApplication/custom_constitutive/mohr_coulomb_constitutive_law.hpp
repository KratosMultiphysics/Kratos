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
// //

#pragma once

// System includes
#include <cmath>

// Project includes
#include "includes/constitutive_law.h"
#include "includes/serializer.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) MohrCoulombConstitutiveLaw : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MohrCoulombConstitutiveLaw);

    MohrCoulombConstitutiveLaw();

    int Check(const Properties&   rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo&  rCurrentProcessInfo) const override;

    ConstitutiveLaw::Pointer Clone() const override;


    using ConstitutiveLaw::GetValue;
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    using ConstitutiveLaw::SetValue;
    void SetValue(const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMohrCoulomb(const Properties& rProp, Vector& rCautchyStressVector);
    [[nodiscard]] Vector CalculatePrincipalStresses(Vector& CauchyStressVector);

protected:
    // Member Variables
    double mStateVariable;

    double CalculateCoulombYieldFunction(Vector& principalStress, double phi, double cohesion);
    double CalculateTensionYieldFunction(Vector& principalStress, double tensionCutOff);

    
    
    private:
    // Serialization
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

}; // Class BilinearCohesive3DLaw

} // namespace Kratos