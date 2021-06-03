//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Alessandro Franci
//  Collaborators:
//
//-------------------------------------------------------------
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/solid_laws/hypoelastic_2D_law.h"
#include "includes/checks.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos {

//********************************CONSTRUCTOR*********************************
//****************************************************************************

Hypoelastic2DLaw::Hypoelastic2DLaw() : PfemSolidConstitutiveLaw() {}

//******************************COPY CONSTRUCTOR******************************
//****************************************************************************

Hypoelastic2DLaw::Hypoelastic2DLaw(const Hypoelastic2DLaw& rOther) : PfemSolidConstitutiveLaw(rOther) {}

//***********************************CLONE************************************
//****************************************************************************

ConstitutiveLaw::Pointer Hypoelastic2DLaw::Clone() const { return Kratos::make_shared<Hypoelastic2DLaw>(*this); }

//*********************************DESTRUCTOR*********************************
//****************************************************************************

Hypoelastic2DLaw::~Hypoelastic2DLaw() {}

ConstitutiveLaw::SizeType Hypoelastic2DLaw::WorkingSpaceDimension() { return 2; }

ConstitutiveLaw::SizeType Hypoelastic2DLaw::GetStrainSize() { return 3; }

void Hypoelastic2DLaw::CalculateMaterialResponseCauchy(Parameters& rParameters) {

    const Vector& r_strain_vector = rParameters.GetStrainVector();
    Vector& r_stress_vector = rParameters.GetStressVector();

    const double young_modulus = this->GetEffectiveYoungModulus(rParameters);
    const double poisson_ratio = this->GetEffectivePoissonRatio(rParameters);
    const double time_step = rParameters.GetProcessInfo()[DELTA_TIME];

    const double second_lame = time_step * young_modulus / (2.0 * (1.0 + poisson_ratio));

    const double strain_trace = r_strain_vector[0] + r_strain_vector[1];

    r_stress_vector[0] += 2.0 * second_lame * (r_strain_vector[0] - strain_trace / 3.0);
    r_stress_vector[1] += 2.0 * second_lame * (r_strain_vector[1] - strain_trace / 3.0);
    r_stress_vector[2] += 2.0 * second_lame * r_strain_vector[2];
}

int Hypoelastic2DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                            const ProcessInfo& rCurrentProcessInfo) {

    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0)
        << "Incorrect or missing YOUNG_MODULUS provided in process info for Hypoelastic2DLaw: "
        << rMaterialProperties[YOUNG_MODULUS] << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[POISSON_RATIO] < 0.0)
        << "Incorrect or missing POISSON_RATIO provided in process info for Hypoelastic2DLaw: "
        << rMaterialProperties[POISSON_RATIO] << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0)
        << "Incorrect or missing DENSITY provided in process info for Hypoelastic2DLaw: "
        << rMaterialProperties[DENSITY] << std::endl;

    return 0;
}

std::string Hypoelastic2DLaw::Info() const { return "Hypoelastic2DLaw"; }

double Hypoelastic2DLaw::GetEffectiveYoungModulus(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    const double effective_young_modulus = r_properties[YOUNG_MODULUS];
    return effective_young_modulus;
}

double Hypoelastic2DLaw::GetEffectivePoissonRatio(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    const double effective_poisson_ratio = r_properties[POISSON_RATIO];
    return effective_poisson_ratio;
}

double Hypoelastic2DLaw::GetEffectiveDensity(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    const double effective_density = r_properties[DENSITY];
    return effective_density;
}

void Hypoelastic2DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PfemSolidConstitutiveLaw)
}

void Hypoelastic2DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PfemSolidConstitutiveLaw)
}

}  // Namespace Kratos
