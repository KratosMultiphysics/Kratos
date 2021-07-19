//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Ruben Zorilla
//  Collaborators:  Massimiliano Zecchetto
//
//-------------------------------------------------------------
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/fluid_laws/newtonian_2D_law.h"
#include "includes/checks.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos {

//********************************CONSTRUCTOR*********************************
//****************************************************************************

Newtonian2DLaw::Newtonian2DLaw() : PfemFluidConstitutiveLaw() {}

//******************************COPY CONSTRUCTOR******************************
//****************************************************************************

Newtonian2DLaw::Newtonian2DLaw(const Newtonian2DLaw& rOther) : PfemFluidConstitutiveLaw(rOther) {}

//***********************************CLONE************************************
//****************************************************************************

ConstitutiveLaw::Pointer Newtonian2DLaw::Clone() const { return Kratos::make_shared<Newtonian2DLaw>(*this); }

//*********************************DESTRUCTOR*********************************
//****************************************************************************

Newtonian2DLaw::~Newtonian2DLaw() {}

ConstitutiveLaw::SizeType Newtonian2DLaw::WorkingSpaceDimension() { return 2; }

ConstitutiveLaw::SizeType Newtonian2DLaw::GetStrainSize() { return 3; }

void Newtonian2DLaw::CalculateMaterialResponseCauchy(Parameters& rValues) {

    const Flags& r_options = rValues.GetOptions();
    const Vector& r_strain_vector = rValues.GetStrainVector();
    Vector& r_stress_vector = rValues.GetStressVector();

    const double effective_dynamic_viscosity = this->GetEffectiveViscosity(rValues);

    const double strain_trace = r_strain_vector[0] + r_strain_vector[1];

    r_stress_vector[0] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[0] - strain_trace / 3.0);
    r_stress_vector[1] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[1] - strain_trace / 3.0);
    r_stress_vector[2] = 2.0 * effective_dynamic_viscosity * r_strain_vector[2];

    if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        this->EffectiveViscousConstitutiveMatrix2D(effective_dynamic_viscosity, rValues.GetConstitutiveMatrix());
    }
}

int Newtonian2DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                          const ProcessInfo& rCurrentProcessInfo) {

    KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.0)
        << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for Newtonian2DLaw: "
        << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[BULK_MODULUS] <= 0.0)
        << "Incorrect or missing BULK_MODULUS provided in process info for Newtonian3DLaw: "
        << rMaterialProperties[BULK_MODULUS] << std::endl;

    return 0;
}

std::string Newtonian2DLaw::Info() const { return "Newtonian2DLaw"; }

double Newtonian2DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    return rParameters.GetMaterialProperties()[DYNAMIC_VISCOSITY];
}

double Newtonian2DLaw::GetEffectiveDensity(ConstitutiveLaw::Parameters& rParameters) const {
    return rParameters.GetMaterialProperties()[DENSITY];
}

void Newtonian2DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
}

void Newtonian2DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
}

}  // Namespace Kratos
