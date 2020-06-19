//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Riccardo Rossi, Alessandro Franci
//  Collaborators:  Massimiliano Zecchetto
//
//-------------------------------------------------------------
//

// System includes
#include <iostream>

// External includes
#include <cmath>

// Project includes
#include "custom_constitutive/fluid_laws/bingham_2D_law.h"
#include "includes/checks.h"
#include "includes/properties.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos {

//********************************CONSTRUCTOR*********************************
//****************************************************************************

Bingham2DLaw::Bingham2DLaw() : PfemFluidConstitutiveLaw() {}

//******************************COPY CONSTRUCTOR******************************
//****************************************************************************

Bingham2DLaw::Bingham2DLaw(const Bingham2DLaw& rOther) : PfemFluidConstitutiveLaw(rOther) {}

//***********************************CLONE************************************
//****************************************************************************

ConstitutiveLaw::Pointer Bingham2DLaw::Clone() const { return Kratos::make_shared<Bingham2DLaw>(*this); }

//*********************************DESTRUCTOR*********************************
//****************************************************************************

Bingham2DLaw::~Bingham2DLaw() {}

ConstitutiveLaw::SizeType Bingham2DLaw::WorkingSpaceDimension() { return 2; }

ConstitutiveLaw::SizeType Bingham2DLaw::GetStrainSize() { return 3; }

void Bingham2DLaw::CalculateMaterialResponseCauchy(Parameters& rValues) {

    Flags& r_options = rValues.GetOptions();

    const Properties& r_properties = rValues.GetMaterialProperties();

    Vector& r_strain_vector = rValues.GetStrainVector();
    Vector& r_stress_vector = rValues.GetStressVector();

    const double dynamic_viscosity = this->GetEffectiveDynamicViscosity(rValues);
    const double yield_shear = this->GetEffectiveYieldShear(rValues);
    const double adaptive_exponent = r_properties[ADAPTIVE_EXPONENT];
    double effective_dynamic_viscosity;

    const double equivalent_strain_rate =
        std::sqrt(2.0 * r_strain_vector[0] * r_strain_vector[0] + 2.0 * r_strain_vector[1] * r_strain_vector[1] +
                  4.0 * r_strain_vector[2] * r_strain_vector[2]);

    // Ensuring that the case of equivalent_strain_rate = 0 is not problematic
    const double tolerance=1e-8;
    if (equivalent_strain_rate < tolerance) {
        effective_dynamic_viscosity = yield_shear * adaptive_exponent;
    } else {
        double regularization = 1.0 - std::exp(-adaptive_exponent * equivalent_strain_rate);
        effective_dynamic_viscosity = dynamic_viscosity + regularization * yield_shear / equivalent_strain_rate;
    }

    const double strain_trace = r_strain_vector[0] + r_strain_vector[1];

    r_stress_vector[0] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[0] - strain_trace / 3.0);
    r_stress_vector[1] = 2.0 * effective_dynamic_viscosity * (r_strain_vector[1] - strain_trace / 3.0);
    r_stress_vector[2] = 2.0 * effective_dynamic_viscosity * r_strain_vector[2];

    if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        this->EffectiveViscousConstitutiveMatrix2D(effective_dynamic_viscosity, rValues.GetConstitutiveMatrix());
    }
}

std::string Bingham2DLaw::Info() const { return "Bingham2DLaw"; }

//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW******************
//*****************************************************************************

int Bingham2DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                        const ProcessInfo& rCurrentProcessInfo) {

    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_VISCOSITY);
    KRATOS_CHECK_VARIABLE_KEY(YIELD_SHEAR);
    KRATOS_CHECK_VARIABLE_KEY(ADAPTIVE_EXPONENT);
    KRATOS_CHECK_VARIABLE_KEY(BULK_MODULUS);

    if (rMaterialProperties[DYNAMIC_VISCOSITY] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for Bingham2DLaw: "
                     << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;
    }

    if (rMaterialProperties[YIELD_SHEAR] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing YIELD_SHEAR provided in process info for Bingham2DLaw: "
                     << rMaterialProperties[YIELD_SHEAR] << std::endl;
    }

    if (rMaterialProperties[ADAPTIVE_EXPONENT] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing ADAPTIVE_EXPONENT provided in process info for Bingham2DLaw: "
                     << rMaterialProperties[ADAPTIVE_EXPONENT] << std::endl;
    }

    if (rMaterialProperties[BULK_MODULUS] <= 0.0) {
        KRATOS_ERROR << "Incorrect or missing BULK_MODULUS provided in process info for Bingham2DLaw: "
                     << rMaterialProperties[BULK_MODULUS] << std::endl;
    }

    return 0;
}

double Bingham2DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    return rParameters.GetConstitutiveMatrix()(2, 2);
}

double Bingham2DLaw::GetEffectiveDensity(ConstitutiveLaw::Parameters& rParameters) const {
    return rParameters.GetMaterialProperties()[DENSITY];
}

double Bingham2DLaw::GetEffectiveDynamicViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    return rParameters.GetMaterialProperties()[DYNAMIC_VISCOSITY];
}

double Bingham2DLaw::GetEffectiveYieldShear(ConstitutiveLaw::Parameters& rParameters) const {
    return rParameters.GetMaterialProperties()[YIELD_SHEAR];
}

void Bingham2DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
}

void Bingham2DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PfemFluidConstitutiveLaw)
}

}  // Namespace Kratos
