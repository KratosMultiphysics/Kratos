//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Massimiliano Zecchetto
//  Collaborators:
//
//-------------------------------------------------------------
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/solid_laws/hypoelastic_temperature_dependent_2D_law.h"
#include "includes/checks.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos {

//********************************CONSTRUCTOR*********************************
//****************************************************************************

HypoelasticTemperatureDependent2DLaw::HypoelasticTemperatureDependent2DLaw() : Hypoelastic2DLaw() {}

//******************************COPY CONSTRUCTOR******************************
//****************************************************************************

HypoelasticTemperatureDependent2DLaw::HypoelasticTemperatureDependent2DLaw(const HypoelasticTemperatureDependent2DLaw& rOther) : Hypoelastic2DLaw(rOther) {}

//***********************************CLONE************************************
//****************************************************************************

ConstitutiveLaw::Pointer HypoelasticTemperatureDependent2DLaw::Clone() const { return Kratos::make_shared<HypoelasticTemperatureDependent2DLaw>(*this); }

//*********************************DESTRUCTOR*********************************
//****************************************************************************

HypoelasticTemperatureDependent2DLaw::~HypoelasticTemperatureDependent2DLaw() {}

int HypoelasticTemperatureDependent2DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                            const ProcessInfo& rCurrentProcessInfo) {

    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0)
        << "Incorrect or missing YOUNG_MODULUS provided in process info for HypoelasticTemperatureDependent2DLaw: "
        << rMaterialProperties[YOUNG_MODULUS] << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[POISSON_RATIO] < 0.0)
        << "Incorrect or missing POISSON_RATIO provided in process info for HypoelasticTemperatureDependent2DLaw: "
        << rMaterialProperties[POISSON_RATIO] << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0)
        << "Incorrect or missing DENSITY provided in process info for HypoelasticTemperatureDependent2DLaw: "
        << rMaterialProperties[DENSITY] << std::endl;

    return 0;
}

std::string HypoelasticTemperatureDependent2DLaw::Info() const { return "HypoelasticTemperatureDependent2DLaw"; }

double HypoelasticTemperatureDependent2DLaw::GetEffectiveYoungModulus(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    double effective_young_modulus;
    if (r_properties.HasTable(TEMPERATURE, YOUNG_MODULUS)) {
        effective_young_modulus = this->GetValueFromTable(TEMPERATURE, YOUNG_MODULUS, rParameters);
    } else {
        effective_young_modulus = r_properties[YOUNG_MODULUS];
    }
    return effective_young_modulus;
}

double HypoelasticTemperatureDependent2DLaw::GetEffectivePoissonRatio(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    double effective_poisson_ratio;
    if (r_properties.HasTable(TEMPERATURE, POISSON_RATIO)) {
        effective_poisson_ratio = this->GetValueFromTable(TEMPERATURE, POISSON_RATIO, rParameters);
    } else {
        effective_poisson_ratio = r_properties[POISSON_RATIO];
    }
    return effective_poisson_ratio;
}

double HypoelasticTemperatureDependent2DLaw::GetEffectiveDensity(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    double effective_density;
    if (r_properties.HasTable(TEMPERATURE, DENSITY)) {
        effective_density = this->GetValueFromTable(TEMPERATURE, DENSITY, rParameters);
    } else {
        effective_density = r_properties[DENSITY];
    }
    return effective_density;
}

void HypoelasticTemperatureDependent2DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Hypoelastic2DLaw)
}

void HypoelasticTemperatureDependent2DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Hypoelastic2DLaw)
}

}  // Namespace Kratos
