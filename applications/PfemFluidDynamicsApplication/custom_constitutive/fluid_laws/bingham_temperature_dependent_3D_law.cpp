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
#include <cmath>

// Project includes
#include "custom_constitutive/fluid_laws/bingham_temperature_dependent_3D_law.h"
#include "includes/checks.h"
#include "includes/properties.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos {

//********************************CONSTRUCTOR*********************************
//****************************************************************************

BinghamTemperatureDependent3DLaw::BinghamTemperatureDependent3DLaw() : Bingham3DLaw() {}

//******************************COPY CONSTRUCTOR******************************
//****************************************************************************

BinghamTemperatureDependent3DLaw::BinghamTemperatureDependent3DLaw(const BinghamTemperatureDependent3DLaw& rOther) : Bingham3DLaw(rOther) {}

//***********************************CLONE************************************
//****************************************************************************

ConstitutiveLaw::Pointer BinghamTemperatureDependent3DLaw::Clone() const { return Kratos::make_shared<BinghamTemperatureDependent3DLaw>(*this); }

//*********************************DESTRUCTOR*********************************
//****************************************************************************

BinghamTemperatureDependent3DLaw::~BinghamTemperatureDependent3DLaw() {}

std::string BinghamTemperatureDependent3DLaw::Info() const { return "BinghamTemperatureDependent3DLaw"; }

//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW******************
//*****************************************************************************

int BinghamTemperatureDependent3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                        const ProcessInfo& rCurrentProcessInfo) {

    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_VISCOSITY);
    KRATOS_CHECK_VARIABLE_KEY(YIELD_SHEAR);
    KRATOS_CHECK_VARIABLE_KEY(ADAPTIVE_EXPONENT);
    KRATOS_CHECK_VARIABLE_KEY(BULK_MODULUS);

    if (rMaterialProperties[DYNAMIC_VISCOSITY] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for BinghamTemperatureDependent3DLaw: "
                     << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;
    }

    if (rMaterialProperties[YIELD_SHEAR] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing YIELD_SHEAR provided in process info for BinghamTemperatureDependent3DLaw: "
                     << rMaterialProperties[YIELD_SHEAR] << std::endl;
    }

    if (rMaterialProperties[ADAPTIVE_EXPONENT] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing ADAPTIVE_EXPONENT provided in process info for BinghamTemperatureDependent3DLaw: "
                     << rMaterialProperties[ADAPTIVE_EXPONENT] << std::endl;
    }

    if (rMaterialProperties[BULK_MODULUS] <= 0.0) {
        KRATOS_ERROR << "Incorrect or missing BULK_MODULUS provided in process info for BinghamTemperatureDependent3DLaw: "
                     << rMaterialProperties[BULK_MODULUS] << std::endl;
    }

    return 0;
}

double BinghamTemperatureDependent3DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    return rParameters.GetConstitutiveMatrix()(2, 2);
}

double BinghamTemperatureDependent3DLaw::GetEffectiveDensity(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    double effective_density;
    if (r_properties.HasTable(TEMPERATURE, DENSITY)) {
        effective_density = this->GetValueFromTable(TEMPERATURE, DENSITY, rParameters);
    } else {
        effective_density = r_properties[DENSITY];
    }
    return effective_density;
}

double BinghamTemperatureDependent3DLaw::GetEffectiveDynamicViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    double effective_viscosity;
    if (r_properties.HasTable(TEMPERATURE, DYNAMIC_VISCOSITY)) {
        effective_viscosity = this->GetValueFromTable(TEMPERATURE, DYNAMIC_VISCOSITY, rParameters);
    } else {
        effective_viscosity = r_properties[DYNAMIC_VISCOSITY];
    }
    return effective_viscosity;
}

double BinghamTemperatureDependent3DLaw::GetEffectiveYieldShear(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    double effective_yield_shear;
    if (r_properties.HasTable(TEMPERATURE, YIELD_SHEAR)) {
        effective_yield_shear = this->GetValueFromTable(TEMPERATURE, YIELD_SHEAR, rParameters);
    } else {
        effective_yield_shear = r_properties[YIELD_SHEAR];
    }
    return effective_yield_shear;
}

void BinghamTemperatureDependent3DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Bingham3DLaw)
}

void BinghamTemperatureDependent3DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Bingham3DLaw)
}

}  // Namespace Kratos
