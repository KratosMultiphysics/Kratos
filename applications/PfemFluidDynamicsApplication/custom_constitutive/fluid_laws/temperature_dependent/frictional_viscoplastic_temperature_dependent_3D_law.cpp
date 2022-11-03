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
#include <cmath>

// Project includes
#include "custom_constitutive/fluid_laws/temperature_dependent/frictional_viscoplastic_temperature_dependent_3D_law.h"
#include "includes/checks.h"
#include "includes/properties.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos {

//********************************CONSTRUCTOR*********************************
//****************************************************************************

FrictionalViscoplasticTemperatureDependent3DLaw::FrictionalViscoplasticTemperatureDependent3DLaw() : FrictionalViscoplastic3DLaw() {}

//******************************COPY CONSTRUCTOR******************************
//****************************************************************************

FrictionalViscoplasticTemperatureDependent3DLaw::FrictionalViscoplasticTemperatureDependent3DLaw(const FrictionalViscoplasticTemperatureDependent3DLaw& rOther) : FrictionalViscoplastic3DLaw(rOther) {}

//***********************************CLONE************************************
//****************************************************************************

ConstitutiveLaw::Pointer FrictionalViscoplasticTemperatureDependent3DLaw::Clone() const { return Kratos::make_shared<FrictionalViscoplasticTemperatureDependent3DLaw>(*this); }

//*********************************DESTRUCTOR*********************************
//****************************************************************************

FrictionalViscoplasticTemperatureDependent3DLaw::~FrictionalViscoplasticTemperatureDependent3DLaw() {}

std::string FrictionalViscoplasticTemperatureDependent3DLaw::Info() const { return "FrictionalViscoplasticTemperatureDependent3DLaw"; }

//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW******************
//*****************************************************************************

int FrictionalViscoplasticTemperatureDependent3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                        const ProcessInfo& rCurrentProcessInfo) const {

    if (rMaterialProperties[DYNAMIC_VISCOSITY] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for FrictionalViscoplasticTemperatureDependent3DLaw: "
                     << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;
    }

    if (rMaterialProperties[INTERNAL_FRICTION_ANGLE] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing INTERNAL_FRICTION_ANGLE provided in process info for FrictionalViscoplasticTemperatureDependent3DLaw: "
                     << rMaterialProperties[INTERNAL_FRICTION_ANGLE] << std::endl;
    }

    if (rMaterialProperties[COHESION] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing COHESION provided in process info for FrictionalViscoplasticTemperatureDependent3DLaw: "
                     << rMaterialProperties[COHESION] << std::endl;
    }

    if (rMaterialProperties[ADAPTIVE_EXPONENT] < 0.0) {
        KRATOS_ERROR << "Incorrect or missing ADAPTIVE_EXPONENT provided in process info for FrictionalViscoplasticTemperatureDependent3DLaw: "
                     << rMaterialProperties[ADAPTIVE_EXPONENT] << std::endl;
    }

    if (rMaterialProperties[BULK_MODULUS] <= 0.0) {
        KRATOS_ERROR << "Incorrect or missing BULK_MODULUS provided in process info for FrictionalViscoplasticTemperatureDependent3DLaw: "
                     << rMaterialProperties[BULK_MODULUS] << std::endl;
    }

    return 0;
}

double FrictionalViscoplasticTemperatureDependent3DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    return rParameters.GetConstitutiveMatrix()(5, 5);
}

double FrictionalViscoplasticTemperatureDependent3DLaw::GetEffectiveDensity(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    double effective_density;
    if (r_properties.HasTable(TEMPERATURE, DENSITY)) {
        effective_density = this->GetValueFromTable(TEMPERATURE, DENSITY, rParameters);
    } else {
        effective_density = r_properties[DENSITY];
    }
    return effective_density;
}

double FrictionalViscoplasticTemperatureDependent3DLaw::GetEffectiveDynamicViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    double effective_viscosity;
    if (r_properties.HasTable(TEMPERATURE, DYNAMIC_VISCOSITY)) {
        effective_viscosity = this->GetValueFromTable(TEMPERATURE, DYNAMIC_VISCOSITY, rParameters);
    } else {
        effective_viscosity = r_properties[DYNAMIC_VISCOSITY];
    }
    return effective_viscosity;
}

double FrictionalViscoplasticTemperatureDependent3DLaw::GetEffectiveFrictionAngle(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    double effective_internal_friction_angle;
    if (r_properties.HasTable(TEMPERATURE, INTERNAL_FRICTION_ANGLE)) {
        effective_internal_friction_angle = this->GetValueFromTable(TEMPERATURE, INTERNAL_FRICTION_ANGLE, rParameters);
    } else {
        effective_internal_friction_angle = r_properties[INTERNAL_FRICTION_ANGLE];
    }
    return effective_internal_friction_angle;
}

double FrictionalViscoplasticTemperatureDependent3DLaw::GetEffectiveCohesion(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    double effective_cohesion;
    if (r_properties.HasTable(TEMPERATURE, COHESION)) {
        effective_cohesion = this->GetValueFromTable(TEMPERATURE, COHESION, rParameters);
    } else {
        effective_cohesion = r_properties[COHESION];
    }
    return effective_cohesion;
}

void FrictionalViscoplasticTemperatureDependent3DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, FrictionalViscoplastic3DLaw)
}

void FrictionalViscoplasticTemperatureDependent3DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, FrictionalViscoplastic3DLaw)
}

}  // Namespace Kratos
