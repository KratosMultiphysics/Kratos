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
#include "custom_constitutive/fluid_laws/newtonian_temperature_dependent_3D_law.h"
#include "includes/checks.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos {

//********************************CONSTRUCTOR*********************************
//****************************************************************************

NewtonianTemperatureDependent3DLaw::NewtonianTemperatureDependent3DLaw() : Newtonian3DLaw() {}

//******************************COPY CONSTRUCTOR******************************
//****************************************************************************

NewtonianTemperatureDependent3DLaw::NewtonianTemperatureDependent3DLaw(const NewtonianTemperatureDependent3DLaw& rOther) : Newtonian3DLaw(rOther) {}

//***********************************CLONE************************************
//****************************************************************************

ConstitutiveLaw::Pointer NewtonianTemperatureDependent3DLaw::Clone() const { return Kratos::make_shared<NewtonianTemperatureDependent3DLaw>(*this); }

//*********************************DESTRUCTOR*********************************
//****************************************************************************

NewtonianTemperatureDependent3DLaw::~NewtonianTemperatureDependent3DLaw() {}

int NewtonianTemperatureDependent3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                          const ProcessInfo& rCurrentProcessInfo) {

    KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.0)
        << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for NewtonianTemperatureDependent3DLaw: "
        << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[BULK_MODULUS] <= 0.0)
        << "Incorrect or missing BULK_MODULUS provided in process info for Newtonian3DLaw: "
        << rMaterialProperties[BULK_MODULUS] << std::endl;

    return 0;
}

std::string NewtonianTemperatureDependent3DLaw::Info() const { return "NewtonianTemperatureDependent3DLaw"; }

double NewtonianTemperatureDependent3DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    double effective_viscosity;
    if (r_properties.HasTable(TEMPERATURE, DYNAMIC_VISCOSITY)) {
        effective_viscosity = this->GetValueFromTable(TEMPERATURE, DYNAMIC_VISCOSITY, rParameters);
    } else {
        effective_viscosity = r_properties[DYNAMIC_VISCOSITY];
    }
    return effective_viscosity;
}

double NewtonianTemperatureDependent3DLaw::GetEffectiveDensity(ConstitutiveLaw::Parameters& rParameters) const {
    const Properties& r_properties = rParameters.GetMaterialProperties();
    double effective_density;
    if (r_properties.HasTable(TEMPERATURE, DENSITY)) {
        effective_density = this->GetValueFromTable(TEMPERATURE, DENSITY, rParameters);
    } else {
        effective_density = r_properties[DENSITY];
    }
    return effective_density;
}

void NewtonianTemperatureDependent3DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Newtonian3DLaw)
}

void NewtonianTemperatureDependent3DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Newtonian3DLaw)
}

}  // Namespace Kratos
