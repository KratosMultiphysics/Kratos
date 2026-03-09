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
#include "custom_constitutive/fluid_laws/temperature_dependent/newtonian_temperature_dependent_3D_law.h"
#include "includes/checks.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

    //********************************CONSTRUCTOR*********************************
    //****************************************************************************

    NewtonianTemperatureDependent3DLaw::NewtonianTemperatureDependent3DLaw() : Newtonian3DLaw() {}

    //******************************COPY CONSTRUCTOR******************************
    //****************************************************************************

    NewtonianTemperatureDependent3DLaw::NewtonianTemperatureDependent3DLaw(const NewtonianTemperatureDependent3DLaw &rOther) : Newtonian3DLaw(rOther) {}

    //***********************************CLONE************************************
    //****************************************************************************

    ConstitutiveLaw::Pointer NewtonianTemperatureDependent3DLaw::Clone() const { return Kratos::make_shared<NewtonianTemperatureDependent3DLaw>(*this); }

    //*********************************DESTRUCTOR*********************************
    //****************************************************************************

    NewtonianTemperatureDependent3DLaw::~NewtonianTemperatureDependent3DLaw() {}

    int NewtonianTemperatureDependent3DLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                                  const ProcessInfo &rCurrentProcessInfo) const
    {

        KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] < 0.0)
            << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for NewtonianTemperatureDependent3DLaw: "
            << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[BULK_MODULUS] < 0.0)
            << "Incorrect or missing BULK_MODULUS provided in process info for Newtonian3DLaw: "
            << rMaterialProperties[BULK_MODULUS] << std::endl;

        return 0;
    }

    std::string NewtonianTemperatureDependent3DLaw::Info() const { return "NewtonianTemperatureDependent3DLaw"; }

    double NewtonianTemperatureDependent3DLaw::GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const
    {
        const Properties &r_properties = rParameters.GetMaterialProperties();
        double effective_variable;
        if (r_properties.HasTable(TEMPERATURE, rVariable))
        {
            effective_variable = this->GetValueFromTable(TEMPERATURE, rVariable, rParameters);
        }
        else
        {
            effective_variable = r_properties[rVariable];
        }
        return effective_variable;
    }
    void NewtonianTemperatureDependent3DLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Newtonian3DLaw)
    }

    void NewtonianTemperatureDependent3DLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Newtonian3DLaw)
    }

} // Namespace Kratos
