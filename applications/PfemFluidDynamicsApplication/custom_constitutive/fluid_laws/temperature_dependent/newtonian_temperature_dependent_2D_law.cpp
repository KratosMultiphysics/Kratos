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
#include "custom_constitutive/fluid_laws/temperature_dependent/newtonian_temperature_dependent_2D_law.h"
#include "includes/checks.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

    //********************************CONSTRUCTOR*********************************
    //****************************************************************************

    NewtonianTemperatureDependent2DLaw::NewtonianTemperatureDependent2DLaw() : Newtonian2DLaw() {}

    //******************************COPY CONSTRUCTOR******************************
    //****************************************************************************

    NewtonianTemperatureDependent2DLaw::NewtonianTemperatureDependent2DLaw(const NewtonianTemperatureDependent2DLaw &rOther) : Newtonian2DLaw(rOther) {}

    //***********************************CLONE************************************
    //****************************************************************************

    ConstitutiveLaw::Pointer NewtonianTemperatureDependent2DLaw::Clone() const { return Kratos::make_shared<NewtonianTemperatureDependent2DLaw>(*this); }

    //*********************************DESTRUCTOR*********************************
    //****************************************************************************

    NewtonianTemperatureDependent2DLaw::~NewtonianTemperatureDependent2DLaw() {}

    int NewtonianTemperatureDependent2DLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                                  const ProcessInfo &rCurrentProcessInfo) const
    {

        KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] < 0.0)
            << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for NewtonianTemperatureDependent2DLaw: "
            << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[BULK_MODULUS] < 0.0)
            << "Incorrect or missing BULK_MODULUS provided in process info for Newtonian3DLaw: "
            << rMaterialProperties[BULK_MODULUS] << std::endl;

        return 0;
    }

    std::string NewtonianTemperatureDependent2DLaw::Info() const { return "NewtonianTemperatureDependent2DLaw"; }

    double NewtonianTemperatureDependent2DLaw::GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const
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

    void NewtonianTemperatureDependent2DLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Newtonian2DLaw)
    }

    void NewtonianTemperatureDependent2DLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Newtonian2DLaw)
    }

} // Namespace Kratos
