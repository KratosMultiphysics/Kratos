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

namespace Kratos
{

    //********************************CONSTRUCTOR*********************************
    //****************************************************************************

    FrictionalViscoplasticTemperatureDependent3DLaw::FrictionalViscoplasticTemperatureDependent3DLaw() : FrictionalViscoplastic3DLaw() {}

    //******************************COPY CONSTRUCTOR******************************
    //****************************************************************************

    FrictionalViscoplasticTemperatureDependent3DLaw::FrictionalViscoplasticTemperatureDependent3DLaw(const FrictionalViscoplasticTemperatureDependent3DLaw &rOther) : FrictionalViscoplastic3DLaw(rOther) {}

    //***********************************CLONE************************************
    //****************************************************************************

    ConstitutiveLaw::Pointer FrictionalViscoplasticTemperatureDependent3DLaw::Clone() const { return Kratos::make_shared<FrictionalViscoplasticTemperatureDependent3DLaw>(*this); }

    //*********************************DESTRUCTOR*********************************
    //****************************************************************************

    FrictionalViscoplasticTemperatureDependent3DLaw::~FrictionalViscoplasticTemperatureDependent3DLaw() {}

    std::string FrictionalViscoplasticTemperatureDependent3DLaw::Info() const { return "FrictionalViscoplasticTemperatureDependent3DLaw"; }

    //******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW******************
    //*****************************************************************************

    int FrictionalViscoplasticTemperatureDependent3DLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                                               const ProcessInfo &rCurrentProcessInfo) const
    {

        KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] < 0.0)
            << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for FrictionalViscoplasticTemperatureDependent3DLaw: "
            << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[INTERNAL_FRICTION_ANGLE] < 0.0)
            << "Incorrect or missing INTERNAL_FRICTION_ANGLE provided in process info for FrictionalViscoplasticTemperatureDependent3DLaw: "
            << rMaterialProperties[INTERNAL_FRICTION_ANGLE] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[COHESION] < 0.0)
            << "Incorrect or missing COHESION provided in process info for FrictionalViscoplasticTemperatureDependent3DLaw: "
            << rMaterialProperties[COHESION] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[ADAPTIVE_EXPONENT] < 0.0)
            << "Incorrect or missing ADAPTIVE_EXPONENT provided in process info for FrictionalViscoplasticTemperatureDependent3DLaw: "
            << rMaterialProperties[ADAPTIVE_EXPONENT] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[BULK_MODULUS] < 0.0)
            << "Incorrect or missing BULK_MODULUS provided in process info for FrictionalViscoplasticTemperatureDependent3DLaw: "
            << rMaterialProperties[BULK_MODULUS] << std::endl;

        return 0;
    }

    double FrictionalViscoplasticTemperatureDependent3DLaw::GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const
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

    void FrictionalViscoplasticTemperatureDependent3DLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, FrictionalViscoplastic3DLaw)
    }

    void FrictionalViscoplasticTemperatureDependent3DLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, FrictionalViscoplastic3DLaw)
    }

} // Namespace Kratos
