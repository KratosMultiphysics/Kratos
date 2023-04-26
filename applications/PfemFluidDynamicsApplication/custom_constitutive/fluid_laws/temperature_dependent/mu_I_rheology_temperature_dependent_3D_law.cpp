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
#include "custom_constitutive/fluid_laws/temperature_dependent/mu_I_rheology_temperature_dependent_3D_law.h"
#include "includes/checks.h"
#include "includes/properties.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

    //********************************CONSTRUCTOR*********************************
    //****************************************************************************

    MuIRheologyTemperatureDependent3DLaw::MuIRheologyTemperatureDependent3DLaw() : MuIRheology3DLaw() {}

    //******************************COPY CONSTRUCTOR******************************
    //****************************************************************************

    MuIRheologyTemperatureDependent3DLaw::MuIRheologyTemperatureDependent3DLaw(const MuIRheologyTemperatureDependent3DLaw &rOther) : MuIRheology3DLaw(rOther) {}

    //***********************************CLONE************************************
    //****************************************************************************

    ConstitutiveLaw::Pointer MuIRheologyTemperatureDependent3DLaw::Clone() const { return Kratos::make_shared<MuIRheologyTemperatureDependent3DLaw>(*this); }

    //*********************************DESTRUCTOR*********************************
    //****************************************************************************

    MuIRheologyTemperatureDependent3DLaw::~MuIRheologyTemperatureDependent3DLaw() {}

    std::string MuIRheologyTemperatureDependent3DLaw::Info() const { return "MuIRheologyTemperatureDependent3DLaw"; }

    //******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW******************
    //*****************************************************************************

    int MuIRheologyTemperatureDependent3DLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                                    const ProcessInfo &rCurrentProcessInfo) const
    {

        KRATOS_ERROR_IF(rMaterialProperties[STATIC_FRICTION] < 0.0)
            << "Incorrect or missing STATIC_FRICTION provided in process info for MuIRheologyTemperatureDependent3DLaw: "
            << rMaterialProperties[STATIC_FRICTION] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_FRICTION] < 0.0)
            << "Incorrect or missing DYNAMIC_FRICTION provided in process info for MuIRheologyTemperatureDependent3DLaw: "
            << rMaterialProperties[DYNAMIC_FRICTION] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[INERTIAL_NUMBER_ZERO] < 0.0)
            << "Incorrect or missing INERTIAL_NUMBER_ZERO provided in process info for MuIRheologyTemperatureDependent3DLaw: "
            << rMaterialProperties[INERTIAL_NUMBER_ZERO] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[GRAIN_DIAMETER] < 0.0)
            << "Incorrect or missing GRAIN_DIAMETER provided in process info for MuIRheologyTemperatureDependent3DLaw: "
            << rMaterialProperties[GRAIN_DIAMETER] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[GRAIN_DENSITY] < 0.0)
            << "Incorrect or missing GRAIN_DENSITY provided in process info for MuIRheologyTemperatureDependent3DLaw: "
            << rMaterialProperties[GRAIN_DENSITY] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[REGULARIZATION_COEFFICIENT] < 0.0)
            << "Incorrect or missing REGULARIZATION_COEFFICIENT provided in process info for MuIRheologyTemperatureDependent3DLaw: "
            << rMaterialProperties[REGULARIZATION_COEFFICIENT] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[BULK_MODULUS] < 0.0)
            << "Incorrect or missing BULK_MODULUS provided in process info for MuIRheologyTemperatureDependent3DLaw: "
            << rMaterialProperties[BULK_MODULUS] << std::endl;

        return 0;
    }

    double MuIRheologyTemperatureDependent3DLaw::GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const
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

    void MuIRheologyTemperatureDependent3DLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MuIRheology3DLaw)
    }

    void MuIRheologyTemperatureDependent3DLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MuIRheology3DLaw)
    }

} // Namespace Kratos
