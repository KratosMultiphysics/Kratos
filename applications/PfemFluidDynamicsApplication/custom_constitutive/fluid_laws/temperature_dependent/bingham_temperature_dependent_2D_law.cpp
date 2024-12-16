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
#include "custom_constitutive/fluid_laws/temperature_dependent/bingham_temperature_dependent_2D_law.h"
#include "includes/checks.h"
#include "includes/properties.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

    //********************************CONSTRUCTOR*********************************
    //****************************************************************************

    BinghamTemperatureDependent2DLaw::BinghamTemperatureDependent2DLaw() : Bingham2DLaw() {}

    //******************************COPY CONSTRUCTOR******************************
    //****************************************************************************

    BinghamTemperatureDependent2DLaw::BinghamTemperatureDependent2DLaw(const BinghamTemperatureDependent2DLaw &rOther) : Bingham2DLaw(rOther) {}

    //***********************************CLONE************************************
    //****************************************************************************

    ConstitutiveLaw::Pointer BinghamTemperatureDependent2DLaw::Clone() const { return Kratos::make_shared<BinghamTemperatureDependent2DLaw>(*this); }

    //*********************************DESTRUCTOR*********************************
    //****************************************************************************

    BinghamTemperatureDependent2DLaw::~BinghamTemperatureDependent2DLaw() {}

    std::string BinghamTemperatureDependent2DLaw::Info() const { return "BinghamTemperatureDependent2DLaw"; }

    //******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW******************
    //*****************************************************************************

    int BinghamTemperatureDependent2DLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                                const ProcessInfo &rCurrentProcessInfo) const
    {

        KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] < 0.0)
            << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for BinghamTemperatureDependent2DLaw: "
            << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[YIELD_SHEAR] < 0.0)
            << "Incorrect or missing YIELD_SHEAR provided in process info for BinghamTemperatureDependent2DLaw: "
            << rMaterialProperties[YIELD_SHEAR] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[ADAPTIVE_EXPONENT] < 0.0)
            << "Incorrect or missing ADAPTIVE_EXPONENT provided in process info for BinghamTemperatureDependent2DLaw: "
            << rMaterialProperties[ADAPTIVE_EXPONENT] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[BULK_MODULUS] < 0.0)
            << "Incorrect or missing BULK_MODULUS provided in process info for BinghamTemperatureDependent2DLaw: "
            << rMaterialProperties[BULK_MODULUS] << std::endl;

        return 0;
    }

    double BinghamTemperatureDependent2DLaw::GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const
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

    void BinghamTemperatureDependent2DLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Bingham2DLaw)
    }

    void BinghamTemperatureDependent2DLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Bingham2DLaw)
    }

} // Namespace Kratos
