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
#include "custom_constitutive/solid_laws/temperature_dependent/hypoelastic_temperature_dependent_3D_law.h"
#include "includes/checks.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

    //********************************CONSTRUCTOR*********************************
    //****************************************************************************

    HypoelasticTemperatureDependent3DLaw::HypoelasticTemperatureDependent3DLaw() : Hypoelastic3DLaw() {}

    //******************************COPY CONSTRUCTOR******************************
    //****************************************************************************

    HypoelasticTemperatureDependent3DLaw::HypoelasticTemperatureDependent3DLaw(const HypoelasticTemperatureDependent3DLaw &rOther) : Hypoelastic3DLaw(rOther) {}

    //***********************************CLONE************************************
    //****************************************************************************

    ConstitutiveLaw::Pointer HypoelasticTemperatureDependent3DLaw::Clone() const { return Kratos::make_shared<HypoelasticTemperatureDependent3DLaw>(*this); }

    //*********************************DESTRUCTOR*********************************
    //****************************************************************************

    HypoelasticTemperatureDependent3DLaw::~HypoelasticTemperatureDependent3DLaw() {}

    int HypoelasticTemperatureDependent3DLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                                    const ProcessInfo &rCurrentProcessInfo) const
    {

        KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] < 0.0)
            << "Incorrect or missing YOUNG_MODULUS provided in process info for HypoelasticTemperatureDependent3DLaw: "
            << rMaterialProperties[YOUNG_MODULUS] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[POISSON_RATIO] < 0.0)
            << "Incorrect or missing POISSON_RATIO provided in process info for HypoelasticTemperatureDependent3DLaw: "
            << rMaterialProperties[POISSON_RATIO] << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0)
            << "Incorrect or missing DENSITY provided in process info for HypoelasticTemperatureDependent3DLaw: "
            << rMaterialProperties[DENSITY] << std::endl;

        return 0;
    }

    std::string HypoelasticTemperatureDependent3DLaw::Info() const { return "HypoelasticTemperatureDependent3DLaw"; }

    double HypoelasticTemperatureDependent3DLaw::GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const
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

    void HypoelasticTemperatureDependent3DLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Hypoelastic3DLaw)
    }

    void HypoelasticTemperatureDependent3DLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Hypoelastic3DLaw)
    }

} // Namespace Kratos
