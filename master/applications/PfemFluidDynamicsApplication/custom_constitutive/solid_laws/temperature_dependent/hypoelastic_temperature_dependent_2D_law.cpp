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
#include "custom_constitutive/solid_laws/temperature_dependent/hypoelastic_temperature_dependent_2D_law.h"
#include "includes/checks.h"
#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

    //********************************CONSTRUCTOR*********************************
    //****************************************************************************

    HypoelasticTemperatureDependent2DLaw::HypoelasticTemperatureDependent2DLaw() : Hypoelastic2DLaw() {}

    //******************************COPY CONSTRUCTOR******************************
    //****************************************************************************

    HypoelasticTemperatureDependent2DLaw::HypoelasticTemperatureDependent2DLaw(const HypoelasticTemperatureDependent2DLaw &rOther) : Hypoelastic2DLaw(rOther) {}

    //***********************************CLONE************************************
    //****************************************************************************

    ConstitutiveLaw::Pointer HypoelasticTemperatureDependent2DLaw::Clone() const { return Kratos::make_shared<HypoelasticTemperatureDependent2DLaw>(*this); }

    //*********************************DESTRUCTOR*********************************
    //****************************************************************************

    HypoelasticTemperatureDependent2DLaw::~HypoelasticTemperatureDependent2DLaw() {}

    int HypoelasticTemperatureDependent2DLaw::Check(const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
                                                    const ProcessInfo &rCurrentProcessInfo) const
    {

        KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] < 0.0)
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

    double HypoelasticTemperatureDependent2DLaw::GetEffectiveMaterialParameter(ConstitutiveLaw::Parameters &rParameters, const Variable<double> &rVariable) const
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

    void HypoelasticTemperatureDependent2DLaw::save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Hypoelastic2DLaw)
    }

    void HypoelasticTemperatureDependent2DLaw::load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Hypoelastic2DLaw)
    }

} // Namespace Kratos
