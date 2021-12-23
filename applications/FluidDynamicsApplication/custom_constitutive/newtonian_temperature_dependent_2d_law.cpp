//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/cfd_variables.h"

// Application includes
#include "custom_constitutive/newtonian_temperature_dependent_2d_law.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NewtonianTemperatureDependent2DLaw::NewtonianTemperatureDependent2DLaw()
    : Newtonian2DLaw()
{}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NewtonianTemperatureDependent2DLaw::NewtonianTemperatureDependent2DLaw(const NewtonianTemperatureDependent2DLaw& rOther)
    : Newtonian2DLaw(rOther)
{}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NewtonianTemperatureDependent2DLaw::Clone() const
{
    return Kratos::make_shared<NewtonianTemperatureDependent2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NewtonianTemperatureDependent2DLaw::~NewtonianTemperatureDependent2DLaw()
{}

int NewtonianTemperatureDependent2DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    // If the viscosity is not table provided, check its value
    if (!rMaterialProperties.HasTable(TEMPERATURE, DYNAMIC_VISCOSITY)) {
        KRATOS_ERROR << "TEMPERATURE - DYNAMICS_VISCOSITY table is missing in NewtonianTemperatureDependent2DLaw" << std::endl;
    }

    return 0;
}

std::string NewtonianTemperatureDependent2DLaw::Info() const
{
    return "NewtonianTemperatureDependent2DLaw";
}

double NewtonianTemperatureDependent2DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters &rParameters) const
{
    const double effective_viscosity = this->GetValueFromTable(TEMPERATURE, DYNAMIC_VISCOSITY, rParameters);
    return effective_viscosity;
}

void NewtonianTemperatureDependent2DLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Newtonian2DLaw )
}

void NewtonianTemperatureDependent2DLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Newtonian2DLaw )
}

} // Namespace Kratos
