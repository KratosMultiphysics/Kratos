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
#include "custom_constitutive/newtonian_temperature_dependent_3d_law.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NewtonianTemperatureDependent3DLaw::NewtonianTemperatureDependent3DLaw()
    : Newtonian3DLaw()
{}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NewtonianTemperatureDependent3DLaw::NewtonianTemperatureDependent3DLaw(const NewtonianTemperatureDependent3DLaw& rOther)
    : Newtonian3DLaw(rOther)
{}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NewtonianTemperatureDependent3DLaw::Clone() const
{
    return Kratos::make_shared<NewtonianTemperatureDependent3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NewtonianTemperatureDependent3DLaw::~NewtonianTemperatureDependent3DLaw()
{}

int NewtonianTemperatureDependent3DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    // If the viscosity is not table provided, check its value
    if (!rMaterialProperties.HasTable(TEMPERATURE, DYNAMIC_VISCOSITY)) {
        KRATOS_ERROR << "TEMPERATURE - DYNAMICS_VISCOSITY table is missing in NewtonianTemperatureDependent3DLaw" << std::endl;
    }

    return 0;
}

std::string NewtonianTemperatureDependent3DLaw::Info() const
{
    return "NewtonianTemperatureDependent3DLaw";
}

double NewtonianTemperatureDependent3DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters &rParameters) const
{
    const double effective_viscosity = this->GetValueFromTable(TEMPERATURE, DYNAMIC_VISCOSITY, rParameters);
    return effective_viscosity;
}

void NewtonianTemperatureDependent3DLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Newtonian3DLaw )
}

void NewtonianTemperatureDependent3DLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Newtonian3DLaw )
}

} // Namespace Kratos
