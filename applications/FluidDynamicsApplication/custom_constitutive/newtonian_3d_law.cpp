//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "custom_constitutive/newtonian_3d_law.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

Newtonian3DLaw::Newtonian3DLaw()
    : FluidConstitutiveLaw()
{}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

Newtonian3DLaw::Newtonian3DLaw(const Newtonian3DLaw& rOther)
    : FluidConstitutiveLaw(rOther)
{}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer Newtonian3DLaw::Clone() const {
    return Kratos::make_shared<Newtonian3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

Newtonian3DLaw::~Newtonian3DLaw() {}

ConstitutiveLaw::SizeType Newtonian3DLaw::WorkingSpaceDimension() {
    return 3;
}

ConstitutiveLaw::SizeType Newtonian3DLaw::GetStrainSize() {
    return 6;
}

void  Newtonian3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues) {
    const Flags& options = rValues.GetOptions();
    
    const Vector& r_strain_rate = rValues.GetStrainVector();
    Vector& r_viscous_stress = rValues.GetStressVector();

    const double mu = ComputeEffectiveViscosity(rValues);

    const double trace = r_strain_rate[0] + r_strain_rate[1] + r_strain_rate[2];
    const double volumetric_part = trace/3.0; // Note: this should be small for an incompressible fluid (it is basically the incompressibility error)

    //computation of stress
    r_viscous_stress[0] = 2.0*mu*(r_strain_rate[0] - volumetric_part);
    r_viscous_stress[1] = 2.0*mu*(r_strain_rate[1] - volumetric_part);
    r_viscous_stress[2] = 2.0*mu*(r_strain_rate[2] - volumetric_part);
    r_viscous_stress[3] = mu*r_strain_rate[3];
    r_viscous_stress[4] = mu*r_strain_rate[4];
    r_viscous_stress[5] = mu*r_strain_rate[5];

    if( options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        this->NewtonianConstitutiveMatrix3D(mu,rValues.GetConstitutiveMatrix());
    }
}

int Newtonian3DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) {
    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_VISCOSITY);

    if( rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.00 ) {
        KRATOS_ERROR << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for Newtonian3DLaw: " << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;
    }

    return 0;
}

std::string Newtonian3DLaw::Info() const {
    return "Newtonian3DLaw";
}

double Newtonian3DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    // We are abusing the fact that C(5,5) = mu
    return rParameters.GetConstitutiveMatrix()(5,5);
}

double Newtonian3DLaw::ComputeEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    return rParameters.GetMaterialProperties()[DYNAMIC_VISCOSITY];
}

void Newtonian3DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

void Newtonian3DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

} // Namespace Kratos
