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
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "custom_constitutive/newtonian_2d_law.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

Newtonian2DLaw::Newtonian2DLaw()
    : FluidConstitutiveLaw()
{}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

Newtonian2DLaw::Newtonian2DLaw(const Newtonian2DLaw& rOther)
    : FluidConstitutiveLaw(rOther)
{}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer Newtonian2DLaw::Clone() const {
    return Kratos::make_shared<Newtonian2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

Newtonian2DLaw::~Newtonian2DLaw() {}

ConstitutiveLaw::SizeType Newtonian2DLaw::WorkingSpaceDimension() {
    return 2;
}

ConstitutiveLaw::SizeType Newtonian2DLaw::GetStrainSize() {
    return 3;
}

void  Newtonian2DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const Flags& options = rValues.GetOptions();
    const Vector& r_strain_rate = rValues.GetStrainVector();
    Vector& r_viscous_stress = rValues.GetStressVector();

    const double mu = this->GetEffectiveViscosity(rValues);

    const double trace = r_strain_rate[0] + r_strain_rate[1];
    const double volumetric_part = trace/3.0; // Note: this should be small for an incompressible fluid (it is basically the incompressibility error)

    //computation of stress
    r_viscous_stress[0] = 2.0*mu*(r_strain_rate[0] - volumetric_part);
    r_viscous_stress[1] = 2.0*mu*(r_strain_rate[1] - volumetric_part);
    r_viscous_stress[2] = mu*r_strain_rate[2];

    if( options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        this->NewtonianConstitutiveMatrix2D(mu,rValues.GetConstitutiveMatrix());
    }
}

int Newtonian2DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Check viscosity value
    KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.0)
        << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for Newtonian2DLaw: " << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

    return 0;
}

std::string Newtonian2DLaw::Info() const {
    return "Newtonian2DLaw";
}

double Newtonian2DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const
{
    const Properties &r_prop = rParameters.GetMaterialProperties();
    const double effective_viscosity = r_prop[DYNAMIC_VISCOSITY];
    return effective_viscosity;
}

void Newtonian2DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

void Newtonian2DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

} // Namespace Kratos
