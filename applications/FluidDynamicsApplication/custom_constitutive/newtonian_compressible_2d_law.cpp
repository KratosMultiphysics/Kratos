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
#include "custom_constitutive/newtonian_compressible_2d_law.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

CompressibleNewtonian2DLaw::CompressibleNewtonian2DLaw()
    : FluidConstitutiveLaw()
{}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

CompressibleNewtonian2DLaw::CompressibleNewtonian2DLaw(const CompressibleNewtonian2DLaw& rOther)
    : FluidConstitutiveLaw(rOther)
{}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer CompressibleNewtonian2DLaw::Clone() const {
    return Kratos::make_shared<CompressibleNewtonian2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

CompressibleNewtonian2DLaw::~CompressibleNewtonian2DLaw() {}

ConstitutiveLaw::SizeType CompressibleNewtonian2DLaw::WorkingSpaceDimension() {
    return 2;
}

ConstitutiveLaw::SizeType CompressibleNewtonian2DLaw::GetStrainSize() {
    return 3;
}

void  CompressibleNewtonian2DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const Flags& options = rValues.GetOptions();
    const Vector& r_strain_rate = rValues.GetStrainVector();
    Vector& r_viscous_stress = rValues.GetStressVector();
    double mu_star = 0.0;
    double beta_star = 0.0;

    const double mu = this->GetEffectiveViscosity(rValues);
    GetArtificialViscosities(beta_star, mu_star, rValues);
    
    const double trace = r_strain_rate[0] + r_strain_rate[1];
    const double volumetric_part = trace/2.0; // Note: this should be small for an incompressible fluid (it is basically the incompressibility error)

    //computation of stress
    r_viscous_stress[0] = 2.0*(mu + mu_star)*(r_strain_rate[0] - volumetric_part) + beta_star*trace;
    r_viscous_stress[1] = 2.0*(mu + mu_star)*(r_strain_rate[1] - volumetric_part) + beta_star*trace;
    r_viscous_stress[2] = (mu + mu_star)*r_strain_rate[2];

    if( options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        this->NewtonianConstitutiveMatrix2D(mu + mu_star,beta_star,rValues.GetConstitutiveMatrix());
    }
}

int CompressibleNewtonian2DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Check viscosity value
    KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.0)
        << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for Newtonian2DLaw: " << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

    return 0;
}

std::string CompressibleNewtonian2DLaw::Info() const {
    return "CompressibleNewtonian2DLaw";
}

double CompressibleNewtonian2DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const
{
    const Properties &r_prop = rParameters.GetMaterialProperties();
    const double effective_viscosity = r_prop[DYNAMIC_VISCOSITY];
    return effective_viscosity;
}

void CompressibleNewtonian2DLaw::GetArtificialViscosities(double& rbeta_star, double& rmu_star, ConstitutiveLaw::Parameters& rParameters) 
{
    const SizeType n_nodes = 3;
    const GeometryType& r_geom = rParameters.GetElementGeometry();
    const array_1d<double,n_nodes>& rN = rParameters.GetShapeFunctionsValues();

    // Compute Gauss pt. interpolation value
    for (unsigned int i = 0; i < n_nodes; ++i){
        rbeta_star += rN[i] * r_geom[i].GetValue(ARTIFICIAL_BULK_VISCOSITY);
        rmu_star   += rN[i] * r_geom[i].GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
    }
}

void CompressibleNewtonian2DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

void CompressibleNewtonian2DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

} // Namespace Kratos
