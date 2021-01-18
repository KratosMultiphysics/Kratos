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
#include "custom_constitutive/newtonian_compressible_3d_law.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

CompressibleNewtonian3DLaw::CompressibleNewtonian3DLaw()
    : FluidConstitutiveLaw()
{}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

CompressibleNewtonian3DLaw::CompressibleNewtonian3DLaw(const CompressibleNewtonian3DLaw& rOther)
    : FluidConstitutiveLaw(rOther)
{}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer CompressibleNewtonian3DLaw::Clone() const {
    return Kratos::make_shared<CompressibleNewtonian3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

CompressibleNewtonian3DLaw::~CompressibleNewtonian3DLaw() {}

ConstitutiveLaw::SizeType CompressibleNewtonian3DLaw::WorkingSpaceDimension() {
    return 3;
}

ConstitutiveLaw::SizeType CompressibleNewtonian3DLaw::GetStrainSize() {
    return 6;
}

void  CompressibleNewtonian3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    const Flags& options = rValues.GetOptions();
    const Vector& r_strain_rate = rValues.GetStrainVector();
    Vector& r_viscous_stress = rValues.GetStressVector();
    double mu_star = 0.0;
    double beta_star = 0.0;
    
    const double mu = this->GetEffectiveViscosity(rValues);
    GetArtificialViscosities(beta_star, mu_star, rValues);

    const double trace = r_strain_rate[0] + r_strain_rate[1] + r_strain_rate[2];
    const double volumetric_part = trace/3.0; // Note: this should be small for an incompressible fluid (it is basically the incompressibility error)

    //computation of stress
    r_viscous_stress[0] = 2.0*(mu + mu_star)*(r_strain_rate[0] - volumetric_part) + beta_star*trace;
    r_viscous_stress[1] = 2.0*(mu + mu_star)*(r_strain_rate[1] - volumetric_part) + beta_star*trace;
    r_viscous_stress[2] = 2.0*(mu + mu_star)*(r_strain_rate[2] - volumetric_part) + beta_star*trace;
    r_viscous_stress[3] = (mu + mu_star)*r_strain_rate[3];
    r_viscous_stress[4] = (mu + mu_star)*r_strain_rate[4];
    r_viscous_stress[5] = (mu + mu_star)*r_strain_rate[5];

    if( options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        this->NewtonianConstitutiveMatrix3D(mu + mu_star,beta_star,rValues.GetConstitutiveMatrix());
        //basetype::    
    }
}

int CompressibleNewtonian3DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{

    // Check viscosity value
    KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.0)
        << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for Newtonian2DLaw: " << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

    for (unsigned int i = 0; i < rElementGeometry.size(); i++) {
        const Node<3>& rNode = rElementGeometry[i];
        KRATOS_ERROR_IF(rNode.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) < 0.0)
            << "ARTIFICIAL_DYNAMIC_VISCOSITY was not correctly assigned to nodes for Constitutive Law.\n";
        KRATOS_ERROR_IF(rNode.GetValue(ARTIFICIAL_BULK_VISCOSITY) < 0.0)
            << "ARTIFICIAL_BULK_VISCOSITY was not correctly assigned to nodes for Constitutive Law.\n";
    }
    return 0;
}

std::string CompressibleNewtonian3DLaw::Info() const {
    return "CompressibleNewtonian3DLaw";
}

double CompressibleNewtonian3DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const
{
    const Properties &r_prop = rParameters.GetMaterialProperties();
    const double effective_viscosity = r_prop[DYNAMIC_VISCOSITY];
    return effective_viscosity;
}

void CompressibleNewtonian3DLaw::GetArtificialViscosities(double& rbeta_star, double& rmu_star, ConstitutiveLaw::Parameters& rParameters) 
{
    const GeometryType& r_geom = rParameters.GetElementGeometry();
    const SizeType n_nodes = r_geom.size();
    const auto& rN = rParameters.GetShapeFunctionsValues();
    
    // Compute Gauss pt. interpolation value
    for (unsigned int i = 0; i < n_nodes; ++i){
        
        rbeta_star += rN[i] * r_geom[i].GetValue(ARTIFICIAL_BULK_VISCOSITY);
        rmu_star   += rN[i] * r_geom[i].GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
    }
}

void CompressibleNewtonian3DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

void CompressibleNewtonian3DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

} // Namespace Kratos
