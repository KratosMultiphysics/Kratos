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
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/bingham_2d_law.h"

#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

Bingham2DLaw::Bingham2DLaw()
    : FluidConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

Bingham2DLaw::Bingham2DLaw(const Bingham2DLaw& rOther)
    : FluidConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer Bingham2DLaw::Clone() const
{
    return Kratos::make_shared<Bingham2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

Bingham2DLaw::~Bingham2DLaw()
{
}

ConstitutiveLaw::SizeType Bingham2DLaw::WorkingSpaceDimension() {
    return 2;
}

ConstitutiveLaw::SizeType Bingham2DLaw::GetStrainSize() const {
    return 3;
}


void  Bingham2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    //CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();

    Vector& S = rValues.GetStrainVector(); //using the short name S to reduce the lenght of the expressions
    Vector& StressVector = rValues.GetStressVector();

    //-----------------------------//

    //1.- Lame constants
    const double mu          = MaterialProperties[DYNAMIC_VISCOSITY];
    // KRATOS_WATCH(mu)

    const double sigma_y    = MaterialProperties[YIELD_STRESS];
    const double m    = MaterialProperties[REGULARIZATION_COEFFICIENT];

    const double gamma_dot = std::sqrt(2.*S[0]*S[0] + 2.*S[1]*S[1] + S[2]*S[2]);

    const double min_gamma_dot = 1e-12;

    //limit the gamma_dot to a minimum so to ensure that the case of gamma_dot=0 is not problematic
    const double g = std::max(gamma_dot, min_gamma_dot);

    double Regularization = 1.0 - std::exp(-m*g);
    const double mu_effective = mu + Regularization * sigma_y / g;
    mMuEffective = mu_effective;
    
    const double trS = S[0]+S[1];
    const double eps_vol = trS/3.0;

    //computation of stress
    StressVector[0] = 2.0*mu_effective*(S[0]);
    StressVector[1] = 2.0*mu_effective*(S[1]);
    StressVector[2] = mu_effective*S[2];

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        this->NewtonianConstitutiveMatrix2D(mu_effective,rValues.GetConstitutiveMatrix());
        // BinghamConstitutiveMatrix2D(g, m, sigma_y, mu_effective, S, rValues.GetConstitutiveMatrix());

    }

}

void Bingham2DLaw::BinghamConstitutiveMatrix2D(
    const double g,
    const double m,
    const double sigma_y, 
    const double mu_effective, 
    Vector& S,
    Matrix& rC) {

    const double d_mu_d_g = -sigma_y / (g * g) * (1.0 - std::exp(-m * g)) + (m * sigma_y / g) * std::exp(-m * g);
    
    rC(0,0) = 4.0 / 3.0 * mu_effective + 4.0 * (2.0 / 3.0 * S[0] - 1.0 / 3.0 * S[1]) * S[0] / g * d_mu_d_g;
    rC(0,1) = -2.0 / 3.0 * mu_effective + 4.0 * (2.0 / 3.0 * S[0] - 1.0 / 3.0 * S[1]) * S[1] / g * d_mu_d_g;
    rC(0,2) = 2.0 * (2.0 / 3.0 * S[0] - 1.0 / 3.0 * S[1]) * S[2] / g * d_mu_d_g;

    rC(1,0) = -2.0 / 3.0 * mu_effective + 4.0 * (2.0 / 3.0 * S[1] - 1.0 / 3.0 * S[0]) * S[0] / g * d_mu_d_g;
    rC(1,1) = 4.0 / 3.0 * mu_effective + 4.0 * (2.0 / 3.0 * S[1] - 1.0 / 3.0 * S[0]) * S[1] / g * d_mu_d_g;
    rC(1,2) = 2.0 * (2.0 / 3.0 * S[1] - 1.0 / 3.0 * S[0]) * S[2] / g * d_mu_d_g;

    rC(2,0) = S[2] * 2.0 * S[0] / g * d_mu_d_g;
    rC(2,1) = S[2] * 2.0 * S[1] / g * d_mu_d_g;
    rC(2,2) = mu_effective + S[2] * S[2] / g * d_mu_d_g;

    
    // rC(0,0) = 2 * mu_effective + 4 * S[0] * d_mu_d_g * S[0] / g; 
    // rC(0,1) = 4 * S[0] * d_mu_d_g * S[1] / g;
    // rC(0,2) = 2 * S[0] * d_mu_d_g * S[2] / g;

    // rC(1,0) = 4 * S[1] * d_mu_d_g * S[0] / g;
    // rC(1,1) = 2 * mu_effective + 4 * S[1] * d_mu_d_g * S[1] / g;
    // rC(1,2) = 2 * S[1] * d_mu_d_g * S[2] / g;

    // rC(2,0) = 2 * S[2] * d_mu_d_g * S[0] / g;
    // rC(2,1) = 2 * S[2] * d_mu_d_g * S[1] / g;
    // rC(2,2) = mu_effective + S[2] * d_mu_d_g * S[2] / g;
}

std::string Bingham2DLaw::Info() const {
    return "Bingham2DLaw";
}


//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
//************************************************************************************

int Bingham2DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if( rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.00 ) {
        KRATOS_ERROR << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for Bingham2DLaw: " << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;
    }

    if( rMaterialProperties[YIELD_STRESS] <= 0.00 ) {
        KRATOS_ERROR << "Incorrect or missing YIELD_STRESS provided in process info for Bingham2DLaw: " << rMaterialProperties[YIELD_STRESS] << std::endl;
    }

    if( rMaterialProperties[REGULARIZATION_COEFFICIENT] <= 0.00 ) {
        KRATOS_ERROR << "Incorrect or missing REGULARIZATION_COEFFICIENT provided in process info for Bingham2DLaw: " << rMaterialProperties[REGULARIZATION_COEFFICIENT] << std::endl;
    }

    return 0;

}


double Bingham2DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    // We are abusing the fact that C(5,5) = mu_effective
    return rParameters.GetConstitutiveMatrix()(5,5);
}

void Bingham2DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

void Bingham2DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

double& Bingham2DLaw::CalculateValue(
    Parameters& rParameters,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    const Properties& MaterialProperties = rParameters.GetMaterialProperties();
    const double sigma_y = MaterialProperties[YIELD_STRESS];
    Vector& StressVector = rParameters.GetStressVector();
    const double tau = std::sqrt(StressVector[0] * StressVector[0] + 
                                 StressVector[1] * StressVector[1] + 
                                 0.5 * StressVector[2] * StressVector[2]);
    double yielded_state = (tau < sigma_y) ? 1.0 : 0.0;

    if (rThisVariable == STRAIN_ENERGY) {
        rValue = yielded_state;
    }
    else if (rThisVariable == MU) {
        rValue = mMuEffective;
    }
    return rValue;
}

} // Namespace Kratos
