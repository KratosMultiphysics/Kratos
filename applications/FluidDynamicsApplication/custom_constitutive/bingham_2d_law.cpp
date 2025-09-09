//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Uxue Chasco
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

void  Bingham2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
   //-----------------------------//

    //Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();

    const Vector& S                  = rValues.GetStrainVector(); //using the short name S to reduce the length of the expressions
    Vector& StressVector                  = rValues.GetStressVector();

    //-----------------------------//

    //1.- Lame constants
    const double mu         = MaterialProperties[DYNAMIC_VISCOSITY];
    const double sigma_y    = MaterialProperties[YIELD_STRESS];
    const double m          = MaterialProperties[REGULARIZATION_COEFFICIENT];

    const double gamma_dot = std::sqrt(2.*S[0]*S[0] + 2.*S[1]*S[1] + S[2]*S[2]);

    const double min_gamma_dot = 1e-12;

    //limit the gamma_dot to a minimum so to ensure that the case of gamma_dot=0 is not problematic
    const double g = std::max(gamma_dot, min_gamma_dot);

    double Regularization = 1.0 - std::exp(-m*g);
    const double mu_effective = mu + Regularization * sigma_y / g;
    const double trS = S[0]+S[1]+S[2];
    const double eps_vol = trS/3.0;

    //computation of stress
    StressVector[0] = 2.0*mu_effective*(S[0] - eps_vol);
    StressVector[1] = 2.0*mu_effective*(S[1] - eps_vol);
    StressVector[2] = mu_effective*S[2];

    if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {
        this->NewtonianConstitutiveMatrix2D(mu_effective, rValues.GetConstitutiveMatrix());
    }
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

} // Namespace Kratos
