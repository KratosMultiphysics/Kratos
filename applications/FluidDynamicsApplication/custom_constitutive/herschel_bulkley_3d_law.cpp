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
#include "custom_constitutive/herschel_bulkley_3d_law.h"

#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HerschelBulkley3DLaw::HerschelBulkley3DLaw()
    : FluidConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HerschelBulkley3DLaw::HerschelBulkley3DLaw(const HerschelBulkley3DLaw& rOther)
    : FluidConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HerschelBulkley3DLaw::Clone() const
{
    return Kratos::make_shared<HerschelBulkley3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HerschelBulkley3DLaw::~HerschelBulkley3DLaw()
{
}

ConstitutiveLaw::SizeType HerschelBulkley3DLaw::WorkingSpaceDimension() {
    return 3;
}

ConstitutiveLaw::SizeType HerschelBulkley3DLaw::GetStrainSize() {
    return 6;
}

void  HerschelBulkley3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    //CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();
    
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();    

    Vector& S                  = rValues.GetStrainVector(); //using the short name S to reduce the lenght of the expressions
    Vector& StressVector                  = rValues.GetStressVector();

    //-----------------------------//
    const double sigma_y    = MaterialProperties[YIELD_STRESS];
    const double m    = MaterialProperties[REGULARIZATION_COEFFICIENT];
    const double k = MaterialProperties[POWER_LAW_K];
    const double n = MaterialProperties[POWER_LAW_N];
        
    const double gamma_dot = std::sqrt(2.*S[0]*S[0] + 2.*S[1]*S[1] + 2.*S[2]*S[2] 
                                + S[3]*S[3] + S[4]*S[4] + S[5]*S[5]);
    
    const double min_gamma_dot = 1e-6; 

    //limit the gamma_dot to a minimum so to ensure that the case of gamma_dot=0 is not problematic
    const double g = std::max(gamma_dot, min_gamma_dot);
    
    const double Regularization = 1.0 - std::exp(-m*g);
    const double mu_effective = k*std::pow(g,n-1) + Regularization * sigma_y / g;
    
    const double trS = S[0]+S[1]+S[2];
    const double eps_vol = trS/3.0;
    
    //computation of stress
    StressVector[0] = 2.0*mu_effective*(S[0] - eps_vol);
    StressVector[1] = 2.0*mu_effective*(S[1] - eps_vol);
    StressVector[2] = 2.0*mu_effective*(S[2] - eps_vol);
    StressVector[3] = mu_effective*S[3];
    StressVector[4] = mu_effective*S[4];
    StressVector[5] = mu_effective*S[5];

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
//         if(gamma_dot < min_gamma_dot)
//         {
        this->NewtonianConstitutiveMatrix3D(mu_effective,rValues.GetConstitutiveMatrix());
//         }
//         else
//         {
// 
//                 const double x0 =             pow(S[3], 2);
//                 const double x1 =             pow(S[4], 2);
//                 const double x2 =             pow(S[5], 2);
//                 const double x3 =             2*pow(S[0], 2) + 2*pow(S[1], 2) + 2*pow(S[2], 2) + x0 + x1 + x2;
//                 const double x4 =             pow(x3, 3);
//                 const double x5 =             1.0/x4;
//                 const double x6 =             sqrt(x3);
//                 const double x7 =             m*x6;
//                 const double x8 =             exp(-x7);
//                 const double x9 =             (4.0L/3.0L)*x5*x8;
//                 const double x10 =             exp(x7);
//                 const double x11 =             mu*x10*x4;
//                 const double x12 =             pow(x3, 5.0L/2.0L);
//                 const double x13 =             x10 - 1;
//                 const double x14 =             sigma_y*x12*x13;
//                 const double x15 =             x11 + x14;
//                 const double x16 =             2*S[0];
//                 const double x17 =             S[1] + S[2] - x16;
//                 const double x18 =             pow(x3, 3.0L/2.0L);
//                 const double x19 =             m*x18 - x13*x3;
//                 const double x20 =             (1.0L/3.0L)*x5*x8;
//                 const double x21 =             2*x11 + 2*x14;
//                 const double x22 =             4*S[1]*sigma_y*x19*x6;
//                 const double x23 =             4*S[2]*sigma_y*x19*x6;
//                 const double x24 =             1.0/x12;
//                 const double x25 =             (2.0L/3.0L)*S[3]*sigma_y*x19*x24*x8;
//                 const double x26 =             (2.0L/3.0L)*S[4]*sigma_y*x19*x24*x8;
//                 const double x27 =             (2.0L/3.0L)*S[5]*sigma_y*x19*x24*x8;
//                 const double x28 =             2*S[1];
//                 const double x29 =             S[0] + S[2] - x28;
//                 const double x30 =             4*S[0]*sigma_y*x19*x6;
//                 const double x31 =             2*S[2];
//                 const double x32 =             S[0] + S[1] - x31;
//                 const double x33 =             1.0/x18;
//                 const double x34 =             -x10 + x7 + 1;
//                 const double x35 =             S[3]*sigma_y*x33*x34*x8;
//                 const double x36 =             x5*x8;
//                 const double x37 =             sigma_y*x19*x6;
//                 const double x38 =             S[4]*x35;
//                 const double x39 =             S[5]*x35;
//                 const double x40 =             S[4]*sigma_y*x33*x34*x8;
//                 const double x41 =             S[5]*x40;
//                 const double x42 =             S[5]*sigma_y*x33*x34*x8;
//                 C(0,0)=x9*(-S[0]*sigma_y*x17*x19*x6 + x15);
//                 C(0,1)=-x20*(x17*x22 + x21);
//                 C(0,2)=-x20*(x17*x23 + x21);
//                 C(0,3)=-x17*x25;
//                 C(0,4)=-x17*x26;
//                 C(0,5)=-x17*x27;
//                 C(1,0)=-x20*(x21 + x29*x30);
//                 C(1,1)=x9*(-S[1]*sigma_y*x19*x29*x6 + x15);
//                 C(1,2)=-x20*(x21 + x23*x29);
//                 C(1,3)=-x25*x29;
//                 C(1,4)=-x26*x29;
//                 C(1,5)=-x27*x29;
//                 C(2,0)=-x20*(x21 + x30*x32);
//                 C(2,1)=-x20*(x21 + x22*x32);
//                 C(2,2)=x9*(-S[2]*sigma_y*x19*x32*x6 + x15);
//                 C(2,3)=-x25*x32;
//                 C(2,4)=-x26*x32;
//                 C(2,5)=-x27*x32;
//                 C(3,0)=x16*x35;
//                 C(3,1)=x28*x35;
//                 C(3,2)=x31*x35;
//                 C(3,3)=x36*(x0*x37 + x15);
//                 C(3,4)=x38;
//                 C(3,5)=x39;
//                 C(4,0)=x16*x40;
//                 C(4,1)=x28*x40;
//                 C(4,2)=x31*x40;
//                 C(4,3)=x38;
//                 C(4,4)=x36*(x1*x37 + x15);
//                 C(4,5)=x41;
//                 C(5,0)=x16*x42;
//                 C(5,1)=x28*x42;
//                 C(5,2)=x31*x42;
//                 C(5,3)=x39;
//                 C(5,4)=x41;
//                 C(5,5)=x36*(x15 + x2*x37);
// 
//         }
            
    }
	  
}


//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
//************************************************************************************

int HerschelBulkley3DLaw::Check(const Properties& rMaterialProperties,
                              const GeometryType& rElementGeometry,
                              const ProcessInfo& rCurrentProcessInfo)
{    
    KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS);
    KRATOS_CHECK_VARIABLE_KEY(REGULARIZATION_COEFFICIENT);
    KRATOS_CHECK_VARIABLE_KEY(POWER_LAW_K);
    KRATOS_CHECK_VARIABLE_KEY(POWER_LAW_N);

    if( rMaterialProperties[YIELD_STRESS] <= 0.00 ) {
        KRATOS_ERROR << "Incorrect or missing YIELD_STRESS provided in process info for HerschelBulkley3DLaw: " << rMaterialProperties[YIELD_STRESS] << std::endl;
    }

    if( rMaterialProperties[REGULARIZATION_COEFFICIENT] <= 0.00 ) {
        KRATOS_ERROR << "Incorrect or missing REGULARIZATION_COEFFICIENT provided in process info for HerschelBulkley3DLaw: " << rMaterialProperties[REGULARIZATION_COEFFICIENT] << std::endl;
    }

    if( rMaterialProperties[POWER_LAW_K] <= 0.00 ) {
        KRATOS_ERROR << "Incorrect or missing POWER_LAW_K provided in process info for HerschelBulkley3DLaw: " << rMaterialProperties[POWER_LAW_K] << std::endl;
    }

    if( rMaterialProperties[POWER_LAW_N] <= 0.00 ) {
        KRATOS_ERROR << "Incorrect or missing POWER_LAW_N provided in process info for HerschelBulkley3DLaw: " << rMaterialProperties[POWER_LAW_N] << std::endl;
    }

    return 0;

}

std::string HerschelBulkley3DLaw::Info() const {
    return "HerschelBulkley3DLaw";
}

double HerschelBulkley3DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    // We are abusing the fact that C(5,5) = mu_effective
    return rParameters.GetConstitutiveMatrix()(5,5);
}

void HerschelBulkley3DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

void HerschelBulkley3DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

} // Namespace Kratos
