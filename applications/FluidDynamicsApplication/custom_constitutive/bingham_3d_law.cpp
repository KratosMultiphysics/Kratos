//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/bingham_3d_law.h"

#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

Bingham3DLaw::Bingham3DLaw()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

Bingham3DLaw::Bingham3DLaw(const Bingham3DLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer Bingham3DLaw::Clone() const
{
    Bingham3DLaw::Pointer p_clone(new Bingham3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

Bingham3DLaw::~Bingham3DLaw()
{
}


void  Bingham3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
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

    //1.- Lame constants
    const double mu          = MaterialProperties[DYNAMIC_VISCOSITY];
    const double sigma_y    = MaterialProperties[YIELD_STRESS];
    const double m    = MaterialProperties[REGULARIZATION_COEFFICIENT];
    
    const double gamma_dot = std::sqrt(2.*S[0]*S[0] + 2.*S[1]*S[1] + 2.*S[2]*S[2] 
                                + S[3]*S[3] + S[4]*S[4] + S[5]*S[5]);
    
    const double min_gamma_dot = 1e-12;

    //limit the gamma_dot to a minimum so to ensure that the case of gamma_dot=0 is not problematic
    const double g = std::max(gamma_dot, min_gamma_dot);
    
    double Regularization = 1.0 - std::exp(-m*g);
    const double mu_effective = mu + Regularization * sigma_y / g;
    
//     KRATOS_WATCH( mu_effective )
    
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
        Matrix& C                  = rValues.GetConstitutiveMatrix();
        
//         if(gamma_dot < min_gamma_dot)
//         {
            noalias(C)  = ZeroMatrix(6,6);
            C(0,0) = 4.0/3.0*mu_effective;  C(0,1) = -2.0/3.0*mu_effective; C(0,2) = -2.0/3.0*mu_effective;
            C(1,2) = -2.0/3.0*mu_effective; C(1,1) = 4.0/3.0*mu_effective;  C(1,2) = -2.0/3.0*mu_effective;
            C(2,2) = -2.0/3.0*mu_effective; C(2,2) = -2.0/3.0*mu_effective; C(2,2) = 4.0/3.0*mu_effective;
            C(3,3) = mu_effective;
            C(4,4) = mu_effective;
            C(5,5) = mu_effective;
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


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void Bingham3DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = 6;

	//Set the spacedimension
	rFeatures.mSpaceDimension = 3;

}

//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
//************************************************************************************

// bool Bingham3DLaw::CheckParameters(Parameters& rValues)
// {
//     return rValues.CheckAllParameters();
// }



int Bingham3DLaw::Check(const Properties& rMaterialProperties,
                              const GeometryType& rElementGeometry,
                              const ProcessInfo& rCurrentProcessInfo)
{

    if(DYNAMIC_VISCOSITY.Key() == 0 || rMaterialProperties[DYNAMIC_VISCOSITY]<= 0.00)
        KRATOS_THROW_ERROR( std::invalid_argument,"DYNAMIC_VISCOSITY has Key zero or invalid value ", "" )

    if(YIELD_STRESS.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"YIELD_STRESS has Key zero invalid value ", "" )

    if(REGULARIZATION_COEFFICIENT.Key() == 0 || rMaterialProperties[REGULARIZATION_COEFFICIENT]<=0.00)
        KRATOS_THROW_ERROR( std::invalid_argument,"REGULARIZATION_COEFFICIENT has Key zero or invalid value ", "" )


    return 0;

}


} // Namespace Kratos
