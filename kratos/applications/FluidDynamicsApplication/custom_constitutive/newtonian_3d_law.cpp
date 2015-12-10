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
#include "custom_constitutive/newtonian_3d_law.h"

#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

Newtonian3DLaw::Newtonian3DLaw()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

Newtonian3DLaw::Newtonian3DLaw(const Newtonian3DLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer Newtonian3DLaw::Clone() const
{
    Newtonian3DLaw::Pointer p_clone(new Newtonian3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

Newtonian3DLaw::~Newtonian3DLaw()
{
}


void  Newtonian3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
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
    const double trS = S[0]+S[1]+S[2];
    const double eps_vol = trS/3.0;
    
    //computation of stress
    StressVector[0] = 2.0*mu*(S[0] - eps_vol);
    StressVector[1] = 2.0*mu*(S[1] - eps_vol);
    StressVector[2] = 2.0*mu*(S[2] - eps_vol);
    StressVector[3] = mu*S[3];
    StressVector[4] = mu*S[4];
    StressVector[5] = mu*S[5];

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        Matrix& C                  = rValues.GetConstitutiveMatrix();
        
        noalias(C)  = ZeroMatrix(6,6);

        C(0,0) = 4.0/3.0*mu;  C(0,1) = -2.0/3.0*mu; C(0,2) = -2.0/3.0*mu;
        C(1,2) = -2.0/3.0*mu; C(1,1) = 4.0/3.0*mu;  C(1,2) = -2.0/3.0*mu;
        C(2,2) = -2.0/3.0*mu; C(2,2) = -2.0/3.0*mu; C(2,2) = 4.0/3.0*mu;
        C(3,3) = mu;
        C(4,4) = mu;
        C(5,5) = mu;
            
    }
	  
}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void Newtonian3DLaw::GetLawFeatures(Features& rFeatures)
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

// bool Newtonian3DLaw::CheckParameters(Parameters& rValues)
// {
//     return rValues.CheckAllParameters();
// }



int Newtonian3DLaw::Check(const Properties& rMaterialProperties,
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
