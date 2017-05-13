// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//
// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/linear_plane_stress.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearPlaneStress::LinearPlaneStress()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LinearPlaneStress::LinearPlaneStress(const LinearPlaneStress& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearPlaneStress::Clone() const
{
    LinearPlaneStress::Pointer p_clone(new LinearPlaneStress(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LinearPlaneStress::~LinearPlaneStress()
{
}

void  LinearPlaneStress::CalculateMaterialResponsePK2 (Parameters& rValues)
{
    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();
     
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();    
    Vector& StrainVector                  = rValues.GetStrainVector();
    Vector& StressVector                  = rValues.GetStressVector();
    const double& E          = MaterialProperties[YOUNG_MODULUS];
    const double& NU    = MaterialProperties[POISSON_RATIO];

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if(Options.IsNot( ConstitutiveLaw::COMPUTE_STRAIN )) //large strains
    {

        //1.-Compute total deformation gradient
        const Matrix& F = rValues.GetDeformationGradientF();
        
        Matrix Etensor = prod(trans(F),F);
        Etensor -= IdentityMatrix(2,2);
        Etensor *= 0.5;
        
        noalias(StrainVector) = MathUtils<double>::StrainTensorToVector(Etensor);
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {
        if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

          Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();
          CalculateElasticMatrix( ConstitutiveMatrix, E, NU );
          noalias(StressVector) = prod(ConstitutiveMatrix, StrainVector);
        }
        else {
            
            CalculateStress( StrainVector, StressVector, E, NU );
        }
      
    }
    else if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

        Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix( ConstitutiveMatrix, E, NU );

    }
    
}

//note that since we are in the hypothesis of small strains we can use the same function for everything
void LinearPlaneStress::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}



//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void LinearPlaneStress::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( PLANE_STRESS_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);
	
	//Set the strain size
	rFeatures.mStrainSize = 3;

	//Set the spacedimension
	rFeatures.mSpaceDimension = 2;

}

int LinearPlaneStress::Check(const Properties& rMaterialProperties,
                              const GeometryType& rElementGeometry,
                              const ProcessInfo& rCurrentProcessInfo)
{

    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
        KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ", "" )

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
        KRATOS_THROW_ERROR( std::invalid_argument,"POISSON_RATIO has Key zero invalid value ", "" )


    if(DENSITY.Key() == 0 || rMaterialProperties[DENSITY]<0.00)
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY has Key zero or invalid value ", "" )


    return 0;

}


} // Namespace Kratos
