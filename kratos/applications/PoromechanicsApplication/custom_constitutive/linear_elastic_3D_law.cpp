//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

/* Project includes */
#include "custom_constitutive/linear_elastic_3D_law.hpp"

#include "poromechanics_application.h"

namespace Kratos
{

//Default Constructor
LinearElastic3DLaw::LinearElastic3DLaw() : ConstitutiveLaw() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
LinearElastic3DLaw::LinearElastic3DLaw(const LinearElastic3DLaw& rOther) : ConstitutiveLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
LinearElastic3DLaw::~LinearElastic3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LinearElastic3DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the spacedimension
	rFeatures.mSpaceDimension = 3;
    
	//Set the strain size
	rFeatures.mStrainSize = 6;
}

//----------------------------------------------------------------------------------------

int LinearElastic3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo)
{
    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
        KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ", "" )

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
        KRATOS_THROW_ERROR( std::invalid_argument,"POISSON_RATIO has Key zero invalid value ", "" )

    if(DENSITY_SOLID.Key() == 0 || rMaterialProperties[DENSITY_SOLID]<0.00)
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY_SOLID has Key zero or invalid value ", "" )

    if(DENSITY_WATER.Key() == 0 || rMaterialProperties[DENSITY_WATER]<0.00)
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY_WATER has Key zero or invalid value ", "" )

    return 0;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer LinearElastic3DLaw::Clone() const
{
    LinearElastic3DLaw::Pointer p_clone(new LinearElastic3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------

void LinearElastic3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LinearElastic3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    rValues.CheckAllParameters();

    Vector& StrainVector = rValues.GetStrainVector();
    Vector& StressVector = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();

    const double& YoungModulus = rValues.GetMaterialProperties()[YOUNG_MODULUS];
    const double& PoissonCoefficient = rValues.GetMaterialProperties()[POISSON_RATIO];

    if( rValues.GetOptions().Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {        
        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
        this->CalculateStress(StressVector,StrainVector,ConstitutiveMatrix);
    }
    else if(rValues.GetOptions().Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ))
    {
        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LinearElastic3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double& LinearElastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    return rValue;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LinearElastic3DLaw::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,const double &rYoungModulus,const double &rPoissonCoefficient )
{
    rConstitutiveMatrix.clear();
    
    rConstitutiveMatrix ( 0 , 0 ) = (rYoungModulus*(1.0-rPoissonCoefficient)/((1.0+rPoissonCoefficient)*(1.0-2.0*rPoissonCoefficient)));
    rConstitutiveMatrix ( 1 , 1 ) = rConstitutiveMatrix ( 0 , 0 );
    rConstitutiveMatrix ( 2 , 2 ) = rConstitutiveMatrix ( 0 , 0 );

    rConstitutiveMatrix ( 3 , 3 ) = rConstitutiveMatrix ( 0 , 0 )*(1-2*rPoissonCoefficient)/(2.0*(1.0-rPoissonCoefficient));
    rConstitutiveMatrix ( 4 , 4 ) = rConstitutiveMatrix ( 3 , 3 );
    rConstitutiveMatrix ( 5 , 5 ) = rConstitutiveMatrix ( 3 , 3 );

    rConstitutiveMatrix ( 0 , 1 ) = rConstitutiveMatrix ( 0 , 0 )*rPoissonCoefficient/(1.0-rPoissonCoefficient);
    rConstitutiveMatrix ( 1 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

    rConstitutiveMatrix ( 0 , 2 ) = rConstitutiveMatrix ( 0 , 1 );
    rConstitutiveMatrix ( 2 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

    rConstitutiveMatrix ( 1 , 2 ) = rConstitutiveMatrix ( 0 , 1 );
    rConstitutiveMatrix ( 2 , 1 ) = rConstitutiveMatrix ( 0 , 1 );
}

void LinearElastic3DLaw::CalculateStress( Vector& rStressVector, const Vector& rStrainVector, const Matrix& rConstitutiveMatrix )
{
    rStressVector = prod(rConstitutiveMatrix,rStrainVector);
}

} // Namespace Kratos
