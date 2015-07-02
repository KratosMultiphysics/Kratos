//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              IPouplana $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/linear_elastic_plastic_3D_law.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearElasticPlastic3DLaw::LinearElasticPlastic3DLaw()
    : HyperElasticPlastic3DLaw()
{

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearElasticPlastic3DLaw::LinearElasticPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
 : HyperElasticPlastic3DLaw(pFlowRule,pYieldCriterion,pHardeningLaw)
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LinearElasticPlastic3DLaw::LinearElasticPlastic3DLaw(const LinearElasticPlastic3DLaw& rOther)
    : HyperElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearElasticPlastic3DLaw::Clone() const
{
    LinearElasticPlastic3DLaw::Pointer p_clone(new LinearElasticPlastic3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LinearElasticPlastic3DLaw::~LinearElasticPlastic3DLaw()
{
}


//*******************************OPERATIONS FROM BASE CLASS***************************
//************************************************************************************

//***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************

double& LinearElasticPlastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if (rThisVariable==DAMAGE_VARIABLE)
    {
        const FlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
        rValue=InternalVariables.DeltaPlasticStrain;
    }

    return( rValue );
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************

void  LinearElasticPlastic3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const ProcessInfo&  CurProcessInfo    = rValues.GetProcessInfo();
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();    
    const GeometryType& DomainGeometry   = rValues.GetElementGeometry();
    
    Vector& StrainVector                  = rValues.GetStrainVector();
    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();
    
    //-----------------------------//

    //0.- Initialize parameters
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize(); //it has to be called at the start

    // Initialize variables from the process information
    ReturnMappingVariables.DeltaTime = CurProcessInfo[DELTA_TIME];
    
    if(CurProcessInfo[IMPLEX] == 1)	
        ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,true);
    else
        ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,false);

    //1.- Lame constants
    const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];

    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
    {
        //only needed 
        const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
		
        //1.-Total Deformation Gradient
        Matrix TotalDeformationGradientF0 = DeformationGradientF;

        //2.-Right Cauchy Green
        Matrix RightCauchyGreen = prod(trans(TotalDeformationGradientF0),TotalDeformationGradientF0);

        //3.-Green-Lagrange Strain: E= 0.5*(FT*F-1)
        this->CalculateGreenLagrangeStrain(RightCauchyGreen,StrainVector);


        if( Options.Is( ConstitutiveLaw::LAST_KNOWN_CONFIGURATION ) || Options.Is( ConstitutiveLaw::FINAL_CONFIGURATION ))
        {
            //2.-Total Deformation Gradient
            Matrix& DeformationGradientF0   = rValues.GetDeformationGradientF0();

            TotalDeformationGradientF0      = prod(DeformationGradientF, DeformationGradientF0);
	  
            //3.-Left Cauchy Green
            Matrix LeftCauchyGreenMatrix = prod(TotalDeformationGradientF0,trans(TotalDeformationGradientF0));
	  
            //4.-Almansi Strain: e= 0.5*(1-invbT*invb)
            this->CalculateAlmansiStrain(LeftCauchyGreenMatrix,StrainVector);
 
            //5.-Pull-back Almansi Strain
            this->TransformStrains(StrainVector,DeformationGradientF,StrainMeasure_Almansi,StrainMeasure_GreenLagrange);
        }
      }

    //2.-Calculate Total PK2 stress
    
    Matrix LinearElasticMatrix = ZeroMatrix(StressVector.size());
    
    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
        LinearElasticMatrix = ConstitutiveMatrix;
        
        this->CalculateStress( StressVector,LinearElasticMatrix,StrainVector,DomainGeometry,ReturnMappingVariables );
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        this->CalculateSecantConstitutiveMatrix( ConstitutiveMatrix, ReturnMappingVariables );
        
        if(ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION))
            this->CalculateTangentConstitutiveMatrix( ConstitutiveMatrix, LinearElasticMatrix, ReturnMappingVariables );
    }

    if( Options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) )
    {
        mpFlowRule->UpdateInternalVariables( ReturnMappingVariables );
    }

}

//************************************************************************************
//************************************************************************************

void LinearElasticPlastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const ProcessInfo&  CurProcessInfo   = rValues.GetProcessInfo();
    const Properties& MaterialProperties = rValues.GetMaterialProperties();   
    const GeometryType& DomainGeometry   = rValues.GetElementGeometry();

    Vector& StrainVector                 = rValues.GetStrainVector();
    Vector& StressVector                 = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix           = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize(); //it has to be called at the start

    // Initialize variables from the process information
    ReturnMappingVariables.DeltaTime = CurProcessInfo[DELTA_TIME];
    
    if(CurProcessInfo[IMPLEX] == 1)	
        ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,true);
    else
        ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,false);

    //1.- Lame constants
    const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];

    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
    {
        //1.-Compute total deformation gradient
        const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();
        Matrix& DeformationGradientF0      = rValues.GetDeformationGradientF0();
        
        Matrix TotalDeformationGradientF0  = prod(DeformationGradientF, DeformationGradientF0);
        TotalDeformationGradientF0         = this->DeformationGradient3D( TotalDeformationGradientF0 );

        //2.-Push-Forward Left Cauchy-Green tensor b to the new configuration
        Matrix LeftCauchyGreenMatrix = prod(TotalDeformationGradientF0,trans(TotalDeformationGradientF0));

        //3.-Almansi Strain: e= 0.5*(1-invFT*invF)
        this->CalculateAlmansiStrain(LeftCauchyGreenMatrix,StrainVector);
    }

    //2.-Calculate Total Kirchhoff stress
    
    Matrix LinearElasticMatrix = ZeroMatrix(StressVector.size());
    
    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
        LinearElasticMatrix = ConstitutiveMatrix;
        
        this->CalculateStress( StressVector,LinearElasticMatrix,StrainVector,DomainGeometry,ReturnMappingVariables );
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        this->CalculateSecantConstitutiveMatrix( ConstitutiveMatrix, ReturnMappingVariables );
        
        if(ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION))
            this->CalculateTangentConstitutiveMatrix( ConstitutiveMatrix, LinearElasticMatrix, ReturnMappingVariables );
    }

    if( Options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) )
    {
        mpFlowRule->UpdateInternalVariables( ReturnMappingVariables );
    }
}


//********************* COMPUTE LINEAR ELASTIC CONSTITUTIVE MATRIX *******************
//************************************************************************************


void LinearElasticPlastic3DLaw::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,const double &rYoungModulus,const double &rPoissonCoefficient )
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


//******************************* COMPUTE STRESS *************************************
//************************************************************************************

void LinearElasticPlastic3DLaw::CalculateStress( Vector& rStressVector, const Matrix& rLinearElasticMatrix, const Vector& rStrainVector,
                                                 const GeometryType& rDomainGeometry, FlowRule::RadialReturnVariables& rReturnMappingVariables )
{
    noalias(rStressVector) = prod(rLinearElasticMatrix, rStrainVector);
    Matrix StressMatrix = MathUtils<double>::StressVectorToTensor(rStressVector);
    rReturnMappingVariables.TrialIsoStressMatrix = StressMatrix;
    
    Matrix StrainMatrix = MathUtils<double>::StrainVectorToTensor(rStrainVector);
    rReturnMappingVariables.StrainMatrix = StrainMatrix;
    
    //Element characteristic size
    double CharacteristicSize = 1.0;
    this->CalculateCharacteristicSize(CharacteristicSize,rDomainGeometry);
    rReturnMappingVariables.CharacteristicSize = CharacteristicSize;
    
    mpFlowRule->CalculateReturnMapping( rReturnMappingVariables, StressMatrix );
    
    rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );
}


//********************** COMPUTE ELEMENT CHARACTERISTIC SIZE *************************
//************************************************************************************

void LinearElasticPlastic3DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& rDomainGeometry )
{

}


//********************** COMPUTE SECANT CONSTITUTIVE MATRIX **************************
//************************************************************************************

void LinearElasticPlastic3DLaw::CalculateSecantConstitutiveMatrix( Matrix& rConstitutiveMatrix, FlowRule::RadialReturnVariables& rReturnMappingVariables )
{

}


//********************** COMPUTE TANGENT CONSTITUTIVE MATRIX **************************
//************************************************************************************

void LinearElasticPlastic3DLaw::CalculateTangentConstitutiveMatrix( Matrix& rConstitutiveMatrix, const Matrix& rLinearElasticMatrix, 
                                                                    FlowRule::RadialReturnVariables& rReturnMappingVariables )
{
    double Alpha = 1.0;
    Matrix TangentConstitutiveMatrix = ZeroMatrix(rConstitutiveMatrix.size1());
    
    mpFlowRule->ComputeElastoPlasticTangentMatrix(rReturnMappingVariables, rLinearElasticMatrix, Alpha, TangentConstitutiveMatrix);
    
    noalias(rConstitutiveMatrix) += TangentConstitutiveMatrix;
}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void LinearElasticPlastic3DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

} // Namespace Kratos
