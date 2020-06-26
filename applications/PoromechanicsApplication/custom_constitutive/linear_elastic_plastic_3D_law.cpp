//
//   Project Name:        KratosPoromechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_elastic_plastic_3D_law.hpp"

#include "poromechanics_application_variables.h"

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
    return Kratos::make_shared<LinearElasticPlastic3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LinearElasticPlastic3DLaw::~LinearElasticPlastic3DLaw()
{
}


//*******************************OPERATIONS FROM BASE CLASS***************************
//************************************************************************************


//******************CALCULATE VALUE: DOUBLE - VECTOR - MATRIX*************************
//************************************************************************************

double& LinearElasticPlastic3DLaw::CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue )
{

  return (this->GetValue(rThisVariable,rValue ));

}

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
    Flags& Options                        = rValues.GetOptions();
    const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    Vector& rStrainVector                 = rValues.GetStrainVector();

    //-----------------------------//

    if(Options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {

        //1.-Deformation Gradient
        const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();

        //2.-Right Cauchy Green
        Matrix RightCauchyGreen = prod(trans(DeformationGradientF),DeformationGradientF);

        //3.-Green-Lagrange Strain: E= 0.5*(FT*F-1)
        this->CalculateGreenLagrangeStrain(RightCauchyGreen,rStrainVector);

    }

    //0.- Initialize parameters
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize(); //it has to be called at the start

    // Initialize variables from the process information
    ReturnMappingVariables.DeltaTime = CurrentProcessInfo[DELTA_TIME];

    if(CurrentProcessInfo[IMPLEX] == 1)
        ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,true);
    else
        ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,false);

    // Strain and Stress matrices
    const unsigned int Dim = this->WorkingSpaceDimension();
    Matrix AuxMatrix(Dim,Dim);
    noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(rStrainVector);
    ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
    noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
    ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);

    // CharacteristicSize
    double CharacteristicSize = 1.0;
    this->CalculateCharacteristicSize(CharacteristicSize,rValues.GetElementGeometry());
    ReturnMappingVariables.CharacteristicSize = CharacteristicSize;

    //1.- Lame constants
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    const unsigned int VoigtSize = rStrainVector.size();
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);

    if(Options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {

        //1.-Deformation Gradient
        const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();

        //2.-Right Cauchy Green
        Matrix RightCauchyGreen = prod(trans(DeformationGradientF),DeformationGradientF);

        //3.-Green-Lagrange Strain: E= 0.5*(FT*F-1)
        this->CalculateGreenLagrangeStrain(RightCauchyGreen,rStrainVector);

    }

    //2.-Calculate Total PK2 stress

    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {
        if(Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS))
        {
            // COMPUTE_CONSTITUTIVE_TENSOR
            Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            Vector EffectiveStressVector(VoigtSize);

            this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,EffectiveStressVector,LinearElasticMatrix,rStrainVector);

            this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
        }
        else
        {
            // COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
            Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            Vector& rStressVector = rValues.GetStressVector();

            this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);

            this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
        }
    }
    else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS) && Options.IsNot( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ))
    {
        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();

        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
    }

    if( Options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) )
    {
        Vector EffectiveStressVector(VoigtSize);

        this->UpdateInternalStateVariables(ReturnMappingVariables,EffectiveStressVector,LinearElasticMatrix,rStrainVector);
        if (Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
        {
            Vector& rStressVector = rValues.GetStressVector();
            this->UpdateStressVector(rStressVector,ReturnMappingVariables,EffectiveStressVector);
        }
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
    Flags& Options                        = rValues.GetOptions();
    const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    Vector& rStrainVector                 = rValues.GetStrainVector();

    //-----------------------------//

    if(Options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
        //1.-Compute total deformation gradient
        const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();

        //2.-Push-Forward Left Cauchy-Green tensor b to the new configuration
        Matrix LeftCauchyGreenMatrix = prod(DeformationGradientF,trans(DeformationGradientF));

        //3.-Almansi Strain: e= 0.5*(1-invFT*invF)
        this->CalculateAlmansiStrain(LeftCauchyGreenMatrix,rStrainVector);
    }

    //0.- Initialize parameters
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize(); //it has to be called at the start

    // Initialize variables from the process information
    ReturnMappingVariables.DeltaTime = CurrentProcessInfo[DELTA_TIME];

    if(CurrentProcessInfo[IMPLEX] == 1)
        ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,true);
    else
        ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,false);

    // Strain and Stress matrices
    const unsigned int Dim = this->WorkingSpaceDimension();
    Matrix AuxMatrix(Dim,Dim);
    noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(rStrainVector);
    ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
    noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
    ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);

    // CharacteristicSize
    double CharacteristicSize = 1.0;
    this->CalculateCharacteristicSize(CharacteristicSize,rValues.GetElementGeometry());
    ReturnMappingVariables.CharacteristicSize = CharacteristicSize;

    //1.- Lame constants
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    const unsigned int VoigtSize = rStrainVector.size();
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);


    if(Options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
        //1.-Compute total deformation gradient
        const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();

        //2.-Push-Forward Left Cauchy-Green tensor b to the new configuration
        Matrix LeftCauchyGreenMatrix = prod(DeformationGradientF,trans(DeformationGradientF));

        //3.-Almansi Strain: e= 0.5*(1-invFT*invF)
        this->CalculateAlmansiStrain(LeftCauchyGreenMatrix,rStrainVector);
    }

    //2.-Calculate Total Kirchhoff stress

    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {
        if(Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS))
        {
            // COMPUTE_CONSTITUTIVE_TENSOR
            Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            Vector EffectiveStressVector(VoigtSize);

            this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,EffectiveStressVector,LinearElasticMatrix,rStrainVector);

            this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
        }
        else
        {
            // COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
            Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            Vector& rStressVector = rValues.GetStressVector();

            this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);

            this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
        }
    }
    else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS) && Options.IsNot( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ))
    {
        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();

        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
    }

    if( Options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) )
    {
        Vector EffectiveStressVector(VoigtSize);

        this->UpdateInternalStateVariables(ReturnMappingVariables,EffectiveStressVector,LinearElasticMatrix,rStrainVector);
        if (Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
        {
            Vector& rStressVector = rValues.GetStressVector();
            this->UpdateStressVector(rStressVector,ReturnMappingVariables,EffectiveStressVector);
        }
    }
}


//********************** COMPUTE ELEMENT CHARACTERISTIC SIZE *************************
//************************************************************************************

void LinearElasticPlastic3DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry )
{

}


//********************* COMPUTE LINEAR ELASTIC CONSTITUTIVE MATRIX *******************
//************************************************************************************


void LinearElasticPlastic3DLaw::CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,const double& YoungModulus,const double& PoissonCoefficient )
{
    rLinearElasticMatrix.clear();

    // 3D linear elastic constitutive matrix
    rLinearElasticMatrix ( 0 , 0 ) = (YoungModulus*(1.0-PoissonCoefficient)/((1.0+PoissonCoefficient)*(1.0-2.0*PoissonCoefficient)));
    rLinearElasticMatrix ( 1 , 1 ) = rLinearElasticMatrix ( 0 , 0 );
    rLinearElasticMatrix ( 2 , 2 ) = rLinearElasticMatrix ( 0 , 0 );

    rLinearElasticMatrix ( 3 , 3 ) = rLinearElasticMatrix ( 0 , 0 )*(1.0-2.0*PoissonCoefficient)/(2.0*(1.0-PoissonCoefficient));
    rLinearElasticMatrix ( 4 , 4 ) = rLinearElasticMatrix ( 3 , 3 );
    rLinearElasticMatrix ( 5 , 5 ) = rLinearElasticMatrix ( 3 , 3 );

    rLinearElasticMatrix ( 0 , 1 ) = rLinearElasticMatrix ( 0 , 0 )*PoissonCoefficient/(1.0-PoissonCoefficient);
    rLinearElasticMatrix ( 1 , 0 ) = rLinearElasticMatrix ( 0 , 1 );

    rLinearElasticMatrix ( 0 , 2 ) = rLinearElasticMatrix ( 0 , 1 );
    rLinearElasticMatrix ( 2 , 0 ) = rLinearElasticMatrix ( 0 , 1 );

    rLinearElasticMatrix ( 1 , 2 ) = rLinearElasticMatrix ( 0 , 1 );
    rLinearElasticMatrix ( 2 , 1 ) = rLinearElasticMatrix ( 0 , 1 );
}


//********************************* RETURN MAPPING ***********************************
//************************************************************************************

void LinearElasticPlastic3DLaw::CalculateReturnMapping( FlowRule::RadialReturnVariables& rReturnMappingVariables, Matrix& rStressMatrix,
                                                        Vector& rStressVector, const Matrix& LinearElasticMatrix, const Vector& StrainVector )
{
    noalias(rStressVector) = prod(LinearElasticMatrix, StrainVector);
    noalias(rReturnMappingVariables.TrialIsoStressMatrix) = MathUtils<double>::StressVectorToTensor(rStressVector);

    mpFlowRule->CalculateReturnMapping( rReturnMappingVariables, rStressMatrix );

    noalias(rStressVector) = MathUtils<double>::StressTensorToVector( rStressMatrix, StrainVector.size() );
}


//*************************** COMPUTE CONSTITUTIVE MATRIX ****************************
//************************************************************************************

void LinearElasticPlastic3DLaw::CalculateConstitutiveTensor( Matrix& rConstitutiveMatrix, FlowRule::RadialReturnVariables& rReturnMappingVariables,
                                                                const Matrix& LinearElasticMatrix )
{
    // Calculate secant component of the constitutive matrix
    noalias(rConstitutiveMatrix) = (1.0 - rReturnMappingVariables.TrialStateFunction)*LinearElasticMatrix; // Csec = (1-d)*Ce

    // Calculate tangent component if necessary
    if(rReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION))
    {
        double Alpha = 1.0;

        mpFlowRule->ComputeElastoPlasticTangentMatrix(rReturnMappingVariables, LinearElasticMatrix, Alpha, rConstitutiveMatrix);
    }
}


//**************************** UPDATE INTERNAL VARIABLES *****************************
//************************************************************************************

void LinearElasticPlastic3DLaw::UpdateInternalStateVariables( FlowRule::RadialReturnVariables& rReturnMappingVariables,Vector& rStressVector,
                                                            const Matrix& LinearElasticMatrix, const Vector& StrainVector )
{
    noalias(rStressVector) = prod(LinearElasticMatrix, StrainVector);
    noalias(rReturnMappingVariables.TrialIsoStressMatrix) = MathUtils<double>::StressVectorToTensor(rStressVector);

    mpFlowRule->UpdateInternalVariables( rReturnMappingVariables );
}

//**************************** UPDATE STRESS VECTOR **********************************
//************************************************************************************

void LinearElasticPlastic3DLaw::UpdateStressVector( Vector& rStressVector, FlowRule::RadialReturnVariables& rReturnMappingVariables,
                                                    const Vector& EffectiveStressVector )
{
    noalias(rStressVector) = (1.0 - rReturnMappingVariables.TrialStateFunction)*EffectiveStressVector; // S = (1-d)*Se
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
