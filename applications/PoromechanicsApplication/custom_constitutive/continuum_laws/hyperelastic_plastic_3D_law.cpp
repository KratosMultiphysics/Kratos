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
#include "custom_constitutive/continuum_laws/hyperelastic_plastic_3D_law.hpp"

#include "poromechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlastic3DLaw::HyperElasticPlastic3DLaw()
    : HyperElastic3DLaw()
{

}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlastic3DLaw::HyperElasticPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : HyperElastic3DLaw()
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  pYieldCriterion;
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticPlastic3DLaw::HyperElasticPlastic3DLaw(const HyperElasticPlastic3DLaw& rOther)
    : HyperElastic3DLaw(rOther)
    ,mElasticLeftCauchyGreen(rOther.mElasticLeftCauchyGreen)
    ,mpYieldCriterion(rOther.mpYieldCriterion)
    ,mpHardeningLaw(rOther.mpHardeningLaw)
{

  mpFlowRule       = rOther.mpFlowRule->Clone();
  //they not contain member variables to be stored and cloned:
  //mpYieldCriterion = rOther.mpYieldCriterion->Clone();
  //mpHardeningLaw   = rOther.mpHardeningLaw->Clone();

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticPlastic3DLaw::Clone() const
{
    return Kratos::make_shared<HyperElasticPlastic3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlastic3DLaw::~HyperElasticPlastic3DLaw()
{
}


//*******************************OPERATIONS FROM BASE CLASS***************************
//************************************************************************************

//***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
//************************************************************************************

bool HyperElasticPlastic3DLaw::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool HyperElasticPlastic3DLaw::Has( const Variable<Vector>& rThisVariable )
{
    return false;
}

bool HyperElasticPlastic3DLaw::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}


//******************CALCULATE VALUE: DOUBLE - VECTOR - MATRIX*************************
//************************************************************************************

double& HyperElasticPlastic3DLaw::CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue )
{

  return (this->GetValue(rThisVariable,rValue ));

}

//***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************

double& HyperElasticPlastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if (rThisVariable == DETERMINANT_F)
    {
        rValue = mDeterminantF0;
    }

    if (rThisVariable == PLASTIC_STRAIN)
    {
        const FlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
        rValue = InternalVariables.EquivalentPlasticStrain;
    }

    if (rThisVariable == DELTA_PLASTIC_STRAIN)
    {
        const FlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
        rValue = InternalVariables.DeltaPlasticStrain;
    }

    return( rValue );
}

Vector& HyperElasticPlastic3DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return( rValue );
}

Matrix& HyperElasticPlastic3DLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return( rValue );
}


//***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************


void HyperElasticPlastic3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                  const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == DETERMINANT_F)
    {
        mDeterminantF0 = rValue;
    }

    if (rThisVariable == PLASTIC_STRAIN)
    {
        mpFlowRule->SetInternalVariables().EquivalentPlasticStrain = rValue;
    }

    if (rThisVariable == DELTA_PLASTIC_STRAIN)
    {
        mpFlowRule->SetInternalVariables().DeltaPlasticStrain = rValue;
    }

}

void HyperElasticPlastic3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                  const ProcessInfo& rCurrentProcessInfo )
{

}

void HyperElasticPlastic3DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                  const ProcessInfo& rCurrentProcessInfo )
{

}



//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************


void HyperElasticPlastic3DLaw::InitializeMaterial( const Properties& rMaterialProperties,
						   const GeometryType& rElementGeometry,
						   const Vector& rShapeFunctionsValues )
{

  HyperElastic3DLaw::InitializeMaterial(rMaterialProperties,rElementGeometry,rShapeFunctionsValues);

  mElasticLeftCauchyGreen       = identity_matrix<double> (3);

  mpHardeningLaw->SetProperties(rMaterialProperties);

  mpFlowRule->InitializeMaterial( mpYieldCriterion, mpHardeningLaw, rMaterialProperties );
  //mpFlowRule->InitializeMaterial( rMaterialProperties );
}


//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************


void  HyperElasticPlastic3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
{

  this->CalculateMaterialResponseKirchhoff (rValues);

  //1.- Obtain parameters
  Flags & Options                    = rValues.GetOptions();

  Vector& StressVector               = rValues.GetStressVector();
  Vector& StrainVector               = rValues.GetStrainVector();

  const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
  const double& DeterminantF         = rValues.GetDeterminantF();

  Matrix& ConstitutiveMatrix         = rValues.GetConstitutiveMatrix();

  //2.-Green-Lagrange Strain:
  if(Options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
      TransformStrains (StrainVector, DeformationGradientF, StrainMeasure_Almansi, StrainMeasure_GreenLagrange);
    }

  //3.-Calculate Total PK2 stress
  if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {
      TransformStresses(StressVector, DeformationGradientF, DeterminantF, StressMeasure_Kirchhoff, StressMeasure_PK2);
    }

  //4.-Calculate PK2 constitutive tensor
  if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
      PullBackConstitutiveMatrix(ConstitutiveMatrix, DeformationGradientF);
    }

}


//************************************************************************************
//************************************************************************************


void HyperElasticPlastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const ProcessInfo& rCurrentProcessInfo = rValues.GetProcessInfo();
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();

    const Matrix& DeformationGradientF    = rValues.GetDeformationGradientF();
    const double& DeterminantF            = rValues.GetDeterminantF();

    const GeometryType&  DomainGeometry   = rValues.GetElementGeometry ();
    const Vector&        ShapeFunctions   = rValues.GetShapeFunctionsValues ();

    Vector& StrainVector                  = rValues.GetStrainVector();
    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.Identity = identity_matrix<double> ( 3 );

    ElasticVariables.SetElementGeometry(DomainGeometry);
    ElasticVariables.SetShapeFunctionsValues(ShapeFunctions);

    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize(); //it has to be called at the start

    // Initialize variables from the process information
    ReturnMappingVariables.DeltaTime = rCurrentProcessInfo[DELTA_TIME];

    if(rCurrentProcessInfo[IMPLEX] == 1)
      ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,true);
    else
      ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,false);

    // Initialize Splited Parts: Isochoric and Volumetric stresses and constitutive tensors
    double voigtsize = StressVector.size();
    VectorSplit SplitStressVector;
    MatrixSplit SplitConstitutiveMatrix;

    //1.- Lame constants
    const double& YoungModulus       = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];

    ElasticVariables.LameLambda      = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    ElasticVariables.LameMu          =  YoungModulus/(2*(1+PoissonCoefficient));

    //1.1- Thermal constants
    if( MaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT) )
      ElasticVariables.ThermalExpansionCoefficient = MaterialProperties[THERMAL_EXPANSION_COEFFICIENT];
    else
      ElasticVariables.ThermalExpansionCoefficient = 0;

    if( MaterialProperties.Has(REFERENCE_TEMPERATURE) )
      ElasticVariables.ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
    else
      ElasticVariables.ReferenceTemperature = 0;


    //2.-Determinant of the Total DeformationGradientF
    ElasticVariables.DeterminantF = DeterminantF;

    //3.-Compute Incremental DeformationGradientF_bar
    double detF = DeterminantF / mDeterminantF0;

    ElasticVariables.J_pow13 = pow(detF,1.0/3.0);

    ElasticVariables.DeformationGradientF = DeformationGradientF;

    ElasticVariables.DeformationGradientF = this->Transform2DTo3D(ElasticVariables.DeformationGradientF);

    ElasticVariables.DeformationGradientF = prod(ElasticVariables.DeformationGradientF, this->mInverseDeformationGradientF0);

    ElasticVariables.DeformationGradientF /= ElasticVariables.J_pow13; //now ElasticVariables.DeformationGradientF is DeformationGradientFbar

    //4.-Left Cauchy-Green tensor b_bar to the new configuration
    ElasticVariables.CauchyGreenMatrix.resize(3,3,false);
    noalias(ElasticVariables.CauchyGreenMatrix) = prod(mElasticLeftCauchyGreen,trans(ElasticVariables.DeformationGradientF));
    ElasticVariables.CauchyGreenMatrix = prod(ElasticVariables.DeformationGradientF,ElasticVariables.CauchyGreenMatrix);


    //5.-Calculate trace of Left Cauchy-Green tensor b_bar
    ElasticVariables.traceCG = 0;
    for( unsigned int i=0; i<3; i++)
    {
       ElasticVariables.traceCG += ElasticVariables.CauchyGreenMatrix( i , i );
    }

    ReturnMappingVariables.LameMu_bar = ElasticVariables.LameMu * ( ElasticVariables.traceCG / 3.0  );

    //4.-Almansi Strain:
    if(Options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
        // e= 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix,StrainVector);
	// correct b_bar to b
	double J_pow23 = pow(ElasticVariables.DeterminantF,2.0/3.0);
	StrainVector /= (J_pow23*J_pow23);
    }


    //5.-Calculate Total Kirchhoff stress
    SplitStressVector.Isochoric.resize(voigtsize,false);
    noalias(SplitStressVector.Isochoric)  = ZeroVector(voigtsize);
    Matrix IsochoricStressMatrix(3,3);
    noalias(IsochoricStressMatrix) = ZeroMatrix(3,3);

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
      this->CalculatePlasticIsochoricStress( ElasticVariables, ReturnMappingVariables, StressMeasure_Kirchhoff, IsochoricStressMatrix, SplitStressVector.Isochoric );

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {

        SplitStressVector.Volumetric.resize(voigtsize,false);
        noalias(SplitStressVector.Volumetric) = ZeroVector(voigtsize);

        ElasticVariables.CauchyGreenMatrix = ElasticVariables.Identity;

        this->CalculateVolumetricStress ( ElasticVariables, SplitStressVector.Volumetric );

        //Kirchhoff Stress:
        StressVector = SplitStressVector.Isochoric + SplitStressVector.Volumetric;

	// if( ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION) ){
	//   std::cout<<" StressVector.Isochoric "<<SplitStressVector.Isochoric<<std::endl;
	//   std::cout<<" StressVector.Volumetric "<<SplitStressVector.Volumetric<<std::endl;
	// }

        if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) )
        {
            StressVector = SplitStressVector.Isochoric;
        }
        else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) )
        {
            StressVector = SplitStressVector.Volumetric;
        }

    }


    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

        if( ReturnMappingVariables.Options.IsNot(FlowRule::RETURN_MAPPING_COMPUTED) )
        {
            KRATOS_ERROR << " ReturnMappingCall was not performed  ...error in the constitutive calculation..." << std::endl;
        }

        //initialize constitutive tensors
        ConstitutiveMatrix.clear();
        SplitConstitutiveMatrix.Isochoric  = ConstitutiveMatrix;
        SplitConstitutiveMatrix.Volumetric = ConstitutiveMatrix;
	SplitConstitutiveMatrix.Plastic    = ConstitutiveMatrix;

        ElasticVariables.CauchyGreenMatrix = ElasticVariables.Identity;

	this->CalculateIsochoricConstitutiveMatrix  ( ElasticVariables, ReturnMappingVariables.TrialIsoStressMatrix, SplitConstitutiveMatrix.Isochoric );

	this->CalculateVolumetricConstitutiveMatrix ( ElasticVariables, SplitConstitutiveMatrix.Volumetric );

	if( ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION) )
	  this->CalculatePlasticConstitutiveMatrix  ( ElasticVariables, ReturnMappingVariables, SplitConstitutiveMatrix.Plastic );


	// std::cout<< " Isochoric Constitutive "<<SplitConstitutiveMatrix.Isochoric<<std::endl;
	// std::cout<< " Volumetric Constitutive "<<SplitConstitutiveMatrix.Volumetric<<std::endl;
	//if( ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION) )
	  //std::cout<< " Plastic Constitutive   "<<SplitConstitutiveMatrix.Plastic<<std::endl;

	ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric + SplitConstitutiveMatrix.Volumetric + SplitConstitutiveMatrix.Plastic;

        if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) )
        {
            ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric + SplitConstitutiveMatrix.Plastic;
        }
        else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) )
        {
            ConstitutiveMatrix = SplitConstitutiveMatrix.Volumetric;
        }
    }



    if( Options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) )
    {
      mpFlowRule->UpdateInternalVariables ( ReturnMappingVariables );

      mElasticLeftCauchyGreen  = ( IsochoricStressMatrix * ( 1.0 / ElasticVariables.LameMu ) );
      mElasticLeftCauchyGreen += ( ElasticVariables.traceCG/3.0) * ElasticVariables.Identity;

    }



    // std::cout<<" StrainVector "<<StrainVector<<std::endl;
    // std::cout<<" StressVector "<<StressVector<<std::endl;
    // std::cout<<" ConstitutiveMatrix "<<ConstitutiveMatrix<<std::endl;


}


//******************************* COMPUTE TOTAL STRESS  ******************************
//************************************************************************************


//******************************* COMPUTE ISOCHORIC STRESS  **************************
//************************************************************************************

void HyperElasticPlastic3DLaw::CalculatePlasticIsochoricStress( MaterialResponseVariables & rElasticVariables,
							 FlowRule::RadialReturnVariables & rReturnMappingVariables,
							 StressMeasure rStressMeasure,
							 Matrix& rIsoStressMatrix,
							 Vector& rIsoStressVector)
{

    //note.- rElasticVariables.traceCG is "traceCG_bar"

    if(rStressMeasure == StressMeasure_PK2)
    {

        //rElasticVariables.CauchyGreenMatrix is InverseRightCauchyGreen

        //2.-Isochoric part of the 2nd Piola Kirchhoff Stress Matrix
        rIsoStressMatrix  = (rElasticVariables.Identity - (rElasticVariables.traceCG/3.0)*rElasticVariables.CauchyGreenMatrix );
        rIsoStressMatrix *= rElasticVariables.LameMu;

        //std::cout<<" PK2 "<<std::endl;
    }

    if(rStressMeasure == StressMeasure_Kirchhoff)
    {

        //rElasticVariables.CauchyGreenMatrix is LeftCauchyGreen

        //2.-Isochoric part of the Kirchhoff Stress Matrix
        rIsoStressMatrix  = (rElasticVariables.CauchyGreenMatrix - (rElasticVariables.traceCG/3.0)*rElasticVariables.Identity );
        rIsoStressMatrix *= rElasticVariables.LameMu;

        //std::cout<<" Kirchooff "<<IsoStressMatrix<<std::endl;

    }


    //thermal effects:
    rReturnMappingVariables.Temperature = this->CalculateDomainTemperature(rElasticVariables, rReturnMappingVariables.Temperature);

    rReturnMappingVariables.TrialIsoStressMatrix = rIsoStressMatrix;

    //std::cout<<" TrialIsoStressMatrix "<<rReturnMappingVariables.TrialIsoStressMatrix<<std::endl;

    mpFlowRule->CalculateReturnMapping( rReturnMappingVariables, rIsoStressMatrix );

    rIsoStressVector = MathUtils<double>::StressTensorToVector( rIsoStressMatrix, rIsoStressVector.size() );

    //std::cout<<" PLASTICITY "<<rElasticVariables.Plasticity<<" rIsoStressVector "<<rIsoStressVector<<std::endl;

}


//***********************COMPUTE PLASTIC CONSTITUTIVE MATRIX**************************
//************************************************************************************

void HyperElasticPlastic3DLaw::CalculatePlasticConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
								   FlowRule::RadialReturnVariables & rReturnMappingVariables,
								   Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    Matrix IsoStressMatrix = rReturnMappingVariables.TrialIsoStressMatrix;

    //std::cout<< " TrialStressMatrix 3D "<<IsoStressMatrix<<std::endl;

    FlowRule::PlasticFactors ScalingFactors;
    mpFlowRule->CalculateScalingFactors( rReturnMappingVariables, ScalingFactors );


    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
		rConstitutiveMatrix( i, j ) = PlasticConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, IsoStressMatrix, ScalingFactors,
                                          this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}




//********************CONSTITUTIVE MATRIX PLASTIC COMPONENT***************************
//************************************************************************************


double& HyperElasticPlastic3DLaw::PlasticConstitutiveComponent(double & rCabcd,
							       const MaterialResponseVariables & rElasticVariables,
							       const Matrix & rIsoStressMatrix,
							       const FlowRule::PlasticFactors & rScalingFactors,
							       const unsigned int& a, const unsigned int& b,
							       const unsigned int& c, const unsigned int& d)
{

    //Plastic part of the algorithmic moduli

    rCabcd  = (1.0/3.0)*(rElasticVariables.CauchyGreenMatrix(a,b)*rElasticVariables.CauchyGreenMatrix(c,d));
    rCabcd -= (0.5*(rElasticVariables.CauchyGreenMatrix(a,c)*rElasticVariables.CauchyGreenMatrix(b,d)+rElasticVariables.CauchyGreenMatrix(a,d)*rElasticVariables.CauchyGreenMatrix(b,c)));

    rCabcd *= rElasticVariables.traceCG * rElasticVariables.LameMu;

    rCabcd += (rElasticVariables.CauchyGreenMatrix(c,d)*rIsoStressMatrix(a,b) + rIsoStressMatrix(c,d)*rElasticVariables.CauchyGreenMatrix(a,b));

    rCabcd *= (-2.0/3.0) * ( (-1) * rScalingFactors.Beta1 );

    rCabcd  -= rScalingFactors.Beta3 * 2.0 * (rElasticVariables.LameMu * (rElasticVariables.traceCG/3.0)) * ( rScalingFactors.Normal(a,b) * rScalingFactors.Normal(c,d) );

    rCabcd  -= rScalingFactors.Beta4 * 2.0 * (rElasticVariables.LameMu * (rElasticVariables.traceCG/3.0)) * ( rScalingFactors.Normal(a,b) * rScalingFactors.Dev_Normal(c,d) );

    return rCabcd;
}



//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElasticPlastic3DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( FINITE_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}


//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
//************************************************************************************

bool HyperElasticPlastic3DLaw::CheckParameters(Parameters& rValues)
{
    return rValues.CheckAllParameters();
}



int HyperElasticPlastic3DLaw::Check(const Properties& rMaterialProperties,
                             const GeometryType& rElementGeometry,
                             const ProcessInfo& rCurrentProcessInfo) const
{

    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
    {
        KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;
    }

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
    {
        KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value " << std::endl;
    }


    if(DENSITY.Key() == 0 || rMaterialProperties[DENSITY]<0.00)
    {
        KRATOS_ERROR << "DENSITY has Key zero or invalid value " << std::endl;
    }


    return 0;

}

} // Namespace Kratos
