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
#include "custom_constitutive/hyperelastic_plastic_3D_law.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlastic3DLaw::HyperElasticPlastic3DLaw()
    : ConstitutiveLaw()
{

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlastic3DLaw::HyperElasticPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : ConstitutiveLaw()
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  pYieldCriterion;
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticPlastic3DLaw::HyperElasticPlastic3DLaw(const HyperElasticPlastic3DLaw& rOther)
    : ConstitutiveLaw(rOther)
  ,mElasticLeftCauchyGreen(rOther.mElasticLeftCauchyGreen)
{

  mpFlowRule       = rOther.mpFlowRule->Clone();
  mpYieldCriterion = rOther.mpYieldCriterion->Clone();
  mpHardeningLaw   = rOther.mpHardeningLaw->Clone();

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticPlastic3DLaw::Clone() const
{
    ConstitutiveLaw::Pointer p_clone(new HyperElasticPlastic3DLaw(*this));
    return p_clone;
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
  // if(rThisVariable == DELTA_PLASTIC_DISSIPATION || rThisVariable == PLASTIC_DISSIPATION )
  //   return true;
  
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


//***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************

double& HyperElasticPlastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
  if (rThisVariable==PLASTIC_STRAIN)
    {
      const FlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
      rValue=InternalVariables.EquivalentPlasticStrain;
    }
  
  if (rThisVariable==DELTA_PLASTIC_STRAIN)
    {
      const FlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
      rValue=InternalVariables.DeltaPlasticStrain;
    }


  // if (rThisVariable==PLASTIC_DISSIPATION)
  //   {
  //     const FlowRule::ThermalVariables& ThermalVariables = mpFlowRule->GetThermalVariables();
  //     rValue=ThermalVariables.PlasticDissipation;
  //   }
  
  // if (rThisVariable==DELTA_PLASTIC_DISSIPATION)
  //   {
  //     const FlowRule::ThermalVariables& ThermalVariables = mpFlowRule->GetThermalVariables();
  //     rValue=ThermalVariables.DeltaPlasticDissipation;
  //   }

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

  mElasticLeftCauchyGreen = identity_matrix<double> (3);

  mpFlowRule->InitializeMaterial( mpYieldCriterion, mpHardeningLaw, rMaterialProperties );

}

//************************************************************************************
//************************************************************************************


void HyperElasticPlastic3DLaw::InitializeSolutionStep( const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry, //this is just to give the array of nodes
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{

}

//************************************************************************************
//************************************************************************************


void HyperElasticPlastic3DLaw::FinalizeSolutionStep( const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry, //this is just to give the array of nodes
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{

}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

Matrix& HyperElasticPlastic3DLaw::DeformationGradient3D (Matrix & Matrix2D)
{
    //Takes a matrix 2x2 and transforms it to a 3x3 adding a 3rd row and a 3rd column with a 1 in the diagonal

    if (Matrix2D.size1() == 2 && Matrix2D.size2() == 2)
    {

        Matrix2D.resize( 3, 3, true);

        Matrix2D( 0 , 2 ) = 0.0;
        Matrix2D( 1 , 2 ) = 0.0;

        Matrix2D( 2 , 0 ) = 0.0;
        Matrix2D( 2 , 1 ) = 0.0;

        Matrix2D( 2 , 2 ) = 1.0;

    }
    else if(Matrix2D.size1() != 3 && Matrix2D.size2() != 3)
    {

        KRATOS_ERROR(std::invalid_argument,"Passed Matrix dimensions in DeformtationGradient3D not correct ","");

    }

    return Matrix2D;

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
  //Matrix& ConstitutiveMatrix         = rValues.GetConstitutiveMatrix();

  //2.-Green-Lagrange Strain:
  if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
    {
      TransformStrains (StrainVector, DeformationGradientF, StrainMeasure_Almansi, StrainMeasure_GreenLagrange);
    }

  //3.-Calculate Total PK2 stress
  if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {
      TransformStresses(StressVector, DeformationGradientF, DeterminantF, StressMeasure_Kirchhoff, StressMeasure_PK2);
    }

  //4.-Calculate PK2 constitutive tensor
  //if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
  //{
      //the tensor is calculated automatically in the initial configuration 
      //because the Options are checked before the Constitutive Matrix calculation

      /* Options.Is( ConstitutiveLaw::INITIAL_CONFIGURATION ) || 
	 Options.Is( ConstitutiveLaw::LAST_KNOWN_CONFIGURATION ) || 
	 Options.Is( ConstitutiveLaw::FINAL_CONFIGURATION ) */

  //}


   // std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
   // std::cout<<" Stress "<<StressVector<<std::endl;

}


//************************************************************************************
//************************************************************************************


void HyperElasticPlastic3DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
{
    this->CalculateMaterialResponsePK2 (rValues);

    Vector& StressVector = rValues.GetStressVector();
    const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
    const double& DeterminantF         = rValues.GetDeterminantF();

    TransformStresses(StressVector,DeformationGradientF,DeterminantF,StressMeasure_PK2,StressMeasure_PK1);
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

    const ProcessInfo&  CurProcessInfo    = rValues.GetProcessInfo();
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();
    const double&   DeterminantF          = rValues.GetDeterminantF();

    const GeometryType&  DomainGeometry   = rValues.GetElementGeometry ();
    const Vector&        ShapeFunctions   = rValues.GetShapeFunctionsValues ();

    Vector& StrainVector                  = rValues.GetStrainVector();

    double& DeterminantF0                 = rValues.GetDeterminantF0();

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.IdentityMatrix = identity_matrix<double> ( 3 );

    ElasticVariables.SetElementGeometry(DomainGeometry);
    ElasticVariables.SetShapeFunctionsValues(ShapeFunctions);

    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize(); //it has to be called at the start

    // Initialize variables from the process information
    ReturnMappingVariables.DeltaTime = CurProcessInfo[DELTA_TIME];
    
    if(CurProcessInfo[IMPLEX] == 1)	
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

    //2.-Determinant of the Total Deformation Gradient
    ElasticVariables.DeterminantF0 = DeterminantF0 * DeterminantF;

    //3.-Compute DeformationGradient_bar
    ElasticVariables.J_pow13 = pow(DeterminantF,1.0/3.0);

    Matrix DeformationGradientFbar = DeformationGradientF;

    DeformationGradientFbar = DeformationGradient3D(DeformationGradientFbar);

    DeformationGradientFbar /= ElasticVariables.J_pow13;


    //4.-Left Cauchy-Green tensor b_bar to the new configuration
    ElasticVariables.CauchyGreenMatrix = prod(mElasticLeftCauchyGreen,trans(DeformationGradientFbar));
    ElasticVariables.CauchyGreenMatrix = prod(DeformationGradientFbar,ElasticVariables.CauchyGreenMatrix);

    //5.-Calculate trace of Left Cauchy-Green tensor b_bar
    ElasticVariables.traceCG = 0;
    for( unsigned int i=0; i<3; i++)
    {
       ElasticVariables.traceCG += ElasticVariables.CauchyGreenMatrix( i , i );
    }

    ReturnMappingVariables.LameMu_bar = ElasticVariables.LameMu * ( ElasticVariables.traceCG / 3.0  );

    //4.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
    {
        // e= 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix,StrainVector);
	// correct b_bar to b
	double J_pow23 = pow(ElasticVariables.DeterminantF0,2.0/3.0);
	StrainVector /= (J_pow23*J_pow23);
    }

 
    //5.-Calculate Total Kirchhoff stress
    SplitStressVector.Isochoric = ZeroVector(voigtsize);
    Matrix IsochoricStressMatrix = ZeroMatrix(3);

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
      this->CalculateIsochoricStress( ElasticVariables, ReturnMappingVariables, StressMeasure_Kirchhoff, IsochoricStressMatrix, SplitStressVector.Isochoric );

    //OPTION 1:
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {

        SplitStressVector.Volumetric = ZeroVector(voigtsize);

        ElasticVariables.CauchyGreenMatrix = ElasticVariables.IdentityMatrix;

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
      
        if( ReturnMappingVariables.Options.Is(FlowRule::NOT_RETURN_MAPPING_COMPUTED) )
	  KRATOS_ERROR(std::logic_error, " ReturnMappingCall was not performed  ...error in the constitutive calculation...","");

        //initialize constitutive tensors
        ConstitutiveMatrix.clear();
        SplitConstitutiveMatrix.Isochoric  = ConstitutiveMatrix;
        SplitConstitutiveMatrix.Volumetric = ConstitutiveMatrix;
	SplitConstitutiveMatrix.Plastic    = ConstitutiveMatrix;

        ElasticVariables.CauchyGreenMatrix = ElasticVariables.IdentityMatrix;

	if ( Options.Is( ConstitutiveLaw::LAST_KNOWN_CONFIGURATION ) || Options.Is( ConstitutiveLaw::INITIAL_CONFIGURATION ) ){
	  

	  Matrix InverseDeformationGradientF ( 3, 3 );
	  double DetInvF=0;
	  MathUtils<double>::InvertMatrix( DeformationGradientF, InverseDeformationGradientF, DetInvF);


	  this->CalculateIsochoricConstitutiveMatrix  ( ElasticVariables, InverseDeformationGradientF, ReturnMappingVariables.TrialIsoStressMatrix, SplitConstitutiveMatrix.Isochoric );
	  
	  this->CalculateVolumetricConstitutiveMatrix ( ElasticVariables, InverseDeformationGradientF, SplitConstitutiveMatrix.Volumetric );

	  
	  if( ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION) ){
	    this->CalculatePlasticConstitutiveMatrix  ( ElasticVariables, InverseDeformationGradientF, ReturnMappingVariables, SplitConstitutiveMatrix.Plastic );
	  }

	}
	else{
	  
	  this->CalculateIsochoricConstitutiveMatrix  ( ElasticVariables, ReturnMappingVariables.TrialIsoStressMatrix, SplitConstitutiveMatrix.Isochoric );
	  
	  this->CalculateVolumetricConstitutiveMatrix ( ElasticVariables, SplitConstitutiveMatrix.Volumetric );
	  
	  if( ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION) )
	    this->CalculatePlasticConstitutiveMatrix  ( ElasticVariables, ReturnMappingVariables, SplitConstitutiveMatrix.Plastic );
	  
	}

	// std::cout<< " Isochoric Constitutive "<<SplitConstitutiveMatrix.Isochoric<<std::endl;
	// std::cout<< " Volumetric Constitutive "<<SplitConstitutiveMatrix.Volumetric<<std::endl;
	//if( ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION) )
	  //std::cout<< " Plastic Constitutive   "<<SplitConstitutiveMatrix.Plastic<<std::endl;
	

        //if( Options.Is(ConstitutiveLaw::TOTAL_TENSOR ) )
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
      mElasticLeftCauchyGreen += ( ElasticVariables.traceCG/3.0) * ElasticVariables.IdentityMatrix;
 
    }

}


//************************************************************************************
//************************************************************************************

void HyperElasticPlastic3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{

    this->CalculateMaterialResponseKirchhoff (rValues);

    Vector& StressVector                = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix          = rValues.GetConstitutiveMatrix();
    double& DeterminantF0               = rValues.GetDeterminantF0();
    const double& DeterminantF          = rValues.GetDeterminantF();

    double detF0 = DeterminantF0 * DeterminantF;

    //Set to cauchy Stress:
    StressVector       /= detF0;
    ConstitutiveMatrix /= detF0;

    //std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    //std::cout<<" Stress "<<StressVector<<std::endl;

}


//***********************************UPDATE*******************************************
//************************************************************************************

void HyperElasticPlastic3DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
{
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponsePK2 (rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    Vector& StressVector                   =rValues.GetStressVector();
    Matrix& DeformationGradientF0          =rValues.GetDeformationGradientF0();
    double& DeterminantF0                  =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& DeterminantF             =rValues.GetDeterminantF();

    //1.-Push-Forward to the updated configuration to be used as a reference in the next step
    TransformStresses(StressVector,DeformationGradientF,DeterminantF,StressMeasure_PK2,StressMeasure_Cauchy);  //Cauchy Stress

    //2.-Update Internal Variables
    DeformationGradientF0  =  prod(DeformationGradientF, DeformationGradientF0);
    DeterminantF0         *=  DeterminantF;
}

//************************************************************************************
//************************************************************************************


void HyperElasticPlastic3DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
{
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponsePK1 (rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    Vector& StressVector                   =rValues.GetStressVector();
    Matrix& DeformationGradientF0          =rValues.GetDeformationGradientF0();
    double& DeterminantF0                  =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& DeterminantF             =rValues.GetDeterminantF();

    //1.-Push-Forward to the updated configuration to be used as a reference in the next step
    TransformStresses(StressVector,DeformationGradientF,DeterminantF,StressMeasure_PK1,StressMeasure_Cauchy);  //increment of Cauchy Stress

    //2.-Update Internal Variables
    DeformationGradientF0  =  prod(DeformationGradientF, DeformationGradientF0);
    DeterminantF0         *=  DeterminantF;

}

//************************************************************************************
//************************************************************************************


void HyperElasticPlastic3DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
{

    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseKirchhoff (rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    Vector& StressVector                   =rValues.GetStressVector();
    Matrix& DeformationGradientF0          =rValues.GetDeformationGradientF0();
    double& DeterminantF0                  =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& DeterminantF             =rValues.GetDeterminantF();

    //1.-Push-Forward to the updated configuration to be used as a reference in the next step
    TransformStresses(StressVector,DeformationGradientF,DeterminantF,StressMeasure_Kirchhoff,StressMeasure_Cauchy);  //increment of Cauchy Stress

    //2.-Update Internal Variables
    DeformationGradientF0  =  prod(DeformationGradientF, DeformationGradientF0);
    DeterminantF0         *=  DeterminantF;

}


//************************************************************************************
//************************************************************************************

void HyperElasticPlastic3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseCauchy (rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    Matrix& DeformationGradientF0          =rValues.GetDeformationGradientF0();
    double& DeterminantF0                  =rValues.GetDeterminantF0();
    const Matrix& DeformationGradientF     =rValues.GetDeformationGradientF();
    const double& DeterminantF             =rValues.GetDeterminantF();

    //2.-Update Internal Variables
    DeformationGradientF0  =  prod(DeformationGradientF, DeformationGradientF0);
    DeterminantF0         *=  DeterminantF;
}

//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HyperElasticPlastic3DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = 0.5 * ( rRightCauchyGreen( 2, 2 ) - 1.00 );
    rStrainVector[3] = rRightCauchyGreen( 0, 1 ); // xy
    rStrainVector[4] = rRightCauchyGreen( 1, 2 ); // yz
    rStrainVector[5] = rRightCauchyGreen( 0, 2 ); // xz

}

//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HyperElasticPlastic3DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{

    // e= 0.5*(1-invbT*invb)

    //Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen ( 3, 3 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 2, 2 ) );
    rStrainVector[3] = - InverseLeftCauchyGreen( 0, 1 ); // xy
    rStrainVector[4] = - InverseLeftCauchyGreen( 1, 2 ); // yz
    rStrainVector[5] = - InverseLeftCauchyGreen( 0, 2 ); // xz
}


//******************************* COMPUTE ISOCHORIC STRESS  **************************
//************************************************************************************

void HyperElasticPlastic3DLaw::CalculateIsochoricStress( MaterialResponseVariables & rElasticVariables,
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
        rIsoStressMatrix  = (rElasticVariables.IdentityMatrix - (rElasticVariables.traceCG/3.0)*rElasticVariables.CauchyGreenMatrix );
        rIsoStressMatrix *= rElasticVariables.LameMu;

        //std::cout<<" PK2 "<<std::endl;
    }

    if(rStressMeasure == StressMeasure_Kirchhoff)
    {

        //rElasticVariables.CauchyGreenMatrix is LeftCauchyGreen

        //2.-Isochoric part of the Kirchhoff Stress Matrix
        rIsoStressMatrix  = (rElasticVariables.CauchyGreenMatrix - (rElasticVariables.traceCG/3.0)*rElasticVariables.IdentityMatrix );
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

//******************************* COMPUTE DOMAIN TEMPERATURE  ************************
//************************************************************************************


double &  HyperElasticPlastic3DLaw::CalculateDomainTemperature (const MaterialResponseVariables & rElasticVariables,
								double & rTemperature)
{
  
    //1.-Temperature from nodes
    // const GeomteryType& DomainGeometry = rElasticVariables.GetElementGeometry();
    // const Vector& ShapeFunctionsValues = rElasticVariables.GetShapeFunctionsValues();
    // const unsigned int number_of_nodes = DomainGeometry.size();
    
    // rTemperature=0;
       
    // for ( unsigned int j = 0; j < number_of_nodes; j++ )
    //   {
    // 	rTemperature += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(TEMPERATURE);
    //   }

    //2.-Temperature not included
    rTemperature = 0;

    return rTemperature;
}


//******************************* COMPUTE DOMAIN PRESSURE  ***************************
//************************************************************************************


double &  HyperElasticPlastic3DLaw::CalculateDomainPressure (const MaterialResponseVariables & rElasticVariables,
							double & rPressure)
{

    double BulkModulus = rElasticVariables.LameLambda + (2.0/3.0) * rElasticVariables.LameMu;

    double auxiliar = 0;
    
    //(J²-1)/2
    //auxiliar =  0.5*(rElasticVariables.DeterminantF0*rElasticVariables.DeterminantF0-1);
    
    //(ln(J))
    auxiliar =  std::log(rElasticVariables.DeterminantF0);

    rPressure = BulkModulus *  auxiliar;

    return rPressure;
}

//************************* COMPUTE DOMAIN PRESSURE FACTORS***************************
//************************************************************************************

Vector&  HyperElasticPlastic3DLaw::CalculateDomainPressureFactors (const MaterialResponseVariables & rElasticVariables,
							      Vector & rFactors)
							      
{

    double BulkModulus = rElasticVariables.LameLambda + (2.0/3.0) * rElasticVariables.LameMu;

    if(rFactors.size()!=3) rFactors.resize(3);

    //(J²-1)/2
    // rFactors[0] =  rElasticVariables.DeterminantF0*rElasticVariables.DeterminantF0;
    // rFactors[1] =  (rElasticVariables.DeterminantF0*rElasticVariables.DeterminantF0-1);
    // rFactors[2] =  BulkModulus*rElasticVariables.DeterminantF0;

    //(ln(J))
    rFactors[0] =  1.0;
    rFactors[1] =  (2.0*std::log(rElasticVariables.DeterminantF0));
    rFactors[2] =  BulkModulus;

    return rFactors;
}



//******************************* COMPUTE VOLUMETRIC STRESS  *************************
//************************************************************************************

void HyperElasticPlastic3DLaw::CalculateVolumetricStress(const MaterialResponseVariables & rElasticVariables,
       Vector& rVolStressVector )
{

    //1.- Declaration
    Matrix VolStressMatrix ( 3 , 3 );

    double Pressure = 0;

    Pressure = this->CalculateDomainPressure (rElasticVariables, Pressure);

    //2.- Volumetric part of the Kirchhoff StressMatrix from nodal pressures
    VolStressMatrix = rElasticVariables.DeterminantF0 * Pressure * rElasticVariables.CauchyGreenMatrix;


    rVolStressVector = MathUtils<double>::StressTensorToVector(VolStressMatrix,rVolStressVector.size());

}

//***********************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX************************
//************************************************************************************

void HyperElasticPlastic3DLaw::CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
        const Matrix & rIsoStressMatrix,
        Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = IsochoricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, rIsoStressMatrix,
                                          this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

//***********************COMPUTE VOLUMETRIC CONSTITUTIVE MATRIX***********************
//************************************************************************************


void HyperElasticPlastic3DLaw::CalculateVolumetricConstitutiveMatrix ( const MaterialResponseVariables & rElasticVariables,
								  Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    Vector Factors = ZeroVector(3);
    Factors = this->CalculateDomainPressureFactors( rElasticVariables, Factors );


    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, Factors,
                                          this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


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
    
//********************CONSTITUTIVE MATRIX ISOCHORIC COMPONENT*************************
//************************************************************************************


double& HyperElasticPlastic3DLaw::IsochoricConstitutiveComponent(double & rCabcd,
        const MaterialResponseVariables & rElasticVariables,
        const Matrix & rIsoStressMatrix,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{

    //Isochoric part of the hyperelastic constitutive tensor component

    rCabcd  = (1.0/3.0)*(rElasticVariables.CauchyGreenMatrix(a,b)*rElasticVariables.CauchyGreenMatrix(c,d));
    rCabcd -= (0.5*(rElasticVariables.CauchyGreenMatrix(a,c)*rElasticVariables.CauchyGreenMatrix(b,d)+rElasticVariables.CauchyGreenMatrix(a,d)*rElasticVariables.CauchyGreenMatrix(b,c)));
    
    rCabcd *= rElasticVariables.traceCG * rElasticVariables.LameMu;
    
    rCabcd += (rElasticVariables.CauchyGreenMatrix(c,d)*rIsoStressMatrix(a,b) + rIsoStressMatrix(c,d)*rElasticVariables.CauchyGreenMatrix(a,b));
    rCabcd *= (-2.0/3.0);

    return rCabcd;
}


//********************CONSTITUTIVE MATRIX VOLUMETRIC COMPONENT************************
//************************************************************************************


double& HyperElasticPlastic3DLaw::VolumetricConstitutiveComponent(double & rCabcd,
        const MaterialResponseVariables & rElasticVariables,
        const Vector & rFactors,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{

    //Volumetric part of the hyperelastic constitutive tensor component: (J²-1)/2  -  (ln(J)/J)

    //1.Volumetric Elastic constitutive tensor component:
    rCabcd  = rFactors[0]*(rElasticVariables.CauchyGreenMatrix(a,b)*rElasticVariables.CauchyGreenMatrix(c,d));
    rCabcd -= rFactors[1]*(0.5*(rElasticVariables.CauchyGreenMatrix(a,c)*rElasticVariables.CauchyGreenMatrix(b,d)+rElasticVariables.CauchyGreenMatrix(a,d)*rElasticVariables.CauchyGreenMatrix(b,c)));
    rCabcd *= rFactors[2];
    

    return rCabcd;
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



//***********************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX PULL-BACK**************
//************************************************************************************

void HyperElasticPlastic3DLaw::CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
								     const Matrix & rInverseDeformationGradientF,
								     const Matrix & rIsoStressMatrix,
								     Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();


    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
	  rConstitutiveMatrix( i, j ) = IsochoricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, rInverseDeformationGradientF, rIsoStressMatrix,
                                          this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

//***********************COMPUTE VOLUMETRIC CONSTITUTIVE MATRIX PULL-BACK*************
//************************************************************************************


void HyperElasticPlastic3DLaw::CalculateVolumetricConstitutiveMatrix ( const MaterialResponseVariables & rElasticVariables,
								       const Matrix & rInverseDeformationGradientF,
								       Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    Vector Factors = ZeroVector(3);
    Factors = this->CalculateDomainPressureFactors( rElasticVariables, Factors );


    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
	  rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, rInverseDeformationGradientF, Factors,
                                          this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}


//***********************COMPUTE PLASTIC CONSTITUTIVE MATRIX PULL-BACK****************
//************************************************************************************

void HyperElasticPlastic3DLaw::CalculatePlasticConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
								   const Matrix & rInverseDeformationGradientF,
								   FlowRule::RadialReturnVariables & rReturnMappingVariables,
								   Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    Matrix IsoStressMatrix = rReturnMappingVariables.TrialIsoStressMatrix;

    FlowRule::PlasticFactors ScalingFactors;
    mpFlowRule->CalculateScalingFactors( rReturnMappingVariables, ScalingFactors );


    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
	  rConstitutiveMatrix( i, j ) = PlasticConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, rInverseDeformationGradientF, IsoStressMatrix, ScalingFactors,
									   this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

//********************CONSTITUTIVE MATRIX ISOCHORIC COMPONENT PULL-BACK***************
//************************************************************************************


double& HyperElasticPlastic3DLaw::IsochoricConstitutiveComponent(double & rCabcd,
								 const MaterialResponseVariables & rElasticVariables,
								 const Matrix & rInverseDeformationGradientF,
								 const Matrix & rIsoStressMatrix,
								 const unsigned int& a, const unsigned int& b,
								 const unsigned int& c, const unsigned int& d)
{
    rCabcd = 0;
    double Cijkl=0;

    unsigned int dimension = rInverseDeformationGradientF.size1();

    //Cabcd
    for(unsigned int j=0; j<dimension; j++)
    {
        for(unsigned int l=0; l<dimension; l++)
        {
            for(unsigned int k=0; k<dimension; k++)
            {
                for(unsigned int i=0; i<dimension; i++)
                {
                    //Cijkl
		  rCabcd +=rInverseDeformationGradientF(a,i)*rInverseDeformationGradientF(b,j)*rInverseDeformationGradientF(c,k)*rInverseDeformationGradientF(d,l)*IsochoricConstitutiveComponent(Cijkl,rElasticVariables,rIsoStressMatrix,i,j,k,l);
                }
            }
        }
    }

    return rCabcd;
}


//********************CONSTITUTIVE MATRIX VOLUMETRIC COMPONENT PULL-BACK**************
//************************************************************************************


double& HyperElasticPlastic3DLaw::VolumetricConstitutiveComponent(double & rCabcd,
								  const MaterialResponseVariables & rElasticVariables,
								  const Matrix & rInverseDeformationGradientF,
								  const Vector & rFactors,
								  const unsigned int& a, const unsigned int& b,
								  const unsigned int& c, const unsigned int& d)
{
    rCabcd = 0;
    double Cijkl=0;

    unsigned int dimension = rInverseDeformationGradientF.size1();

    //Cabcd
    for(unsigned int j=0; j<dimension; j++)
    {
        for(unsigned int l=0; l<dimension; l++)
        {
            for(unsigned int k=0; k<dimension; k++)
            {
                for(unsigned int i=0; i<dimension; i++)
                {
                    //Cijkl
		  rCabcd +=rInverseDeformationGradientF(a,i)*rInverseDeformationGradientF(b,j)*rInverseDeformationGradientF(c,k)*rInverseDeformationGradientF(d,l)*VolumetricConstitutiveComponent(Cijkl,rElasticVariables,rFactors,i,j,k,l);
                }
            }
        }
    }

    return rCabcd;
}


//********************CONSTITUTIVE MATRIX ISOCHORIC COMPONENT PULL-BACK***************
//************************************************************************************


double& HyperElasticPlastic3DLaw::PlasticConstitutiveComponent(double & rCabcd,
							       const MaterialResponseVariables & rElasticVariables,
							       const Matrix & rInverseDeformationGradientF,
							       const Matrix & rIsoStressMatrix,
							       const FlowRule::PlasticFactors & rScalingFactors,					
							       const unsigned int& a, const unsigned int& b,
							       const unsigned int& c, const unsigned int& d)
{
     rCabcd = 0;
    double Cijkl=0;

    unsigned int dimension = rInverseDeformationGradientF.size1();

    //Cabcd
    for(unsigned int j=0; j<dimension; j++)
    {
        for(unsigned int l=0; l<dimension; l++)
        {
            for(unsigned int k=0; k<dimension; k++)
            {
                for(unsigned int i=0; i<dimension; i++)
                {
                    //Cijkl
		  rCabcd +=rInverseDeformationGradientF(a,i)*rInverseDeformationGradientF(b,j)*rInverseDeformationGradientF(c,k)*rInverseDeformationGradientF(d,l)*PlasticConstitutiveComponent(Cijkl,rElasticVariables,rIsoStressMatrix,rScalingFactors,i,j,k,l);
                }
            }
        }
    }

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
                             const ProcessInfo& rCurrentProcessInfo)
{

    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
        KRATOS_ERROR(std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ","");

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
        KRATOS_ERROR(std::invalid_argument,"POISSON_RATIO has Key zero invalid value ","");


    if(DENSITY.Key() == 0 || rMaterialProperties[DENSITY]<0.00)
        KRATOS_ERROR(std::invalid_argument,"DENSITY has Key zero or invalid value ","");


    return 0;

}

} // Namespace Kratos
