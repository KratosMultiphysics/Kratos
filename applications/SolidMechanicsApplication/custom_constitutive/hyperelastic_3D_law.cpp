//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include <cmath>

// Project includes
#include "custom_constitutive/hyperelastic_3D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElastic3DLaw::HyperElastic3DLaw()
    : ConstitutiveLaw()
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElastic3DLaw::HyperElastic3DLaw(const HyperElastic3DLaw& rOther)
    : ConstitutiveLaw(rOther)
    ,mInverseDeformationGradientF0(rOther.mInverseDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
    ,mStrainEnergy(rOther.mStrainEnergy)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElastic3DLaw::Clone() const
{
    return Kratos::make_shared<HyperElastic3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElastic3DLaw::~HyperElastic3DLaw()
{
}


//*******************************OPERATIONS FROM BASE CLASS***************************
//************************************************************************************

//***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
//************************************************************************************

bool HyperElastic3DLaw::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool HyperElastic3DLaw::Has( const Variable<Vector>& rThisVariable )
{
    return false;
}

bool HyperElastic3DLaw::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}


//******************CALCULATE VALUE: DOUBLE - VECTOR - MATRIX*************************
//************************************************************************************

double& HyperElastic3DLaw::CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue )
{

  return (this->GetValue(rThisVariable,rValue ));

}


//***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************

double& HyperElastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
  if (rThisVariable == STRAIN_ENERGY)
  {
    rValue = mStrainEnergy;
  }
  else{
    rValue = 0;
  }


    return( rValue );
}

Vector& HyperElastic3DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return( rValue );
}

Matrix& HyperElastic3DLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return( rValue );
}


//***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************


void HyperElastic3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                  const ProcessInfo& rCurrentProcessInfo )
{

  if (rThisVariable == DETERMINANT_F)
    {
      mDeterminantF0 = rValue;
    }
}

void HyperElastic3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                  const ProcessInfo& rCurrentProcessInfo )
{

}

void HyperElastic3DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                  const ProcessInfo& rCurrentProcessInfo )
{

}



//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************


void HyperElastic3DLaw::InitializeMaterial( const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues )
{
  mDeterminantF0                = 1;
  mInverseDeformationGradientF0 = identity_matrix<double> (3);
  mStrainEnergy                 = 0;

}

//************************************************************************************
//************************************************************************************


void HyperElastic3DLaw::InitializeSolutionStep( const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry, //this is just to give the array of nodes
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{

}

//************************************************************************************
//************************************************************************************


void HyperElastic3DLaw::FinalizeSolutionStep( const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry, //this is just to give the array of nodes
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{

}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************



//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************


void  HyperElastic3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
{

  //-----------------------------//
  //a.-Check if the constitutive parameters are passed correctly to the law calculation
  CheckParameters(rValues);
  mStrainEnergy = 0.0; //When it is not calculated, a zero will be returned

  //b.- Get Values to compute the constitutive law:
  Flags &Options=rValues.GetOptions();

  const Properties& MaterialProperties  = rValues.GetMaterialProperties();
  const Matrix& DeformationGradientF    = rValues.GetDeformationGradientF();
  const double& DeterminantF            = rValues.GetDeterminantF();

  Vector& StrainVector                  = rValues.GetStrainVector();
  Vector& StressVector                  = rValues.GetStressVector();
  Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

  //-----------------------------//

  //0.- Initialize parameters
  MaterialResponseVariables ElasticVariables;
  ElasticVariables.Identity = identity_matrix<double> ( 3 );

  //1.- Lame constants
  const double& YoungModulus        = MaterialProperties[YOUNG_MODULUS];
  const double& PoissonCoefficient  = MaterialProperties[POISSON_RATIO];

  ElasticVariables.LameLambda       = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
  ElasticVariables.LameMu           =  YoungModulus/(2*(1+PoissonCoefficient));

  //2.- Thermal constants
  if( MaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT) )
    ElasticVariables.ThermalExpansionCoefficient = MaterialProperties[THERMAL_EXPANSION_COEFFICIENT];
  else
    ElasticVariables.ThermalExpansionCoefficient = 0;

  if( MaterialProperties.Has(REFERENCE_TEMPERATURE) )
    ElasticVariables.ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
  else
    ElasticVariables.ReferenceTemperature = 0;


  //3.-DeformationGradient Tensor 3D
  ElasticVariables.DeformationGradientF = DeformationGradientF;
  ElasticVariables.DeformationGradientF = Transform2DTo3D( ElasticVariables.DeformationGradientF );

  //4.-Determinant of the Total Deformation Gradient
  ElasticVariables.DeterminantF = DeterminantF;

  //5.-Right Cauchy Green tensor C
  Matrix RightCauchyGreen = prod(trans( ElasticVariables.DeformationGradientF),  ElasticVariables.DeformationGradientF);

  //6.-Inverse of the Right Cauchy-Green tensor C: (stored in the CauchyGreenMatrix)
  ElasticVariables.traceCG = 0;
  ElasticVariables.CauchyGreenMatrix.resize(3,3,false);
  MathUtils<double>::InvertMatrix( RightCauchyGreen, ElasticVariables.CauchyGreenMatrix, ElasticVariables.traceCG);

  //7.-Green-Lagrange Strain:
  if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
      this->CalculateGreenLagrangeStrain(RightCauchyGreen, StrainVector);
    }

  //8.-Calculate Total PK2 stress
  if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {

      this->CalculateStress( ElasticVariables, StressMeasure_PK2, StressVector );
    }

  //9.-Calculate Constitutive Matrix related to Total PK2 stress
  if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
      this->CalculateConstitutiveMatrix ( ElasticVariables, ConstitutiveMatrix );
    }

  if( Options.Is( ConstitutiveLaw::COMPUTE_STRAIN_ENERGY ) )
    {

      double ln_J = std::log(ElasticVariables.DeterminantF);
      double trace_C = 0.0;

      for (unsigned int i = 0; i<RightCauchyGreen.size1();i++)
	{
	  trace_C += RightCauchyGreen(i,i);
	}

      mStrainEnergy =  0.5*ElasticVariables.LameLambda*ln_J*ln_J - ElasticVariables.LameMu*ln_J + 0.5*ElasticVariables.LameMu*(trace_C-3); //see Belytschko page 239


    }

  // std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
  // std::cout<<" Stress "<<StressVector<<std::endl;

  //-----------------------------//

}


//************************************************************************************
//************************************************************************************


void HyperElastic3DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
{
    this->CalculateMaterialResponsePK2 (rValues);

    Vector& StressVector               = rValues.GetStressVector();
    const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
    const double& DeterminantF         = rValues.GetDeterminantF();

    TransformStresses(StressVector,DeformationGradientF,DeterminantF,StressMeasure_PK2,StressMeasure_PK1);
}

//************************************************************************************
//************************************************************************************


void HyperElastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();
    const double&   DeterminantF          = rValues.GetDeterminantF();

    Vector& StrainVector                  = rValues.GetStrainVector();
    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.Identity = identity_matrix<double> ( 3 );

    //1.- Lame constants
    const double& YoungModulus       = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];

    ElasticVariables.LameLambda      = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    ElasticVariables.LameMu          =  YoungModulus/(2*(1+PoissonCoefficient));

    //2.- Thermal constants
    if( MaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT) )
      ElasticVariables.ThermalExpansionCoefficient = MaterialProperties[THERMAL_EXPANSION_COEFFICIENT];
    else
      ElasticVariables.ThermalExpansionCoefficient = 0;

    if( MaterialProperties.Has(REFERENCE_TEMPERATURE) )
      ElasticVariables.ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
    else
      ElasticVariables.ReferenceTemperature = 0;

    //3.-DeformationGradient Tensor 3D
    ElasticVariables.DeformationGradientF = DeformationGradientF;
    ElasticVariables.DeformationGradientF = Transform2DTo3D( ElasticVariables.DeformationGradientF );

    //4.-Determinant of the Total Deformation Gradient
    ElasticVariables.DeterminantF         = DeterminantF;

    //5.-Left Cauchy Green tensor b: (stored in the CauchyGreenMatrix)
    ElasticVariables.CauchyGreenMatrix.resize(3,3,false);
    noalias(ElasticVariables.CauchyGreenMatrix) = prod(ElasticVariables.DeformationGradientF,trans(ElasticVariables.DeformationGradientF));

    for( unsigned int i=0; i<3; i++)
    {
       ElasticVariables.traceCG += ElasticVariables.CauchyGreenMatrix( i , i );
    }

    //6.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
        // e= 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix,StrainVector);
    }


    //7.-Calculate Total kirchhoff stress
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {
        this->CalculateStress( ElasticVariables, StressMeasure_Kirchhoff, StressVector );

    }

    //8.-Calculate Constitutive Matrix related to Total Kirchhoff stress
    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        ElasticVariables.CauchyGreenMatrix = ElasticVariables.Identity;
        this->CalculateConstitutiveMatrix ( ElasticVariables, ConstitutiveMatrix );
    }


    // std::cout<<" StrainVector "<<StrainVector<<std::endl;
    // std::cout<<" StressVector "<<StressVector<<std::endl;
    // std::cout<<" ConstitutiveMatrix "<<ConstitutiveMatrix<<std::endl;


    //-----------------------------//

}


//************************************************************************************
//************************************************************************************

void HyperElastic3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{

    this->CalculateMaterialResponseKirchhoff (rValues);

    const double& DeterminantF          = rValues.GetDeterminantF();
    Vector& StressVector                = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix          = rValues.GetConstitutiveMatrix();

     //Set to cauchy Stress:
    StressVector       /= DeterminantF;
    ConstitutiveMatrix /= DeterminantF;

}


//***********************************UPDATE*******************************************
//************************************************************************************

void HyperElastic3DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
{

    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponsePK2 (rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    UpdateInternalVariables( rValues );
}

//************************************************************************************
//************************************************************************************


void HyperElastic3DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
{

    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponsePK1 (rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    UpdateInternalVariables( rValues );
}

//************************************************************************************
//************************************************************************************


void HyperElastic3DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
{

    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseKirchhoff (rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    UpdateInternalVariables( rValues );
}


//************************************************************************************
//************************************************************************************

void HyperElastic3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseCauchy (rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    UpdateInternalVariables( rValues );
}


//************************************************************************************
//************************************************************************************

void HyperElastic3DLaw::UpdateInternalVariables(Parameters& rValues)
{
    const Matrix& DeformationGradientF    = rValues.GetDeformationGradientF();
    const double& DeterminantF            = rValues.GetDeterminantF();

    Matrix DeformationGradientF0          = DeformationGradientF;
    DeformationGradientF0 = Transform2DTo3D(DeformationGradientF0);
    MathUtils<double>::InvertMatrix( DeformationGradientF0, this->mInverseDeformationGradientF0, mDeterminantF0);
    mDeterminantF0 = DeterminantF; //special treatment of the determinant
}


//***********************COMPUTE TOTAL STRAIN VECTOR**********************************
//************************************************************************************

void HyperElastic3DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{

    //E= 0.5*(FT*F-1) or E = 0.5*(C-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = 0.5 * ( rRightCauchyGreen( 2, 2 ) - 1.00 );
    rStrainVector[3] = rRightCauchyGreen( 0, 1 ); // xy
    rStrainVector[4] = rRightCauchyGreen( 1, 2 ); // yz
    rStrainVector[5] = rRightCauchyGreen( 0, 2 ); // xz

    // Matrix StrainMatrix(3,3);
    // noalias(StrainMatrix) = ZeroMatrix(3,3);
    // CalculateAlmansiStrain( rRightCauchyGreen, rStrainMatrix );
    // rStrainVector = MathUtils<double>::StrainTensorToVector( StrainMatrix, rStrainVector.size() );

}


//***********************COMPUTE TOTAL STRAIN MATRIX**********************************
//************************************************************************************

// void HyperElastic3DLaw::CalculateAlmansiStrain( const Matrix & rRightCauchyGreen,
// 						Matrix& rStrainMatrix )
// {

//     //E= 0.5*(FT*F-1)
//     Matrix Identity = identity_matrix<double> ( 3 );
//     rStrainMatrix = 0.5 * prod( rRightCauchyGreen - Identity);

// }

//***********************COMPUTE TOTAL STRAIN VECTOR**********************************
//************************************************************************************

void HyperElastic3DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
						Vector& rStrainVector )
{

    // e = 0.5*(1-invFT*invF) or e = 0.5*(1-inv(b))

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


    // Matrix StrainMatrix(3,3);
    // noalias(StrainMatrix) = ZeroMatrix(3,3);
    // CalculateAlmansiStrain( rLeftCauchyGreen, rStrainMatrix );
    // rStrainVector = MathUtils<double>::StrainTensorToVector( StrainMatrix, rStrainVector.size() );

}


//***********************COMPUTE TOTAL STRAIN MATRIX**********************************
//************************************************************************************

// void HyperElastic3DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
// 						Matrix& rStrainMatrix )
// {

//     // e = 0.5*(1-invFT*invF) or e = 0.5*(1-inv(b))

//     //Calculating the inverse of the jacobian
//     Matrix InverseLeftCauchyGreen ( 3, 3 );
//     double det_b=0;
//     MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

//     Matrix Identity = identity_matrix<double> ( 3 );
//     rStrainMatrix = 0.5 * prod( Identity - InverseLeftCauchyGreen);

// }


//******************************* COMPUTE DOMAIN TEMPERATURE  ************************
//************************************************************************************


double &  HyperElastic3DLaw::CalculateDomainTemperature (const MaterialResponseVariables & rElasticVariables,
								double & rTemperature)
{

    //1.-Temperature from nodes
    const GeometryType& DomainGeometry = rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues = rElasticVariables.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();

    rTemperature=0;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
      {
	if( DomainGeometry[j].SolutionStepsDataHas(TEMPERATURE) )
	  rTemperature += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(TEMPERATURE);
      }

    //2.-Temperature not included
    //rTemperature = 0;

    return rTemperature;
}


//***************************** COMPUTE VOLUMETRIC FACTOR ****************************
//************************************************************************************


double & HyperElastic3DLaw::CalculateVolumetricFactor (const MaterialResponseVariables & rElasticVariables,
						       double & rFactor)
{
  //(J²-1)/2
  //rFactor = 0.5*(rElasticVariables.DeterminantF*rElasticVariables.DeterminantF-1);

  //(ln(J))
  rFactor = std::log(rElasticVariables.DeterminantF);

  return rFactor;

}


//***************************** COMPUTE VOLUMETRIC PRESSURE  *************************
//************************************************************************************


double & HyperElastic3DLaw::CalculateVolumetricPressure (const MaterialResponseVariables & rElasticVariables,
							 double & rPressure)
{

    double BulkModulus = rElasticVariables.LameLambda + (2.0/3.0) * rElasticVariables.LameMu;

    //Mechanical volumetric factor:
    double Factor = 0;
    Factor = this->CalculateVolumetricFactor( rElasticVariables, Factor );


    //Thermal volumetric factor:
    double DeltaTemperature     = 0;
    double CurrentTemperature   = 0;

    CurrentTemperature = this->CalculateDomainTemperature(rElasticVariables, CurrentTemperature);
    DeltaTemperature   = CurrentTemperature - rElasticVariables.ReferenceTemperature;

    Factor            += 3.0 * rElasticVariables.ThermalExpansionCoefficient * ( (1.0 - std::log(rElasticVariables.DeterminantF)) / (rElasticVariables.DeterminantF) ) * DeltaTemperature;


    rPressure = BulkModulus * Factor;

    return rPressure;
}


//************************* COMPUTE VOLUMETRIC PRESSURE FACTORS***********************
//************************************************************************************

Vector&  HyperElastic3DLaw::CalculateVolumetricPressureFactors (const MaterialResponseVariables & rElasticVariables,
							       Vector & rFactors)

{

    double BulkModulus = rElasticVariables.LameLambda + (2.0/3.0) * rElasticVariables.LameMu;

    if(rFactors.size()!=3) rFactors.resize(3);

    //(J²-1)/2
    // rFactors[0] =  rElasticVariables.DeterminantF*rElasticVariables.DeterminantF;
    // rFactors[1] =  (rElasticVariables.DeterminantF*rElasticVariables.DeterminantF-1);
    // rFactors[2] =  BulkModulus*rElasticVariables.DeterminantF;

    //(ln(J))
    rFactors[0] =  1.0;
    rFactors[1] =  (2.0*std::log(rElasticVariables.DeterminantF));
    rFactors[2] =  BulkModulus;

    return rFactors;
}

//******************************* COMPUTE TOTAL STRESS  ******************************
//************************************************************************************

void HyperElastic3DLaw::CalculateStress( const MaterialResponseVariables & rElasticVariables,
					 StressMeasure rStressMeasure,
					 Vector& rStressVector )
{

    Matrix StressMatrix( 3, 3 );

    //1.- Temporary and selected law

    double Factor = 0;
    Factor = this->CalculateVolumetricFactor( rElasticVariables, Factor );

    //(J²-1)/2
    //double Factor = 0.5*(rElasticVariables.DeterminantF*rElasticVariables.DeterminantF-1);

    //(ln(J))
    //double Factor = (std::log(rElasticVariables.DeterminantF));


    if(rStressMeasure == StressMeasure_PK2)  // the description corresponds to the neohookean material in Belytschko nonlinear finite elements, pag 239
    {
        //rElasticVariables.CauchyGreenMatrix is InverseRightCauchyGreen C^-1

        //2.-2nd Piola Kirchhoff Stress Matrix
        StressMatrix  = rElasticVariables.LameLambda * Factor * rElasticVariables.CauchyGreenMatrix;
        StressMatrix += rElasticVariables.LameMu * ( rElasticVariables.Identity - rElasticVariables.CauchyGreenMatrix );

    }

    if(rStressMeasure == StressMeasure_Kirchhoff) // the description corresponds to the neohookean material in Belytschko nonlinear finite elements, pag 239
    {
        //rElasticVariables.CauchyGreenMatrix is LeftCauchyGreen B

        //2.-Kirchhoff Stress Matrix
        StressMatrix  = rElasticVariables.LameLambda * Factor * rElasticVariables.Identity;

        StressMatrix += rElasticVariables.LameMu * ( rElasticVariables.CauchyGreenMatrix - rElasticVariables.Identity );

   }

    rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );

}



//******************************* COMPUTE ISOCHORIC STRESS  **************************
//************************************************************************************
void HyperElastic3DLaw::CalculateIsochoricStress( const MaterialResponseVariables & rElasticVariables,
						  StressMeasure rStressMeasure,
						  Vector& rIsoStressVector )
{

    //1.-Identity build
    Matrix IsoStressMatrix ( 3, 3 );

    //note.- rElasticVariables.traceCG is "traceCG"

    if(rStressMeasure == StressMeasure_PK2)
    {

        //rElasticVariables.CauchyGreenMatrix is InverseRightCauchyGreen

        //2.-Incompressible part of the 2nd Piola Kirchhoff Stress Matrix
        IsoStressMatrix  = (rElasticVariables.Identity - (rElasticVariables.traceCG/3.0)*rElasticVariables.CauchyGreenMatrix );
        IsoStressMatrix *= rElasticVariables.LameMu*pow(rElasticVariables.DeterminantF,(-2.0/3.0));

        //std::cout<<" PK2 "<<std::endl;
    }

    if(rStressMeasure == StressMeasure_Kirchhoff)
    {

        //rElasticVariables.CauchyGreenMatrix is LeftCauchyGreen

        //2.-Incompressible part of the Kirchhoff Stress Matrix
        IsoStressMatrix  = (rElasticVariables.CauchyGreenMatrix - (rElasticVariables.traceCG/3.0)*rElasticVariables.Identity );
        IsoStressMatrix *= rElasticVariables.LameMu*pow(rElasticVariables.DeterminantF,(-2.0/3.0));

        //std::cout<<" Kirchooff "<<std::endl;

    }

    rIsoStressVector = MathUtils<double>::StressTensorToVector(IsoStressMatrix,rIsoStressVector.size());

}


//******************************* COMPUTE VOLUMETRIC STRESS  *************************
//************************************************************************************

void HyperElastic3DLaw::CalculateVolumetricStress(const MaterialResponseVariables & rElasticVariables,
						  Vector& rVolStressVector )
{

    //1.- Declaration
    Matrix VolStressMatrix ( 3 , 3 );

    double Pressure = 0;

    Pressure = this->CalculateVolumetricPressure (rElasticVariables, Pressure);

    //2.- Volumetric part of the Kirchhoff StressMatrix from nodal pressures
    VolStressMatrix = rElasticVariables.DeterminantF * Pressure * rElasticVariables.CauchyGreenMatrix;


    rVolStressVector = MathUtils<double>::StressTensorToVector(VolStressMatrix,rVolStressVector.size());

}


//***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
//************************************************************************************

void HyperElastic3DLaw::CalculateConstitutiveMatrix ( const MaterialResponseVariables& rElasticVariables,
						      Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = ConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables,
                                          this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

//***********************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX************************
//************************************************************************************

void HyperElastic3DLaw::CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
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


void HyperElastic3DLaw::CalculateVolumetricConstitutiveMatrix ( const MaterialResponseVariables & rElasticVariables,
								Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    Vector Factors(3);
    noalias(Factors) = ZeroVector(3);
    Factors = this->CalculateVolumetricPressureFactors( rElasticVariables, Factors );


    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, Factors,
                                          this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

//***********************CONSTITUTIVE TENSOR COMPONENTS*******************************
//************************************************************************************


double& HyperElastic3DLaw::ConstitutiveComponent(double & rCabcd,
						 const MaterialResponseVariables& rElasticVariables,
						 const unsigned int& a, const unsigned int& b,
						 const unsigned int& c, const unsigned int& d)
{

    //1.- Temporary and selected law
    Vector Factors(3);
    noalias(Factors) = ZeroVector(3);
    Factors = this->CalculateVolumetricPressureFactors( rElasticVariables, Factors );

    double auxiliar1 = Factors[0];
    double auxiliar2 = Factors[1];

    //(J²-1)/2
    //double auxiliar1 =  rElasticVariables.DeterminantF*rElasticVariables.DeterminantF;
    //double auxiliar2 =  (rElasticVariables.DeterminantF*rElasticVariables.DeterminantF-1);

    //(ln(J))
    //double auxiliar1 =  1.0;
    //double auxiliar2 =  (2.0*std::log(rElasticVariables.DeterminantF));


    //1.Elastic constitutive tensor component:
    rCabcd =(rElasticVariables.LameLambda*auxiliar1*rElasticVariables.CauchyGreenMatrix(a,b)*rElasticVariables.CauchyGreenMatrix(c,d));
    rCabcd+=((2*rElasticVariables.LameMu-rElasticVariables.LameLambda*auxiliar2)*0.5*(rElasticVariables.CauchyGreenMatrix(a,c)*rElasticVariables.CauchyGreenMatrix(b,d)+rElasticVariables.CauchyGreenMatrix(a,d)*rElasticVariables.CauchyGreenMatrix(b,c)));

    return rCabcd;
}


//********************CONSTITUTIVE MATRIX ISOCHORIC COMPONENT*************************
//************************************************************************************


double& HyperElastic3DLaw::IsochoricConstitutiveComponent(double & rCabcd,
							  const MaterialResponseVariables & rElasticVariables,
							  const Matrix & rIsoStressMatrix,
							  const unsigned int& a, const unsigned int& b,
							  const unsigned int& c, const unsigned int& d)
{

    //Isochoric part of the hyperelastic constitutive tensor component

    //note.- rElasticVariables.traceCG is "traceCG_bar"

    rCabcd  = (1.0/3.0)*(rElasticVariables.CauchyGreenMatrix(a,b)*rElasticVariables.CauchyGreenMatrix(c,d));

    rCabcd -= (0.5*(rElasticVariables.CauchyGreenMatrix(a,c)*rElasticVariables.CauchyGreenMatrix(b,d)+rElasticVariables.CauchyGreenMatrix(a,d)*rElasticVariables.CauchyGreenMatrix(b,c)));

    rCabcd *= rElasticVariables.traceCG * rElasticVariables.LameMu;

    rCabcd += (rElasticVariables.CauchyGreenMatrix(c,d)*rIsoStressMatrix(a,b) + rIsoStressMatrix(c,d)*rElasticVariables.CauchyGreenMatrix(a,b));

    rCabcd *= (-2.0/3.0);


    return rCabcd;
}


//********************CONSTITUTIVE MATRIX VOLUMETRIC COMPONENT************************
//************************************************************************************


double& HyperElastic3DLaw::VolumetricConstitutiveComponent(double & rCabcd,
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

//*************************CONSTITUTIVE LAW PARTICULAR UTILITIES**********************
//************************************************************************************

/**
 * Takes a matrix 2x2 and transforms it to a 3x3 adding a 3rd row and a 3rd column with a 1 in the diagonal
 * if the matrix passed is 3D is does nothing
 * if the matrix passed is bigger or smaller throws an error
 * @param rMatrix : usually the DeformationGradientF
 */
Matrix& HyperElastic3DLaw::Transform2DTo3D (Matrix& rMatrix)
{


    if (rMatrix.size1() == 2 && rMatrix.size2() == 2)
    {

        rMatrix.resize( 3, 3, true);

        rMatrix( 0 , 2 ) = 0.0;
        rMatrix( 1 , 2 ) = 0.0;

        rMatrix( 2 , 0 ) = 0.0;
        rMatrix( 2 , 1 ) = 0.0;

        rMatrix( 2 , 2 ) = 1.0;

    }
    else if(rMatrix.size1() != 3 && rMatrix.size2() != 3)
    {

        KRATOS_THROW_ERROR( std::invalid_argument,"Matrix Dimensions are not correct ", "" )

    }

    return rMatrix;
}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElastic3DLaw::GetLawFeatures(Features& rFeatures)
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

bool HyperElastic3DLaw::CheckParameters(Parameters& rValues)
{
    return rValues.CheckAllParameters();
}



int HyperElastic3DLaw::Check(const Properties& rMaterialProperties,
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
