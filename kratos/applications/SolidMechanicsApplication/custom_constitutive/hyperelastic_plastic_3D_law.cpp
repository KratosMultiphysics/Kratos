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

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticPlastic3DLaw::HyperElasticPlastic3DLaw(const HyperElasticPlastic3DLaw& rOther)
    : ConstitutiveLaw()
	,mElasticLeftCauchyGreen(rOther.mElasticLeftCauchyGreen)
	,mpFlowRule(rOther.mpFlowRule)
	,mpHardeningLaw(rOther.mpHardeningLaw)
{

	mpYieldCriterion =  rOther.mpYieldCriterion.clone();
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticPlastic3DLaw::Clone() const
{
    HyperElasticPlastic3DLaw::Pointer p_clone(new HyperElasticPlastic3DLaw(*this));
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


void HyperElasticPlastic3DLaw::InitializeMaterial( const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues )
{
  mpFlowRule.InitializeMaterial( mpYieldCriterion, mpHardeningLaw, props );
  mElasticLeftCauchyGreen = identity_matrix<double> (3);
}

//************************************************************************************
//************************************************************************************


void HyperElasticPlastic3DLaw::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues,
        const ProcessInfo& CurrentProcessInfo)
{

}

//************************************************************************************
//************************************************************************************


void HyperElasticPlastic3DLaw::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues,
        const ProcessInfo& CurrentProcessInfo)
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

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const Matrix& DeformationGradientF    = rValues.GetDeformationGradientF();
    const double& DeterminantF            = rValues.GetDeterminantF();

    Vector& StrainVector                  = rValues.GetStrainVector();
    Matrix& DeformationGradientF0         = rValues.GetDeformationGradientF0();
    double& DeterminantF0                 = rValues.GetDeterminantF0();

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.IdentityMatrix = identity_matrix<double> ( 3 );

    //1.- Lame constants
    const double& YoungModulus        = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient  = MaterialProperties[POISSON_RATIO];

    ElasticVariables.LameLambda      = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    ElasticVariables.LameMu          =  YoungModulus/(2*(1+PoissonCoefficient));

    //-----------------------------//
    //OPTION 1: ( initial configuration )
    if( Options.Is( ConstitutiveLaw::INITIAL_CONFIGURATION ) )
    {

        //2.-Total Deformation Gradient
        Matrix TotalDeformationGradientF0  = DeformationGradientF;
        TotalDeformationGradientF0         = DeformationGradient3D( TotalDeformationGradientF0 );

        //3.-Determinant of the Total Deformation Gradient
        ElasticVariables.DeterminantF0 = DeterminantF0 * DeterminantF;

        //4.-Right Cauchy Green
        Matrix RightCauchyGreen = prod(trans(TotalDeformationGradientF0),TotalDeformationGradientF0);

        //5.-Inverse of the Right Cauchy-Green tensor C: (stored in the CauchyGreenMatrix)
        double trace_C = 0;
        ElasticVariables.CauchyGreenMatrix( 3, 3 );
        MathUtils<double>::InvertMatrix( RightCauchyGreen, ElasticVariables.CauchyGreenMatrix, trace_C);

        //6.-Green-Lagrange Strain:
        if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
        {
            this->CalculateGreenLagrangeStrain(RightCauchyGreen, StrainVector);
        }

        //7.-Calculate Total PK2 stress
        if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
        {

            this->CalculateStress( ElasticVariables, StressMeasure_PK2, StressVector );
        }

        if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
        {

            this->CalculateConstitutiveMatrix ( ElasticVariables, ConstitutiveMatrix );
        }

    }

    //-----------------------------//
    //OPTION 2: ( last known configuration : updated lagrangian approach only )
    if( Options.Is( ConstitutiveLaw::LAST_KNOWN_CONFIGURATION )  || Options.Is( ConstitutiveLaw::FINAL_CONFIGURATION ))
    {

        //Determinant of the Total Deformation Gradient
        ElasticVariables.DeterminantF0 = DeterminantF0 * DeterminantF;

        if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
        {

            //Left Cauchy-Green tensor b
            Matrix TotalDeformationGradientF0  = prod(DeformationGradientF, DeformationGradientF0);
            TotalDeformationGradientF0         = DeformationGradient3D( TotalDeformationGradientF0 );
            ElasticVariables.CauchyGreenMatrix = prod(TotalDeformationGradientF0,trans(TotalDeformationGradientF0));

            //Almansi Strain:
            if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
            {
                // e= 0.5*(1-invbT*invb)
                this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix,StrainVector);
            }


            this->CalculateStress( ElasticVariables, StressMeasure_Kirchhoff, StressVector );

            TransformStresses(StressVector, DeformationGradientF, DeterminantF, StressMeasure_Kirchhoff, StressMeasure_PK2); //2nd PK Stress in the last known configuration
        }

        if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
        {

            Matrix InverseDeformationGradientF ( 3, 3 );
            double DetInvF=0;
            MathUtils<double>::InvertMatrix( DeformationGradientF, InverseDeformationGradientF, DetInvF);

            ElasticVariables.CauchyGreenMatrix = ElasticVariables.IdentityMatrix;

            this->CalculateConstitutiveMatrix ( ElasticVariables, InverseDeformationGradientF, ConstitutiveMatrix );

            ConstitutiveMatrix *= DeterminantF;
        }

    }


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

    Vector& StrainVector                  = rValues.GetStrainVector();
    Matrix& DeformationGradientF0         = rValues.GetDeformationGradientF0();
    double& DeterminantF0                 = rValues.GetDeterminantF0();

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.IdentityMatrix = identity_matrix<double> ( 3 );

    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.TimeStep = CurProcessInfo[DELTA_TIME];

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

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
      this->CalculateIsochoricStress( ElasticVariables, ReturnMappingVariables, StressMeasure_Kirchhoff, SplitStressVector.Isochoric );


    //OPTION 1:
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {

        SplitStressVector.Volumetric = ZeroVector(voigtsize);

        ElasticVariables.CauchyGreenMatrix = ElasticVariables.IdentityMatrix;

        this->CalculateVolumetricStress ( ElasticVariables, DomainGeometry, ShapeFunctions, SplitStressVector.Volumetric );

        //Kirchhoff Stress:
        StressVector = SplitStressVector.Isochoric + SplitStressVector.Volumetric;

        //std::cout<<" StressVector.Isochoric "<<SplitStressVector.Isochoric<<std::endl;
        //std::cout<<" StressVector.Volumetric "<<SplitStressVector.Volumetric<<std::endl;

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

        //initialize constitutive tensors
        ConstitutiveMatrix.clear();
        SplitConstitutiveMatrix.Isochoric  = ConstitutiveMatrix;
        SplitConstitutiveMatrix.Volumetric = ConstitutiveMatrix;
	SplitConstitutiveMatrix.Plastic    = ConstitutiveMatrix;

        ElasticVariables.CauchyGreenMatrix = ElasticVariables.IdentityMatrix;

        this->CalculateIsochoricConstitutiveMatrix ( ElasticVariables, ReturnMappingVariables.TrialIsoStressVector, SplitConstitutiveMatrix.Isochoric );

        this->CalculateVolumetricConstitutiveMatrix ( ElasticVariables, DomainGeometry, ShapeFunctions, SplitConstitutiveMatrix.Volumetric );

	this->CalculatePlasticConstitutiveMatrix ( ElasticVariables, ReturnMappingVariables, SplitConstitutiveMatrix.Plastic );

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

      mElasticLeftCauchyGreen  = ( SplitStressVector.Isochoric * ( 1.0 / ElasticVariables.LameMu ) );
      mElasticLeftCauchyGreen += ( ElasticVariables.traceCG/3.0) * ElasticVariables.IdentityMatrix );
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
    Flags &Options=rValues.GetOptions();

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


//******************************* COMPUTE ISOCHORIC STRESS  **************************
//************************************************************************************



void HyperElasticPlastic3DLaw::CalculateIsochoricStress( const MaterialResponseVariables & rElasticVariables,
							 FlowRule::RadialReturnVariables & rReturnMappingVariables,
							 StressMeasure rStressMeasure,
							 Vector& rIsoStressVector )
{

    //1.-Identity build
    Matrix IsoStressMatrix ( 3, 3 );

    //note.- rElasticVariables.traceCG is "traceCG_bar"

    if(rStressMeasure == StressMeasure_PK2)
    {

        //rElasticVariables.CauchyGreenMatrix is InverseRightCauchyGreen

        //2.-Isochoric part of the 2nd Piola Kirchhoff Stress Matrix
        IsoStressMatrix  = (rElasticVariables.IdentityMatrix - (rElasticVariables.traceCG/3.0)*rElasticVariables.CauchyGreenMatrix );
        IsoStressMatrix *= rElasticVariables.LameMu;

        //std::cout<<" PK2 "<<std::endl;
    }

    if(rStressMeasure == StressMeasure_Kirchhoff)
    {

        //rElasticVariables.CauchyGreenMatrix is LeftCauchyGreen

        //2.-Isochoric part of the Kirchhoff Stress Matrix
        IsoStressMatrix  = (rElasticVariables.CauchyGreenMatrix - (rElasticVariables.traceCG/3.0)*rElasticVariables.IdentityMatrix );
        IsoStressMatrix *= rElasticVariables.LameMu;

        //std::cout<<" Kirchooff "<<std::endl;

    }

    rReturnMappingVariables.TrialIsoStressVector = MathUtils<double>::StressTensorToVector(IsoStressMatrix,rIsoStressVector.size());


    rElasticVariables.Plasticity = mpFlowRule.CalculateReturnMapping( rReturnMappingVariables, IsoStressMatrix );

    
    rIsoStressVector = MathUtils<double>::StressTensorToVector(IsoStressMatrix,rIsoStressVector.size());

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
    // rFactor[0] =  rElasticVariables.DeterminantF0*rElasticVariables.DeterminantF0;
    // rFactor[1] =  (rElasticVariables.DeterminantF0*rElasticVariables.DeterminantF0-1);
    // rFactor[2] =  BulkModulus;

    //(ln(J))
    rFactor[0] =  1.0;
    rFactor[1] =  (2.0*std::log(rElasticVariables.DeterminantF0));
    rFactor[2] =  BulkModulus;

    return rFactors;
}



//******************************* COMPUTE VOLUMETRIC STRESS  *************************
//************************************************************************************

void HyperElasticPlastic3DLaw::CalculateVolumetricStress(const MaterialResponseVariables & rElasticVariables,
        const GeometryType& rDomainGeometry,
        const Vector & rShapeFunctions,
        Vector& rVolStressVector )
{

    //1.- Declaration
    Matrix VolStressMatrix ( 3 , 3 );

    double Pressure = 0;

    Pressure = CalculateDomainPressure (rElasticVariables, Pressure);

    //2.- Volumetric part of the Kirchhoff StressMatrix from nodal pressures
    VolStressMatrix = rElasticVariables.DeterminantF0 * Pressure * rElasticVariables.CauchyGreenMatrix;


    rVolStressVector = MathUtils<double>::StressTensorToVector(VolStressMatrix,rVolStressVector.size());

}

//***********************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX************************
//************************************************************************************

void HyperElasticPlastic3DLaw::CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
        const Vector & rIsoStressVector,
        Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    Matrix IsoStressMatrix = MathUtils<double>::StressVectorToTensor( rIsoStressVector );

    static const unsigned int msIndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };


    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = IsochoricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, IsoStressMatrix,
                                          msIndexVoigt3D[i][0], msIndexVoigt3D[i][1], msIndexVoigt3D[j][0], msIndexVoigt3D[j][1]);
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
    Factors = CalculateDomainPressureFactors( rElasticVariables, Factors );

    static const unsigned int msIndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, Factors,
                                          msIndexVoigt3D[i][0], msIndexVoigt3D[i][1], msIndexVoigt3D[j][0], msIndexVoigt3D[j][1]);
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

    Matrix IsoStressMatrix = MathUtils<double>::StressVectorToTensor( rIsoStressVector );

    FlowRule::PlasticFactors ScalingFactors;
    mpFlowRule.CalculateScalingFactors( rReturnMappingVariables, ScalingFactors );

    static const unsigned int msIndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
		rConstitutiveMatrix( i, j ) = PlasticConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, IsoStressMatrix, ScalingFactors,
                                          msIndexVoigt3D[i][0], msIndexVoigt3D[i][1], msIndexVoigt3D[j][0], msIndexVoigt3D[j][1]);
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

    //Isochoric part of the hyperelastic constitutive tensor component: (J²-1)/2  -  (ln(J)/J)

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


//********************CONSTITUTIVE MATRIX ISOCHORIC COMPONENT*************************
//************************************************************************************


double& HyperElasticPlastic3DLaw::PlasticConstitutiveComponent(double & rCabcd,
        const MaterialResponseVariables & rElasticVariables,
        const Matrix & rIsoStressMatrix,
        const FlowRule::PlasticFactors & rScalingFactors,					
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{

    //Isochoric part of the hyperelastic constitutive tensor component: (J²-1)/2  -  (ln(J)/J)

    rCabcd  = (1.0/3.0)*(rElasticVariables.CauchyGreenMatrix(a,b)*rElasticVariables.CauchyGreenMatrix(c,d));
    rCabcd -= (0.5*(rElasticVariables.CauchyGreenMatrix(a,c)*rElasticVariables.CauchyGreenMatrix(b,d)+rElasticVariables.CauchyGreenMatrix(a,d)*rElasticVariables.CauchyGreenMatrix(b,c)));
    
    rCabcd *= rElasticVariables.traceCG * rElasticVariables.LameMu;

    rCabcd += (rElasticVariables.CauchyGreenMatrix(c,d)*rIsoStressMatrix(a,b) + rIsoStressMatrix(c,d)*rElasticVariables.CauchyGreenMatrix(a,b));

    rCabcd *= (-2.0/3.0) * ( (-1) * rScalingFactors.Beta1 );
    
    Cabcd  -= rScalingFactors.Beta3 * 2.0 * (rElasticVariables.LameMu * (rElasticVariables.traceCG/3.0)) * ( rScalingFactors.Normal(a,b) * rScalingFactors.Normal(c,d) );

    Cabcd  -= rScalingFactors.Beta4 * 2.0 * (rElasticVariables.LameMu * (rElasticVariables.traceCG/3.0)) * ( rScalingFactors.Normal(a,b) * rScalingFactors.Dev_Normal(c,d) );

    return rCabcd;
}

//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
//************************************************************************************

bool HyperElasticPlastic3DLaw::CheckParameters(Parameters& rValues)
{
    return rValues.CheckAllParameters();
}



int HyperElasticPlastic3DLaw::Check(const Properties& rProperties,
                             const GeometryType& rGeometry,
                             const ProcessInfo& rCurrentProcessInfo)
{

    if(YOUNG_MODULUS.Key() == 0 || rProperties[YOUNG_MODULUS]<= 0.00)
        KRATOS_ERROR(std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ","");

    const double& nu = rProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    if(POISSON_RATIO.Key() == 0 || check==true)
        KRATOS_ERROR(std::invalid_argument,"POISSON_RATIO has Key zero invalid value ","");


    if(DENSITY.Key() == 0 || rProperties[DENSITY]<0.00)
        KRATOS_ERROR(std::invalid_argument,"DENSITY has Key zero or invalid value ","");


    return 0;

}

} // Namespace Kratos
