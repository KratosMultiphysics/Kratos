//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//
//  References:      This class is adapted from applications/SolidMechanicsApplication/custom_constitutive/hyperelastic_3D_law.cpp


// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_3D_law.hpp"

#include "particle_mechanics_application_variables.h"

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
        mDeterminantF0 = rValue;
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
  mInverseDeformationGradientF0 = IdentityMatrix(3);
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

//*****************************MATERIAL RESPONSES*************************************
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
    ElasticVariables.Identity = IdentityMatrix( 3 );

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

}


//************************************************************************************
//************************************************************************************

void HyperElastic3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    this->CalculateMaterialResponseKirchhoff (rValues);

    const double& DeterminantF  = rValues.GetDeterminantF();
    Vector& StressVector        = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix  = rValues.GetConstitutiveMatrix();

    //Set to cauchy Stress:
    StressVector       /= DeterminantF;
    ConstitutiveMatrix /= DeterminantF;

}


//***********************************UPDATE*******************************************
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

}

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

}


//******************************* COMPUTE DOMAIN TEMPERATURE  ************************
//************************************************************************************


double &  HyperElastic3DLaw::CalculateDomainTemperature (const MaterialResponseVariables & rElasticVariables,
								double & rTemperature)
{

    // Temperature from nodes
    const GeometryType& DomainGeometry = rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues = rElasticVariables.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();

    // Default Value
    rTemperature=0;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        if( DomainGeometry[j].SolutionStepsDataHas(TEMPERATURE) )
            rTemperature += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(TEMPERATURE);
    }

    return rTemperature;
}


//***************************** COMPUTE VOLUMETRIC FACTOR ****************************
//************************************************************************************


double & HyperElastic3DLaw::CalculateVolumetricFactor (const MaterialResponseVariables & rElasticVariables,
						       double & rFactor)
{
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


//******************************* COMPUTE VOLUMETRIC STRESS  *************************
//************************************************************************************

void HyperElastic3DLaw::CalculateVolumetricStress(const MaterialResponseVariables & rElasticVariables,
						  Vector& rVolStressVector )
{

    Matrix VolStressMatrix ( 3 , 3 );

    double Pressure = 0;
    Pressure = this->CalculateVolumetricPressure (rElasticVariables, Pressure);

    VolStressMatrix = rElasticVariables.DeterminantF * Pressure * rElasticVariables.CauchyGreenMatrix;
    rVolStressVector = MathUtils<double>::StressTensorToVector(VolStressMatrix,rVolStressVector.size());

}


//******************************* COMPUTE ISOCHORIC STRESS  **************************
//************************************************************************************
void HyperElastic3DLaw::CalculateIsochoricStress( const MaterialResponseVariables & rElasticVariables,
						  StressMeasure rStressMeasure,
						  Vector& rIsoStressVector )
{

    Matrix IsoStressMatrix ( 3, 3 );

    //note.- rElasticVariables.traceCG is "traceCG"
    if(rStressMeasure == StressMeasure_PK2)
    {
        //rElasticVariables.CauchyGreenMatrix is InverseRightCauchyGreen
        //2.-Incompressible part of the 2nd Piola Kirchhoff Stress Matrix
        IsoStressMatrix  = (rElasticVariables.Identity - (rElasticVariables.traceCG/3.0)*rElasticVariables.CauchyGreenMatrix );
        IsoStressMatrix *= rElasticVariables.LameMu*pow(rElasticVariables.DeterminantF,(-2.0/3.0));
    }

    if(rStressMeasure == StressMeasure_Kirchhoff)
    {
        //rElasticVariables.CauchyGreenMatrix is LeftCauchyGreen
        //2.-Incompressible part of the Kirchhoff Stress Matrix
        IsoStressMatrix  = (rElasticVariables.CauchyGreenMatrix - (rElasticVariables.traceCG/3.0)*rElasticVariables.Identity );
        IsoStressMatrix *= rElasticVariables.LameMu*pow(rElasticVariables.DeterminantF,(-2.0/3.0));
    }

    rIsoStressVector = MathUtils<double>::StressTensorToVector(IsoStressMatrix,rIsoStressVector.size());
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

    //1.Elastic constitutive tensor component:
    rCabcd =(rElasticVariables.LameLambda*auxiliar1*rElasticVariables.CauchyGreenMatrix(a,b)*rElasticVariables.CauchyGreenMatrix(c,d));
    rCabcd+=((2*rElasticVariables.LameMu-rElasticVariables.LameLambda*auxiliar2)*0.5*(rElasticVariables.CauchyGreenMatrix(a,c)*rElasticVariables.CauchyGreenMatrix(b,d)+rElasticVariables.CauchyGreenMatrix(a,d)*rElasticVariables.CauchyGreenMatrix(b,c)));

    return rCabcd;
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
        KRATOS_ERROR << "Matrix Dimensions are not correct !" << std::endl;
    }

    return rMatrix;
}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElastic3DLaw::GetLawFeatures(Features& rFeatures)
{
    // Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( FINITE_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	// Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	// Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	// Set the spacedimension
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

    KRATOS_ERROR_IF (YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS]<= 0.00) << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const double tolerance = 10.e-7;
    const bool check = bool( (nu > 0.5-tolerance ) || (nu < (-1.0 + tolerance)) );

    KRATOS_ERROR_IF (POISSON_RATIO.Key() == 0 || check==true) << "POISSON_RATIO has Key zero invalid value " << std::endl;
    KRATOS_ERROR_IF (DENSITY.Key() == 0 || rMaterialProperties[DENSITY]<0.00) << "DENSITY has Key zero or invalid value " << std::endl;

    return 0;

}

} // Namespace Kratos
