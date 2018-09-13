//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_plastic_3d_law.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************
HenckyElasticPlastic3DLaw::HenckyElasticPlastic3DLaw()
    : HyperElastic3DLaw()
{

}



HenckyElasticPlastic3DLaw::HenckyElasticPlastic3DLaw(MPMFlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : HyperElastic3DLaw( )
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyElasticPlastic3DLaw::HenckyElasticPlastic3DLaw(const HenckyElasticPlastic3DLaw&  rOther)
    : HyperElastic3DLaw(rOther)
    ,mElasticLeftCauchyGreen(rOther.mElasticLeftCauchyGreen)
    ,mpYieldCriterion(rOther.mpYieldCriterion)
    ,mpHardeningLaw(rOther.mpHardeningLaw)
{
    mpMPMFlowRule       = rOther.mpMPMFlowRule->Clone();
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyElasticPlastic3DLaw::Clone() const
{
    HenckyElasticPlastic3DLaw::Pointer p_clone(new HenckyElasticPlastic3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlastic3DLaw::~HenckyElasticPlastic3DLaw()
{
}
//***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
//************************************************************************************

bool HenckyElasticPlastic3DLaw::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool HenckyElasticPlastic3DLaw::Has( const Variable<Vector>& rThisVariable )
{
    return false;
}

bool HenckyElasticPlastic3DLaw::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}


//***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************

double& HenckyElasticPlastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if (rThisVariable==DELTA_PLASTIC_STRAIN)
    {
        const MPMFlowRule::InternalVariables& InternalVariables = mpMPMFlowRule->GetInternalVariables();
        rValue=InternalVariables.DeltaPlasticStrain;
    }
    
    if (rThisVariable==PLASTIC_STRAIN)
    {
        const MPMFlowRule::InternalVariables& InternalVariables = mpMPMFlowRule->GetInternalVariables();
        rValue=InternalVariables.EquivalentPlasticStrain;
    }

    if (rThisVariable==MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN)
    {
        const MPMFlowRule::InternalVariables& InternalVariables = mpMPMFlowRule->GetInternalVariables();
        rValue=InternalVariables.DeltaPlasticVolumetricStrain;
    }

    if (rThisVariable==MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN)
    {
        const MPMFlowRule::InternalVariables& InternalVariables = mpMPMFlowRule->GetInternalVariables();
        rValue=InternalVariables.AccumulatedPlasticVolumetricStrain;
    }

    if (rThisVariable==MP_DELTA_PLASTIC_DEVIATORIC_STRAIN)
    {
        const MPMFlowRule::InternalVariables& InternalVariables = mpMPMFlowRule->GetInternalVariables();
        rValue=InternalVariables.DeltaPlasticDeviatoricStrain;
    }

    if (rThisVariable==MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN)
    {
        const MPMFlowRule::InternalVariables& InternalVariables = mpMPMFlowRule->GetInternalVariables();
        rValue=InternalVariables.AccumulatedPlasticDeviatoricStrain;
    }

    if (rThisVariable==MIU)
    {
        rValue=mPlasticRegion;
    }
    return( rValue );
}

Vector& HenckyElasticPlastic3DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    return( rValue );
}

Matrix& HenckyElasticPlastic3DLaw::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return( rValue );
}


//***********************SET VALUE: DOUBLE - VECTOR - MATRIX**************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
        const ProcessInfo& rCurrentProcessInfo )
{

}

void HenckyElasticPlastic3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo )
{

}

void HenckyElasticPlastic3DLaw::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
        const ProcessInfo& rCurrentProcessInfo )
{

}


//************************************************************************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::InitializeMaterial(const Properties& rProps,
        const GeometryType& rGeom,
        const Vector& rShapeFunctionsValues)
{
    mDeterminantF0                = 1;
    mInverseDeformationGradientF0 = identity_matrix<double> (3);
    mElasticLeftCauchyGreen       = identity_matrix<double> (3);

    mPlasticRegion = 0;

    mpMPMFlowRule->InitializeMaterial(mpYieldCriterion, mpHardeningLaw, rProps);
}
//************************************************************************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::InitializeSolutionStep( const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry, //this is just to give the array of nodes
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    mpHardeningLaw->SetProperties(rMaterialProperties);
}

//************************************************************************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::FinalizeSolutionStep( const Properties& rMaterialProperties,
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


void HenckyElasticPlastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    // Preliminary processes
    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();

    Matrix DeformationGradientF     = rValues.GetDeformationGradientF();
    double DeterminantF             = rValues.GetDeterminantF();

    const GeometryType&  DomainGeometry    = rValues.GetElementGeometry ();
    const Vector&        ShapeFunctions    = rValues.GetShapeFunctionsValues ();

    Vector& StrainVector                   = rValues.GetStrainVector();
    Vector& StressVector                   = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix             = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    PlasticMaterialResponseVariables PlasticVariables;
    ElasticVariables.Identity = identity_matrix<double> ( 3 );

    ElasticVariables.SetElementGeometry(DomainGeometry);
    ElasticVariables.SetShapeFunctionsValues(ShapeFunctions);
        
    MPMFlowRule::RadialReturnVariables ReturnMappingVariables;
    // ReturnMappingVariables.initialize(); //it has to be called at the start
    ReturnMappingVariables.clear();

    // Initialize variables from the process information
    ReturnMappingVariables.DeltaTime = CurrentProcessInfo[DELTA_TIME];

    if(CurrentProcessInfo[IMPLEX] == 1)
        ReturnMappingVariables.Options.Set(MPMFlowRule::IMPLEX_ACTIVE,true);
    else
        ReturnMappingVariables.Options.Set(MPMFlowRule::IMPLEX_ACTIVE,false);

    //1.-Determinant of the Total Deformation Gradient -- detF
    ElasticVariables.DeterminantF = DeterminantF;

    //2.-Compute Incremental DeformationGradient (in 3D) -- F
    ElasticVariables.DeformationGradientF = DeformationGradientF;
    ElasticVariables.DeformationGradientF = Transform2DTo3D(ElasticVariables.DeformationGradientF);
    ElasticVariables.DeformationGradientF = prod(ElasticVariables.DeformationGradientF,mInverseDeformationGradientF0);
     
    //3.-Left Cauchy Green tensor b: (stored in the CauchyGreenMatrix) -- B = FF^T
    ElasticVariables.CauchyGreenMatrix = prod(mElasticLeftCauchyGreen,trans(ElasticVariables.DeformationGradientF));
    ElasticVariables.CauchyGreenMatrix = prod(ElasticVariables.DeformationGradientF,ElasticVariables.CauchyGreenMatrix);

    //4.-Compute the inverse of trial left stretch tensor V
    this->CalculateLeftStretchTensor(PlasticVariables.TrialLeftStretchTensor, ElasticVariables.CauchyGreenMatrix);
    double detV = MathUtils<double>::Det(PlasticVariables.TrialLeftStretchTensor);
    MathUtils<double>::InvertMatrix( PlasticVariables.TrialLeftStretchTensor, PlasticVariables.InverseTrialLeftStretchTensor, detV);

    //5.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
        // Almansi Strain -- E = 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix, StrainVector);
    }

    //6.-Calculate Total Kirchhoff stress
    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

        Matrix StressMatrix           = ZeroMatrix(3,3);
        Vector HenckyMainStrainVector = ZeroVector(3);

        // Hencky Strain -- E = 0.5 * ln(C)
        this->CalculateHenckyMainStrain(ElasticVariables.CauchyGreenMatrix, ReturnMappingVariables, HenckyMainStrainVector);

        ReturnMappingVariables.StrainMatrix = ZeroMatrix(3,3);
        ReturnMappingVariables.TrialIsoStressMatrix = ZeroMatrix(3,3);
        
        Matrix HenckyMainStrainMatrix = ZeroMatrix(3,3);
        for (unsigned int i = 0; i<3; ++i)
            HenckyMainStrainMatrix(i,i) = HenckyMainStrainVector[i];

        Matrix NewElasticLeftCauchyGreen = HenckyMainStrainMatrix;

        this->CalculatePrincipalStressTrial(ElasticVariables, rValues, ReturnMappingVariables, NewElasticLeftCauchyGreen, StressMatrix );

        //Attention!! 
        /*  When I call the return mapping function NewElasticLeftCauchyGreen represents the Hencky strain in matrix form.
            When the return mapping is finished NewElasticLeftCauchyGreen is the NEW elastic left cauchy green tensor.
            If and only if GetElasticLeftCachyGreen is a protected member of the flow rule that I am using.
            Otherwise a public member of the flow rule base class has to be added as in this case.*/
        mpMPMFlowRule->CalculateReturnMapping( ReturnMappingVariables, ElasticVariables.DeformationGradientF, StressMatrix, NewElasticLeftCauchyGreen);

        mPlasticRegion = 0;
        if( ReturnMappingVariables.Options.Is(MPMFlowRule::PLASTIC_REGION) )
        {
            mPlasticRegion = mpMPMFlowRule->GetPlasticRegion();
        }
        
        this->CorrectDomainPressure( StressMatrix, ElasticVariables);

        // Stress vector updated
        StressVector = MathUtils<double>::StressTensorToVector(StressMatrix, StressVector.size());
        
        // Calculate Constitutive Matrix related to Total Kirchhoff stress -- Dep
        ConstitutiveMatrix.clear();
        Matrix AuxConstitutiveMatrix = ZeroMatrix(6,6);

        double alfa = 0.0;
        this->CalculateElastoPlasticTangentMatrix( ReturnMappingVariables, ElasticVariables.CauchyGreenMatrix, alfa, AuxConstitutiveMatrix, ElasticVariables);
        ConstitutiveMatrix = this->SetConstitutiveMatrixToAppropiateDimension(ConstitutiveMatrix, AuxConstitutiveMatrix);

    }

    //7.-Update the variables at the end of iteration
    if( Options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) )
    {
         mpMPMFlowRule->UpdateInternalVariables ( ReturnMappingVariables );
        
        // Update final left cauchy green B_(n+1)
        mElasticLeftCauchyGreen = mpMPMFlowRule->GetElasticLeftCauchyGreen(ReturnMappingVariables);

        // Copying the update DeformationGradientF to elasticVariables
        ElasticVariables.DeformationGradientF = DeformationGradientF;
        ElasticVariables.DeformationGradientF = Transform2DTo3D(ElasticVariables.DeformationGradientF);
        
        // Update DeterminantF0
        mDeterminantF0 = DeterminantF;

        // Compute Inverse
        MathUtils<double>::InvertMatrix( ElasticVariables.DeformationGradientF, mInverseDeformationGradientF0, mDeterminantF0);
    }
    else{
        // Update necessary kinematics variable after return mapping
        this->CorrectKinematics( PlasticVariables, rValues, ReturnMappingVariables, DeterminantF, DeformationGradientF);
    }

}

void HenckyElasticPlastic3DLaw::CalculatePrincipalStressTrial(const MaterialResponseVariables & rElasticVariables, Parameters& rValues, 
    const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix)
{

    mpMPMFlowRule->CalculatePrincipalStressTrial(rReturnMappingVariables, rNewElasticLeftCauchyGreen, rStressMatrix);

}

//************************************************************************************
//************************************************************************************

void HenckyElasticPlastic3DLaw::CalculateLeftStretchTensor(Matrix& rLeftStretchTensor, const Matrix& rCauchyGreenMatrix)
{
    rLeftStretchTensor = identity_matrix<double> (3);

    Matrix EigenVectors  = ZeroMatrix(3,3);
    Vector EigenValues   = ZeroVector(3);
    Matrix SquaredEigenValuesMatrix = ZeroMatrix(3,3);

    double tol = 1e-9;
    int iter = 100;

    SolidMechanicsMathUtilities<double>::EigenVectors(rCauchyGreenMatrix, EigenVectors, EigenValues, tol, iter);

    // Squaring the Eigenvalue of B to obtain the Eigenvalue of V
    for (int i = 0; i < 3; i++)
        SquaredEigenValuesMatrix(i,i) = std::sqrt(EigenValues(i));

    noalias(rLeftStretchTensor) = prod(SquaredEigenValuesMatrix,EigenVectors);
    rLeftStretchTensor = prod(trans(EigenVectors),rLeftStretchTensor);

}

//************************************************************************************
//************************************************************************************

void HenckyElasticPlastic3DLaw::GetDomainPressure( double& rPressure, const MaterialResponseVariables& rElasticVariables)
{
    rPressure = 0.0;
    const GeometryType&  DomainGeometry =  rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues  =  rElasticVariables.GetShapeFunctionsValues();

    const unsigned int number_of_nodes  =  DomainGeometry.size();

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        rPressure += ShapeFunctionsValues[j] * DomainGeometry[j].FastGetSolutionStepValue(PRESSURE); 
    }
}
//************************************************************************************
//************************************************************************************

void HenckyElasticPlastic3DLaw::CorrectDomainPressure( Matrix& rStressMatrix, const MaterialResponseVariables & rElasticVariables)
{

}

//************************************************************************************
//************************************************************************************

void HenckyElasticPlastic3DLaw::CorrectKinematics(const PlasticMaterialResponseVariables& rPlasticVariables, Parameters & rValues, MPMFlowRule::RadialReturnVariables rReturnMappingVariables, double& rDeterminantF, Matrix& rDeformationGradientF )
{
    // Update Elastic Left Cuachy Green B^e_(k+1) 
    Matrix ElasticLeftCauchyGreen = mpMPMFlowRule->GetElasticLeftCauchyGreen(rReturnMappingVariables);

    // Compute left stretch tensor V^e_(k+1) 
    Matrix LeftStretchTensor;
    this->CalculateLeftStretchTensor(LeftStretchTensor, ElasticLeftCauchyGreen);

    // Update Deformation Gradient F^e_(k+1) = V^e_(k+1) R^e_(k+1)
    rDeformationGradientF = Transform2DTo3D(rDeformationGradientF);
    rDeformationGradientF = prod(rPlasticVariables.InverseTrialLeftStretchTensor, rDeformationGradientF);
    rDeformationGradientF = prod(LeftStretchTensor, rDeformationGradientF);
    rDeformationGradientF = this->SetMatrixToAppropriateDimension(rDeformationGradientF);

    // Update Determinant of Deformation Gradient F
    rDeterminantF        = MathUtils<double>::Det(rDeformationGradientF);

    // Set Deformation Gradient and Determinant
    rValues.SetDeformationGradientF(rDeformationGradientF);
    rValues.SetDeterminantF(rDeterminantF);

    KRATOS_ERROR_IF(rDeterminantF <= 0) << "HenckyElasticPlastic3DLaw::CorrectKinematics:: Updated DetF <= 0! " << rDeterminantF << std::endl;

}
//************************************************************************************
//************************************************************************************

void HenckyElasticPlastic3DLaw::CalculateElastoPlasticTangentMatrix( const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& rAlpha, Matrix& rElastoPlasticTangentMatrix, const MaterialResponseVariables& rElasticVariables )
{
    mpMPMFlowRule->ComputeElastoPlasticTangentMatrix( rReturnMappingVariables,  rNewElasticLeftCauchyGreen, rAlpha, rElastoPlasticTangentMatrix);
}

//************************************************************************************
//************************************************************************************

Matrix HenckyElasticPlastic3DLaw::SetMatrixToAppropriateDimension(Matrix& rMatrix){
    if(rMatrix.size1() == 3){
        return rMatrix;
    }
    else if (rMatrix.size1() == 2) { 
        rMatrix = Transform2DTo3D(rMatrix); 
        return rMatrix;
    }
    else KRATOS_ERROR << "Wrong size of input matrix:: Conversion is unknown!" << std::endl;
}

Vector HenckyElasticPlastic3DLaw::SetStressMatrixToAppropiateVectorDimension(Vector& rStressVector, const Matrix& rStressMatrix )
{
    rStressVector(0) = rStressMatrix(0,0);
    rStressVector(1) = rStressMatrix(1,1);
    rStressVector(2) = rStressMatrix(2,2);
    rStressVector(3) = rStressMatrix(0,1);
    rStressVector(4) = rStressMatrix(1,2);
    rStressVector(5) = rStressMatrix(0,2);
    return rStressVector;
}

Matrix HenckyElasticPlastic3DLaw::SetConstitutiveMatrixToAppropiateDimension(Matrix& rConstitutiveMatrix, const Matrix& rElastoPlasticTangentMatrix)
{
    rConstitutiveMatrix = rElastoPlasticTangentMatrix;
    return rConstitutiveMatrix;
}

/** CalculateEigenbases is a function which calculate the matrix with the eigenbases of the
* three principal directions
*/
Matrix HenckyElasticPlastic3DLaw::CalculateEigenbases(const MPMFlowRule::RadialReturnVariables& rReturnMappingVariables, Matrix& rEigenbasesMatrix)
{
    //1- compute eigendirections
    Vector N1 = ZeroVector(3);
    Vector N2 = ZeroVector(3);
    Vector N3 = ZeroVector(3);
    
    for (unsigned int i = 0; i<3; ++i)
    {
        N1(i) = rReturnMappingVariables.MainDirections(i,0);
        N2(i) = rReturnMappingVariables.MainDirections(i,1);
        N3(i) = rReturnMappingVariables.MainDirections(i,2);
    }

    //2- compute the eigenbases matrix
    Matrix M1 = ZeroMatrix(3,3);
    Matrix M2 = ZeroMatrix(3,3);
    Matrix M3 = ZeroMatrix(3,3);
    M1 = MathUtils<double>::TensorProduct3(N1, N1);
    M2 = MathUtils<double>::TensorProduct3(N2, N2);
    M3 = MathUtils<double>::TensorProduct3(N3, N3);

    // EigenbasesMatrix[3x9]
    for(unsigned int i = 0; i<3; i++)
    {
        for(unsigned int j = 0; j<3; j++)
        {
            int index_j = j+3;
            int index_j2 = j+6;

            rEigenbasesMatrix(i,j) = M1(i,j);
            rEigenbasesMatrix(i, index_j) = M2(i,j);
            rEigenbasesMatrix(i, index_j2) = M3(i,j);
        }
    }
    return rEigenbasesMatrix;
}

void HenckyElasticPlastic3DLaw::MyTensorProduct ( const Matrix& rMA, const Matrix& rMB,
        Matrix& rEigenbasesProductMatrix)
{
    rEigenbasesProductMatrix.clear();

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rEigenbasesProductMatrix( i, j ) = TensorComponent(rEigenbasesProductMatrix( i, j ), rMA, rMB,
                                               this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

double& HenckyElasticPlastic3DLaw::TensorComponent(double & rCabcd,
        const Matrix& rMA, const Matrix& rMB,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{
    rCabcd =(rMA(a,b)*rMB(c,d));

    return rCabcd;
}

void HenckyElasticPlastic3DLaw::MyTensorProduct2 ( const Matrix& rMA, const Matrix& rMB,
        Matrix& rEigenbasesProductMatrix)
{
    rEigenbasesProductMatrix.clear();

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rEigenbasesProductMatrix( i, j ) = TensorComponent2(rEigenbasesProductMatrix( i, j ), rMA, rMB,
                                               this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }
}

double& HenckyElasticPlastic3DLaw::TensorComponent2(double & rCabcd,
        const Matrix& rMA, const Matrix& rMB,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{
    rCabcd =(rMA(a,b)*rMB(c,d) + rMB(a,b)*rMA(c,d));

    return rCabcd;
}

void HenckyElasticPlastic3DLaw::MyTensorProduct3 ( const Matrix& rMA, Matrix& rEigenbasesProductMatrix)
{

    rEigenbasesProductMatrix.clear();

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rEigenbasesProductMatrix( i, j ) = TensorComponent3(rEigenbasesProductMatrix( i, j ), rMA,
                                               this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

double& HenckyElasticPlastic3DLaw::TensorComponent3(double & rCabcd,
        const Matrix& rMA,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{
    rCabcd =rMA(a,b)*rMA(c,d) - 0.5*(( rMA(a,c)*rMA(b,d))+( rMA(a,d)*rMA(b,c)));

    return rCabcd;
}

void HenckyElasticPlastic3DLaw::MyTensorProduct4 ( const Matrix& rMA, Matrix& rEigenbasesProductMatrix)
{
    rEigenbasesProductMatrix.clear();

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rEigenbasesProductMatrix( i, j ) = TensorComponent4(rEigenbasesProductMatrix( i, j ), rMA,
                                               this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }
    }
}

double& HenckyElasticPlastic3DLaw::TensorComponent4(double & rCabcd,
        const Matrix& rMA,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{
    rCabcd = 0.5 *(rMA(a,c)*rMA(b,d) + rMA(a,d)*rMA(b,c)) - rMA(a,b)*rMA(c,d) ;

    return rCabcd;
}

//************************************************************************************
//************************************************************************************


Vector& HenckyElasticPlastic3DLaw::GetStressVectorFromMatrix(const Matrix& rStressMatrix,
        Vector& rMainStress,
        const Matrix& rEigenVectors)
{
    Matrix auxMatrix = ZeroMatrix(3,3);
    auxMatrix = prod( rStressMatrix, trans(rEigenVectors));
    auxMatrix = prod( (rEigenVectors), auxMatrix);

    rMainStress = ZeroVector(3);
    for (unsigned int i = 0; i<3; i++)
        rMainStress[i] = auxMatrix(i,i);

    return rMainStress;
}



//************************************************************************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::CalculateHenckyMainStrain(const Matrix& rCauchyGreenMatrix,
        MPMFlowRule::RadialReturnVariables& rReturnMappingVariables,
        Vector& rMainStrain)
{
    Matrix EigenVectors  = ZeroMatrix(3,3);
    Vector EigenValues   = ZeroVector(3);

    double tol = 1e-9;
    int iter = 100;

    SolidMechanicsMathUtilities<double>::EigenVectors(rCauchyGreenMatrix, EigenVectors, EigenValues, tol, iter);
    rReturnMappingVariables.MainDirections     = EigenVectors;

    for (unsigned int i = 0; i<3; ++i)
        rMainStrain[i] = 0.50 * std::log(EigenValues[i]);
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HenckyElasticPlastic3DLaw::GetLawFeatures(Features& rFeatures)
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

bool HenckyElasticPlastic3DLaw::CheckParameters(Parameters& rValues)
{
    return rValues.CheckAllParameters();
}

int HenckyElasticPlastic3DLaw::Check(const Properties& rMaterialProperties,
                                     const GeometryType& rElementGeometry,
                                     const ProcessInfo& rCurrentProcessInfo)
{

    // Verify Positive Density
    KRATOS_ERROR_IF(DENSITY.Key() == 0 || rMaterialProperties[DENSITY]<0.00) << "DENSITY has Key zero or invalid value " << std::endl;

    return 0;
}

} // namespace Kratos
