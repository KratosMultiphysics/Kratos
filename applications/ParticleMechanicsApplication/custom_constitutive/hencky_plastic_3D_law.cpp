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
#include <cmath>

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_plastic_3D_law.hpp"
#include "custom_utilities/particle_mechanics_math_utilities.h"
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
    if (rThisVariable==MP_DELTA_PLASTIC_STRAIN)
    {
        const MPMFlowRule::InternalVariables& InternalVariables = mpMPMFlowRule->GetInternalVariables();
        rValue=InternalVariables.DeltaPlasticStrain;
    }

    if (rThisVariable==MP_EQUIVALENT_PLASTIC_STRAIN)
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
    mInverseDeformationGradientF0 = IdentityMatrix(3);
    mElasticLeftCauchyGreen       = IdentityMatrix(3);

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
    Flags &options=rValues.GetOptions();

    const ProcessInfo& current_process_info = rValues.GetProcessInfo();

    const Matrix deformation_gradient_F     = rValues.GetDeformationGradientF();
    const double determinant_F              = rValues.GetDeterminantF();

    const GeometryType&  domain_geometry    = rValues.GetElementGeometry ();
    const Vector&        shape_functions    = rValues.GetShapeFunctionsValues ();

    Vector& strain_vector                   = rValues.GetStrainVector();
    Vector& stress_vector                   = rValues.GetStressVector();
    Matrix& constitutive_matrix             = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    PlasticMaterialResponseVariables PlasticVariables;
    ElasticVariables.Identity = IdentityMatrix(3);

    ElasticVariables.SetElementGeometry(domain_geometry);
    ElasticVariables.SetShapeFunctionsValues(shape_functions);

    MPMFlowRule::RadialReturnVariables ReturnMappingVariables;
    // ReturnMappingVariables.initialize(); //it has to be called at the start
    ReturnMappingVariables.clear();

    // Initialize variables from the process information
    ReturnMappingVariables.DeltaTime = current_process_info[DELTA_TIME];

    if(current_process_info[IMPLEX] == 1)
        ReturnMappingVariables.Options.Set(MPMFlowRule::IMPLEX_ACTIVE,true);
    else
        ReturnMappingVariables.Options.Set(MPMFlowRule::IMPLEX_ACTIVE,false);

    //1.-Determinant of the Total Deformation Gradient -- detF
    ElasticVariables.DeterminantF = determinant_F;

    //2.-Compute Incremental DeformationGradient (in 3D) -- F
    ElasticVariables.DeformationGradientF = deformation_gradient_F;
    ElasticVariables.DeformationGradientF = Transform2DTo3D(ElasticVariables.DeformationGradientF);
    ElasticVariables.DeformationGradientF = prod(ElasticVariables.DeformationGradientF,mInverseDeformationGradientF0);

    //3.-Left Cauchy Green tensor b: (stored in the CauchyGreenMatrix) -- B = FF^T
    ElasticVariables.CauchyGreenMatrix = prod(mElasticLeftCauchyGreen,trans(ElasticVariables.DeformationGradientF));
    ElasticVariables.CauchyGreenMatrix = prod(ElasticVariables.DeformationGradientF,ElasticVariables.CauchyGreenMatrix);

    //4.-Almansi Strain:
    if(options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
        // Almansi Strain -- E = 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix, strain_vector);
    }

    //5.-Calculate Total Kirchhoff stress
    if( options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

        Matrix stress_matrix             = ZeroMatrix(3);
        Vector hencky_main_strain_vector = ZeroVector(3);

        // Hencky Strain -- E = 0.5 * ln(C)
        this->CalculateHenckyMainStrain(ElasticVariables.CauchyGreenMatrix, ReturnMappingVariables, hencky_main_strain_vector);

        ReturnMappingVariables.StrainMatrix = ZeroMatrix(3);
        ReturnMappingVariables.TrialIsoStressMatrix = ZeroMatrix(3);

        Matrix hencky_main_strain_matrix = ZeroMatrix(3);
        for (unsigned int i = 0; i<3; ++i)
            hencky_main_strain_matrix(i,i) = hencky_main_strain_vector[i];

        this->CalculatePrincipalStressTrial(ElasticVariables, rValues, ReturnMappingVariables, hencky_main_strain_matrix, stress_matrix );

        //Attention!!
        /*  When I call the return mapping function NewElasticLeftCauchyGreen represents the Hencky strain in matrix form.
            When the return mapping is finished NewElasticLeftCauchyGreen is the NEW elastic left cauchy green tensor.
            If and only if GetElasticLeftCachyGreen is a protected member of the flow rule that I am using.
            Otherwise a public member of the flow rule base class has to be added as in this case.*/
        mpMPMFlowRule->CalculateReturnMapping( ReturnMappingVariables, ElasticVariables.DeformationGradientF, stress_matrix, hencky_main_strain_matrix);

        mPlasticRegion = 0;
        if( ReturnMappingVariables.Options.Is(MPMFlowRule::PLASTIC_REGION) )
        {
            mPlasticRegion = mpMPMFlowRule->GetPlasticRegion();
        }

        this->CorrectDomainPressure( stress_matrix, ElasticVariables);

        // Stress vector updated
        stress_vector = MathUtils<double>::StressTensorToVector(stress_matrix, stress_vector.size());

        // Calculate Constitutive Matrix related to Total Kirchhoff stress -- Dep
        constitutive_matrix.clear();
        Matrix aux_constitutive_matrix = ZeroMatrix(6);

        double alfa = 0.0;
        this->CalculateElastoPlasticTangentMatrix( ReturnMappingVariables, ElasticVariables.CauchyGreenMatrix, alfa, aux_constitutive_matrix, ElasticVariables);
        constitutive_matrix = this->SetConstitutiveMatrixToAppropiateDimension(constitutive_matrix, aux_constitutive_matrix);

    }

    //6.-Update the variables at the end of iteration
    if( options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) )
    {
         mpMPMFlowRule->UpdateInternalVariables ( ReturnMappingVariables );

        // Update final left cauchy green B_(n+1)
        mElasticLeftCauchyGreen = mpMPMFlowRule->GetElasticLeftCauchyGreen(ReturnMappingVariables);

        // Copying the update DeformationGradientF to elasticVariables
        ElasticVariables.DeformationGradientF = deformation_gradient_F;
        ElasticVariables.DeformationGradientF = Transform2DTo3D(ElasticVariables.DeformationGradientF);

        // Update DeterminantF0
        mDeterminantF0 = determinant_F;

        // Compute Inverse
        MathUtils<double>::InvertMatrix( ElasticVariables.DeformationGradientF, mInverseDeformationGradientF0, mDeterminantF0);
    }

}

void HenckyElasticPlastic3DLaw::CalculatePrincipalStressTrial(const MaterialResponseVariables & rElasticVariables, Parameters& rValues,
    const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix)
{

    mpMPMFlowRule->CalculatePrincipalStressTrial(rReturnMappingVariables, rNewElasticLeftCauchyGreen, rStressMatrix);

}


//************************************************************************************
//************************************************************************************

void HenckyElasticPlastic3DLaw::GetDomainPressure( double& rPressure, const MaterialResponseVariables& rElasticVariables)
{
    rPressure = 0.0;
    const GeometryType&  domain_geometry =  rElasticVariables.GetElementGeometry();
    const Vector& shape_functions        =  rElasticVariables.GetShapeFunctionsValues();

    const unsigned int number_of_nodes  =  domain_geometry.size();

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        rPressure += shape_functions[j] * domain_geometry[j].FastGetSolutionStepValue(PRESSURE);
    }
}
//************************************************************************************
//************************************************************************************

void HenckyElasticPlastic3DLaw::CorrectDomainPressure( Matrix& rStressMatrix, const MaterialResponseVariables & rElasticVariables)
{

}

//************************************************************************************
//************************************************************************************

void HenckyElasticPlastic3DLaw::CalculateElastoPlasticTangentMatrix( const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& rAlpha, Matrix& rElastoPlasticTangentMatrix, const MaterialResponseVariables& rElasticVariables )
{
    mpMPMFlowRule->ComputeElastoPlasticTangentMatrix( rReturnMappingVariables,  rNewElasticLeftCauchyGreen, rAlpha, rElastoPlasticTangentMatrix);
}

//************************************************************************************
//************************************************************************************

Vector HenckyElasticPlastic3DLaw::SetStressMatrixToAppropiateVectorDimension(Vector& rStressVector, const Matrix& rStressMatrix )
{
    rStressVector[0] = rStressMatrix(0,0);
    rStressVector[1] = rStressMatrix(1,1);
    rStressVector[2] = rStressMatrix(2,2);
    rStressVector[3] = rStressMatrix(0,1);
    rStressVector[4] = rStressMatrix(1,2);
    rStressVector[5] = rStressMatrix(0,2);
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
        N1[i] = rReturnMappingVariables.MainDirections(i,0);
        N2[i] = rReturnMappingVariables.MainDirections(i,1);
        N3[i] = rReturnMappingVariables.MainDirections(i,2);
    }

    //2- compute the eigenbases matrix
    Matrix M1 = ZeroMatrix(3);
    Matrix M2 = ZeroMatrix(3);
    Matrix M3 = ZeroMatrix(3);
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
    Matrix auxMatrix = ZeroMatrix(3);
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
    Matrix eigen_vectors  = ZeroMatrix(3);
    Vector eigen_values   = ZeroVector(3);

    double tol = 1e-9;
    int iter = 100;

    ParticleMechanicsMathUtilities<double>::EigenVectors(rCauchyGreenMatrix, eigen_vectors, eigen_values, tol, iter);
    rReturnMappingVariables.MainDirections     = eigen_vectors;

    for (unsigned int i = 0; i<3; ++i)
        rMainStrain[i] = 0.50 * std::log(eigen_values[i]);
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


//************************************************************************************
//************************************************************************************

int HenckyElasticPlastic3DLaw::Check(const Properties& rMaterialProperties,
                                     const GeometryType& rElementGeometry,
                                     const ProcessInfo& rCurrentProcessInfo)
{

    // Verify Positive Density
    KRATOS_ERROR_IF(DENSITY.Key() == 0 || rMaterialProperties[DENSITY]<0.00) << "DENSITY has Key zero or invalid value " << std::endl;

    return 0;
}

} // namespace Kratos
