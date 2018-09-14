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


// System includes
#include <iostream>
#include <cmath>

// External includes
#include "includes/ublas_interface.h"
#include "includes/mat_variables.h"

// Project includes
#include "custom_constitutive/flow_rules/borja_cam_clay_plastic_flow_rule.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

#include "particle_mechanics_application.h"
#include "custom_utilities/mpm_stress_principal_invariants_utility.h"

namespace Kratos
{



//************ CONSTRUCTOR ***********
BorjaCamClayPlasticFlowRule::BorjaCamClayPlasticFlowRule()
    :MPMFlowRule()
{
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

BorjaCamClayPlasticFlowRule::BorjaCamClayPlasticFlowRule(YieldCriterionPointer pYieldCriterion)
    :MPMFlowRule(pYieldCriterion)
{

}

//********* ASSIGMENT OPERATOR
BorjaCamClayPlasticFlowRule& BorjaCamClayPlasticFlowRule::operator=(BorjaCamClayPlasticFlowRule const& rOther)
{
    MPMFlowRule::operator=(rOther);
    return *this;

}



//********** COPY CONSTRUCTOR *********
BorjaCamClayPlasticFlowRule::BorjaCamClayPlasticFlowRule(BorjaCamClayPlasticFlowRule const& rOther)
    :MPMFlowRule(rOther)
{
}

//*******   CLONE ********
MPMFlowRule::Pointer BorjaCamClayPlasticFlowRule::Clone() const
{
    MPMFlowRule::Pointer p_clone(new BorjaCamClayPlasticFlowRule(*this));
    return p_clone;
}



// ********** DESTRUCTOR **************
BorjaCamClayPlasticFlowRule::~BorjaCamClayPlasticFlowRule()
{
}

void BorjaCamClayPlasticFlowRule::InitializeMaterial(YieldCriterionPointer& pYieldCriterionPointer, HardeningLawPointer& pHardeningPointer, const Properties& rProp)
{
    MPMFlowRule::InitializeMaterial(pYieldCriterionPointer, pHardeningPointer, rProp);

    mElasticPrincipalStrain = ZeroVector(3);
    mPlasticPrincipalStrain = ZeroVector(3);

    mPrincipalStressUpdated = ZeroVector(3);
    mLargeStrainBool = true;
    mRegion = 0;

    // Used to calculate Omega
    mInitialVolumetricStrain = 0.0;

    // Saved state function and its derivative
    mStateFunction = 0.0;
    mStateFunctionFirstDerivative  = ZeroVector(3);
    mStateFunctionSecondDerivative = ZeroVector(6);

    this->InitializeMaterialParameters();
}

// Initiate Material Parameters which are allowed to change
void BorjaCamClayPlasticFlowRule::InitializeMaterialParameters(){
    const double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    const double OtherSlope    = mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];

    mMaterialParameters.PreconsolidationPressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    mMaterialParameters.PlasticHardeningModulus  = mMaterialParameters.PreconsolidationPressure/ (OtherSlope-SwellingSlope);
    mMaterialParameters.ConsistencyParameter     = 0.0;
}


bool BorjaCamClayPlasticFlowRule::CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, 
    Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{
    bool PlasticityActive = false;
    rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);
    
    Vector PrincipalStress = ZeroVector(3);
    Vector MainStrain      = ZeroVector(3);
    
    for (unsigned int i = 0; i<3; ++i)
        MainStrain[i] = rNewElasticLeftCauchyGreen(i,i);

    for(unsigned int i=0; i<3; i++)
    {
        // the rStressMatrix is the precomputed principal stress or trial principal stress
        PrincipalStress(i) = rStressMatrix(i,i);
    }

    // Sorting Principal Stress and Strain - "0" is the largest one and "2" is the lowest one
    MPMStressPrincipalInvariantsUtility::SortPrincipalStress(PrincipalStress, MainStrain, rReturnMappingVariables.MainDirections);

    // Assigning to local variables
    mElasticPrincipalStrain = MainStrain;

    // Check for the yield Condition -- calling the yield criterion
    rReturnMappingVariables.TrialStateFunction = 0.0;
    rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PrincipalStress, 0.0, mMaterialParameters.PreconsolidationPressure);
    
    // If yield is reached, do return mapping
    if (rReturnMappingVariables.TrialStateFunction <= 0.0)
    {
        mRegion = 0;
        mPrincipalStressUpdated = PrincipalStress;
        PlasticityActive = false;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);

        this->UpdateStateVariables(mPrincipalStressUpdated);
        
    }
    else
    {
        unsigned int Region = 0;
        Vector PrincipalStressUpdated = ZeroVector(3);

        // Perform return mapping to the yield surface: Will update mElasticPrincipalStrain, Region, and PrincipalStressUpdated
        bool converged = this->CalculateConsistencyCondition(rReturnMappingVariables, PrincipalStress, mElasticPrincipalStrain, Region, PrincipalStressUpdated);
        KRATOS_ERROR_IF(!converged) << "Warning:: Constitutive Law does not converge! "<<std::endl;

        mRegion = Region;
        mPrincipalStressUpdated = PrincipalStressUpdated;

        PlasticityActive = true;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,true);
    }

    // rStressMatrix is the matrix of the updated stress in cartesian configuration -- this function perform back transformation
    this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, mPrincipalStressUpdated, rStressMatrix);
 
    // Delta plastic strain
    mPlasticPrincipalStrain = MainStrain - mElasticPrincipalStrain ;

    // We're saving the updated info in terms of principal strain and stress in these matrix
    // these information will be used for the evaluation of the second contribution in the
    // consistent tangent matrix
    for (unsigned int i=0; i<3; i++)
    {
        rReturnMappingVariables.StrainMatrix(i,i) = mElasticPrincipalStrain(i);
        rReturnMappingVariables.TrialIsoStressMatrix(i,i) = mPrincipalStressUpdated(i);
    }

    rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED,true);

    return PlasticityActive;

}

bool BorjaCamClayPlasticFlowRule::CalculateConsistencyCondition(RadialReturnVariables& rReturnMappingVariables, Vector& rPrincipalStress, 
    Vector& rPrincipalStrain, unsigned int& region, Vector& rPrincipalStressUpdated)
{

    // Calculate stress return in principal stress
    // The flow rule is written for associated plasticity assumption
    // Refer to paper by (Borja, 1998) for the theoretical description
    bool converged = false;

    // Initiate Main Matrices and Vectors for Newton Iteration
    Vector UnknownVector, DeltaUnknownVector;
    Vector RHSVector = ZeroVector(3);
    Matrix LHSMatrix = ZeroMatrix(3,3);
    Matrix InverseLHSMatrix = ZeroMatrix(3,3);
    
    // Initiate iterator variable
    unsigned int counter = 0;
    unsigned int maxcounter = 20;
    double InitialNormResidual;
    const double tolerance = 5e-03;
    const double norm_tolerance = 5e-012;

    // Initial UnknownVector
    UnknownVector = ZeroVector(3);
    DeltaUnknownVector = ZeroVector(3);
    double TrialVolumetricStrain, TrialDeviatoricStrain;
    Vector TrialDeviatoricStrainVector;
    this->CalculateStrainInvariantsFromPrincipalStrain(rPrincipalStrain, TrialVolumetricStrain, TrialDeviatoricStrain, TrialDeviatoricStrainVector);
    UnknownVector(0) = TrialVolumetricStrain;
    UnknownVector(1) = TrialDeviatoricStrain;

    // Constant direction of deviatoric strain vector -- with check if TrialDeviatoricStrain == 0
    Vector DirectionStrainVector = ZeroVector(3);
    if(std::abs(TrialDeviatoricStrain) > 1.e-9) 
        DirectionStrainVector = std::sqrt(2.0/3.0) * TrialDeviatoricStrainVector / TrialDeviatoricStrain;

    // Initialize additional temporary Matrices and Vectors;
    mStateFunction = 0.0;
    mStateFunctionFirstDerivative  = ZeroVector(3);
    mStateFunctionSecondDerivative = ZeroVector(6);
    Vector PrincipalStressVector = rPrincipalStress;
    Vector PrincipalStrainVector = ZeroVector(3);

    // Initialize Material parameters
    const double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    const double OtherSlope    = mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];

    // Begin Newton Iteration
    while(!converged)
    {
        // Calculate state function and derivatives
        double rAlpha  = TrialVolumetricStrain - UnknownVector(0);
        mStateFunction = mpYieldCriterion->CalculateYieldCondition(mStateFunction, PrincipalStressVector, rAlpha, mMaterialParameters.PreconsolidationPressure);
        mpYieldCriterion->CalculateYieldFunctionDerivative(PrincipalStressVector, mStateFunctionFirstDerivative, rAlpha, mMaterialParameters.PreconsolidationPressure);
        mpYieldCriterion->CalculateYieldFunctionSecondDerivative(PrincipalStressVector, mStateFunctionSecondDerivative);

        // Calculate plastic hardening modulus K_p
        double K_p = mpYieldCriterion->GetHardeningLaw().CalculateHardening(K_p, rAlpha, mMaterialParameters.PreconsolidationPressure);
        K_p *= 1.0 / (OtherSlope-SwellingSlope);

        // Calculate RHS Vector
        RHSVector(0) = UnknownVector(0) - TrialVolumetricStrain + UnknownVector(3) * mStateFunctionFirstDerivative(0);
        RHSVector(1) = UnknownVector(1) - TrialDeviatoricStrain + UnknownVector(3) * mStateFunctionFirstDerivative(1);;
        RHSVector(2) = mStateFunction;

        // Calculate RHS Norm (Residual Norm)
        if (counter == 0) InitialNormResidual = norm_2(RHSVector);
        double CurrentNormResidual = norm_2(RHSVector);
        double NormRatio = CurrentNormResidual/InitialNormResidual;

        // Calculate LHS Matrix
        this->CalculateLHSMatrix(LHSMatrix, PrincipalStressVector, UnknownVector, K_p);

        // Compute Inverse LHS Matrix
        double detLHS = MathUtils<double>::Det(LHSMatrix);
        MathUtils<double>::InvertMatrix( LHSMatrix, InverseLHSMatrix, detLHS);

        // Update DeltaUnknownVector
        DeltaUnknownVector = prod(InverseLHSMatrix, RHSVector);

        // Update Unknown Vector: x^(k+1) = x^(k) + dx^(k+1)
        counter += 1;
        UnknownVector += DeltaUnknownVector;

        // Update PrincipalStrainVector from UnknownVector
        this->CalculatePrincipalStrainFromStrainInvariants(PrincipalStrainVector, UnknownVector(0), UnknownVector(1), DirectionStrainVector);

        // Update PrincipalStressVector from the new UnknownVector
        this->CalculatePrincipalStressVector(PrincipalStrainVector, PrincipalStressVector);

        // Weighted residual Convergence criteria - to exit Newton iteration loop
        if( std::abs(mStateFunction) <= tolerance || std::abs(NormRatio) <= norm_tolerance || counter == maxcounter)
        {    
            // These updates are done since the following variables will be used during ComputeElastoPlasticTangentMatrix
                rAlpha  = TrialVolumetricStrain - UnknownVector(0);
            this->UpdateStateVariables(PrincipalStressVector, rAlpha, UnknownVector(2));
                
            // Update rPrincipalStressUpdated and rPrincipalStrain as the final elastic principal stress and strain
            rPrincipalStressUpdated = PrincipalStressVector;
            rPrincipalStrain = PrincipalStrainVector;

            region = 1;
            converged = true;
        }
    }

    return converged;
}

void BorjaCamClayPlasticFlowRule::CalculateLHSMatrix(Matrix& rLHSMatrix, const Vector& rPrincipalStressVector, const Vector& rUnknownVector, const double& rK_p)
{
    // Reset Zero
    rLHSMatrix = ZeroMatrix(3,3);

    // Compute ElasticMatrix D^e
    Matrix ElasticMatrixDe = ZeroMatrix(2,2);
    this->ComputeElasticMatrix_2X2(rPrincipalStressVector, rUnknownVector(0), rUnknownVector(1), ElasticMatrixDe);

    // Compute Hessian Matrix H
    Matrix HessianMatrixH = ZeroMatrix(2,2);
    this->CalculateHessianMatrix_2x2(HessianMatrixH);

    // Compute GMatrix
    Matrix GMatrix = prod(HessianMatrixH, ElasticMatrixDe); 

    // Arrange rLHSMatrix
    rLHSMatrix(0,0) = -1.0 * (1.0 + rUnknownVector(2) * (GMatrix(0,0) + rK_p * mStateFunctionSecondDerivative(5)));
    rLHSMatrix(0,1) = -1.0 * (rUnknownVector(2) * GMatrix(0,1));
    rLHSMatrix(0,2) = -1.0 * (mStateFunctionFirstDerivative(0));

    rLHSMatrix(1,0) = -1.0 * (rUnknownVector(2) * (GMatrix(1,0) + rK_p * mStateFunctionSecondDerivative(4)));
    rLHSMatrix(1,1) = -1.0 * (1.0 + rUnknownVector(2) * GMatrix(1,1));
    rLHSMatrix(1,2) = -1.0 * (mStateFunctionFirstDerivative(1));

    rLHSMatrix(2,0) = -1.0 * (ElasticMatrixDe(0,0) * mStateFunctionFirstDerivative(0) + ElasticMatrixDe(1,0) * mStateFunctionFirstDerivative(1) + rK_p * mStateFunctionFirstDerivative(2));
    rLHSMatrix(2,1) = -1.0 * (ElasticMatrixDe(0,1) * mStateFunctionFirstDerivative(0) + ElasticMatrixDe(1,1) * mStateFunctionFirstDerivative(1));
    rLHSMatrix(2,2) = 0.0;

}

void BorjaCamClayPlasticFlowRule::CalculateHessianMatrix_2x2(Matrix& rHessianMatrix)
{
    // Material parameters
    const double ShearM = mpYieldCriterion->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];

    // Assemble matrix
    rHessianMatrix(0,0) = 2.0;
    rHessianMatrix(1,1) = 2.0 / pow(ShearM, 2.0);
    rHessianMatrix(0,1) = 0.0;
    rHessianMatrix(1,0) = 0.0;

}

void BorjaCamClayPlasticFlowRule::ComputeElasticMatrix_2X2(const Vector& rPrincipalStressVector, const double& rVolumetricStrain, const double& rDeviatoricStrain, Matrix& rElasticMatrix)
{
    // Material parameters
    const double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    const double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];

    double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    const double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ReferencePressure /= OCR;

    const double ConstantShearModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[INITIAL_SHEAR_MODULUS];
    const double ShearModulus = AlphaShear * ReferencePressure * std::exp( -(rVolumetricStrain - mInitialVolumetricStrain) / SwellingSlope);

    // Decompose principalstress
    double MeanStressP, DeviatoricQ;
    MPMStressPrincipalInvariantsUtility::CalculateStressInvariants( rPrincipalStressVector, MeanStressP, DeviatoricQ);
    DeviatoricQ *= std::sqrt(3.0); //Q = sqrt(3) * J2

    // Assemble matrix
    rElasticMatrix(0,0) = - MeanStressP / SwellingSlope;
    rElasticMatrix(1,1) = 3.0 * (ConstantShearModulus - ShearModulus);
    rElasticMatrix(0,1) = 3.0 * ShearModulus * rDeviatoricStrain / SwellingSlope;
    rElasticMatrix(1,0) = rElasticMatrix(0,1);
}

void BorjaCamClayPlasticFlowRule::ComputePlasticMatrix_2X2(const Vector& rPrincipalStressVector, const double& rVolumetricStrain, const double& rDeviatoricStrain, const Matrix& rElasticMatrix, Matrix& rPlasticMatrix)
{
    // Initialize used temporary matrices and vectors
    Vector a = ZeroVector(2);
    Matrix b = ZeroMatrix(2,2);
    Vector c = ZeroVector(2);
    Vector d = ZeroVector(2);
    double e;

    // Compute Hessian Matrix H
    Matrix HessianMatrixH = ZeroMatrix(2,2);
    this->CalculateHessianMatrix_2x2(HessianMatrixH);

    // Compute GMatrix
    Matrix GMatrix = prod(HessianMatrixH, rElasticMatrix); 

    // Initiate coefficient parameters
    const double K_p      = mMaterialParameters.PlasticHardeningModulus;
    const double K_ptrial = -K_p;
    const double DeltaPhi = mMaterialParameters.ConsistencyParameter;

    // Construct Matrix b
    b(0,0) = 1.0 + DeltaPhi * ( GMatrix(0,0) + K_p * mStateFunctionSecondDerivative(5) );
    b(0,1) = DeltaPhi * GMatrix(0,1);
    b(1,0) = DeltaPhi * ( GMatrix(1,0) + K_p * mStateFunctionSecondDerivative(4) );
    b(1,1) = 1.0 + DeltaPhi * GMatrix(1,1);
    double detb = MathUtils<double>::Det(b);

    // Construct Vector c
    c(0) = 1.0 - DeltaPhi * K_ptrial * mStateFunctionSecondDerivative(5);
    c(1) = - DeltaPhi * K_ptrial * mStateFunctionSecondDerivative(4);

    // Construct Vector d
    d(0) = rElasticMatrix(0,0) * mStateFunctionFirstDerivative(0) + rElasticMatrix(1,0) * mStateFunctionFirstDerivative(1) + K_p * mStateFunctionFirstDerivative(2);
    d(1) = rElasticMatrix(0,1) * mStateFunctionFirstDerivative(0) + rElasticMatrix(1,1) * mStateFunctionFirstDerivative(1);

    // Construct Coefficient e
    e  = d(0) * ( b(1,1) * mStateFunctionFirstDerivative(0) - b(0,1) * mStateFunctionFirstDerivative(1) ) 
       + d(1) * ( b(0,0) * mStateFunctionFirstDerivative(1) - b(1,0) * mStateFunctionFirstDerivative(0) );

    // Construct Vector a
    a(0) = ( d(0) * (b(1,1) * c(0) - b(0,1) * c(1)) + d(1) * (b(0,0) * c(1) - b(1,0) * c(0)) + detb * K_ptrial * mStateFunctionFirstDerivative(2) );
    a(1) = std::sqrt(2.0/3.0) * ( d(1) * b(0,0) - d(0) * b(0,1) );
    
    // Check if e == 0
    if (std::abs(e) < 1.e-9){ a *= 1.0 / 1.e-9; }
    else{ a *= 1.0 / e; }

    // Arrange rPlasticMatrix from all the constructed variables
    rPlasticMatrix(0,0) = b(1,1) * (c(0) - a(0) * mStateFunctionFirstDerivative(0)) - b(0,1) * (c(1) - a(0) * mStateFunctionFirstDerivative(1));
    rPlasticMatrix(0,1) = b(0,1) * (-1.0 + std::sqrt(3.0/2.0) * a(1) * mStateFunctionFirstDerivative(1)) - std::sqrt(3.0/2.0) * b(1,1) * a(1) * mStateFunctionFirstDerivative(0);
    rPlasticMatrix(1,0) = b(0,0) * (c(1) - a(0) * mStateFunctionFirstDerivative(1)) - b(1,0) * (c(0) - a(0) * mStateFunctionFirstDerivative(0));
    rPlasticMatrix(1,1) = b(0,0) * (1.0  - std::sqrt(3.0/2.0) * a(1) * mStateFunctionFirstDerivative(1)) + std::sqrt(3.0/2.0) * b(1,0) * a(1) * mStateFunctionFirstDerivative(0);
    rPlasticMatrix     *= 1.0 / detb; 

}

// Compute Trial elastic principal stress matrix from Trial elastic principal strain matrix
void BorjaCamClayPlasticFlowRule::CalculatePrincipalStressTrial(const RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix)
{
    Vector MainStrain      = ZeroVector(3);
    Vector PrincipalStress = ZeroVector(3);

    for (unsigned int i = 0; i<3; ++i)
    {
        MainStrain[i] = rNewElasticLeftCauchyGreen(i,i);
    }

    this->CalculatePrincipalStressVector(MainStrain, PrincipalStress);

    // Evalute the Kirchhoff principal stress
    for(unsigned int i=0; i<3; i++)
    {
        rStressMatrix(i,i) = PrincipalStress(i);
    }

}


// Function to compute Principal Stress Vector from Principal Strain Vector
void BorjaCamClayPlasticFlowRule::CalculatePrincipalStressVector(Vector& rPrincipalStrain, Vector& rPrincipalStress)
{
    // Calculate volumetric and deviatoric strains from princial strain
    double VolumetricStrain, DeviatoricStrain;
    Vector DeviatoricStrainVector;
    this->CalculateStrainInvariantsFromPrincipalStrain(rPrincipalStrain, VolumetricStrain, DeviatoricStrain, DeviatoricStrainVector);

    // Calculate MeanStressP and DeviatoricStressQ
    double MeanStressP;
    this->CalculateMeanStress(VolumetricStrain, DeviatoricStrain, MeanStressP);
    this->CalculateDeviatoricStress(VolumetricStrain, DeviatoricStrainVector, rPrincipalStress);

    // Combine into PrincipalStress
    for (unsigned int i = 0; i<3; ++i)
        rPrincipalStress[i] += MeanStressP;

}


// Function calculate Mean Stress from given VolumetricStrain and DeviatoricStrain
void BorjaCamClayPlasticFlowRule::CalculateMeanStress(const double& rVolumetricStrain, const double& rDeviatoricStrain, double& rMeanStress)
{
    const double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    const double AlphaShear    = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];

    double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    const double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ReferencePressure /= OCR;

    rMeanStress = ReferencePressure * std::exp( -(rVolumetricStrain - mInitialVolumetricStrain) / SwellingSlope) * (1.0 + 1.5 * AlphaShear * pow(rDeviatoricStrain, 2) / SwellingSlope);

}


// Function calculate DeviatoricStressVector from given VolumetricStrain and DeviatoricStrainVector
void BorjaCamClayPlasticFlowRule::CalculateDeviatoricStress(const double& rVolumetricStrain, const Vector & rDeviatoricStrainVector, Vector& rDeviatoricStress)
{
    double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    const double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ReferencePressure /= OCR;    
    
    const double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    const double AlphaShear    = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];
    const double ConstantShearModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[INITIAL_SHEAR_MODULUS];

    rDeviatoricStress = rDeviatoricStrainVector;
    const double ShearModulus = AlphaShear * ReferencePressure * std::exp( -(rVolumetricStrain - mInitialVolumetricStrain) / SwellingSlope);
    rDeviatoricStress *= 2.0 * ( ConstantShearModulus - ShearModulus );

}

// Function which returns principal strains from volumetric and deviatoric strain components
void BorjaCamClayPlasticFlowRule::CalculatePrincipalStrainFromStrainInvariants(Vector& rPrincipalStrain, const double& rVolumetricStrain, const double& rDeviatoricStrain, const Vector& rDirectionVector)
{
    rPrincipalStrain = ZeroVector(3);
    
    for (unsigned int i = 0; i<3; ++i)
    {
        rPrincipalStrain[i] += 1.0/3.0 * rVolumetricStrain;
    }
    rPrincipalStrain += std::sqrt(3.0/2.0) * rDeviatoricStrain * rDirectionVector;
}


// Function which returns volumetric and deviatoric strain components from principal strain
void BorjaCamClayPlasticFlowRule::CalculateStrainInvariantsFromPrincipalStrain(const Vector& rPrincipalStrain, double& rVolumetricStrain, double& rDeviatoricStrain, Vector& rDeviatoricStrainVector)
{
    rDeviatoricStrainVector = rPrincipalStrain;
    rVolumetricStrain = sum(rPrincipalStrain);
    for (unsigned int i = 0; i<3; ++i)
    {
        rDeviatoricStrainVector[i] -= 1.0/3.0 * rVolumetricStrain;
    }
    rDeviatoricStrain = std::sqrt(2.0/3.0) * norm_2(rDeviatoricStrainVector);
}

// Function to return matrix from principal space to normal space
void BorjaCamClayPlasticFlowRule::ReturnStressFromPrincipalAxis(const Matrix& rEigenVectors, const Vector& rPrincipalStress, Matrix& rStressMatrix)
{
    rStressMatrix = ZeroMatrix(3,3); 
    Vector auxN = ZeroVector(3);
    Matrix auxM = ZeroMatrix(3,3);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
            auxN(j) = rEigenVectors(j,i);
        auxM = MathUtils<double>::TensorProduct3(auxN, auxN);
        rStressMatrix += rPrincipalStress(i)*auxM;
    }
}

// Function that compute the consistent tangent stiffness matrix (in normal space) considering both elastic and elasto-plastic case
void BorjaCamClayPlasticFlowRule::ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& alfa, Matrix& rConsistMatrix)
{
    // Prepare PrincipalStressVector and its invariants
    Vector PrincipalStressVector = mPrincipalStressUpdated;
    double MeanStressP, DeviatoricQ;
    MPMStressPrincipalInvariantsUtility::CalculateStressInvariants( PrincipalStressVector, MeanStressP, DeviatoricQ);
    DeviatoricQ *= std::sqrt(3.0); //Q = sqrt(3) * J2

    // Compute StrainComponents and direction -- with check if DeviatoricStrain == 0
    double VolumetricStrain, DeviatoricStrain;
    Vector DeviatoricStrainVector;
    this->CalculateStrainInvariantsFromPrincipalStrain(mElasticPrincipalStrain, VolumetricStrain, DeviatoricStrain, DeviatoricStrainVector);
    
    Vector DirectionVector = ZeroVector(3);
    if (std::abs(DeviatoricStrain) > 1.e-9)
        DirectionVector = std::sqrt(2.0/3.0) * DeviatoricStrainVector / DeviatoricStrain;

    // Compute ElasticMatrix (2x2) D^e
    Matrix ElasticMatrixDe = ZeroMatrix(2,2);
    this->ComputeElasticMatrix_2X2(PrincipalStressVector, VolumetricStrain, DeviatoricStrain, ElasticMatrixDe);

    // Compute PlasticMatrix (2x2) D^p
    Matrix PlasticMatrixDp = IdentityMatrix(2);
    if (rReturnMappingVariables.Options.Is(MPMFlowRule::PLASTIC_REGION))
    {
        this->ComputePlasticMatrix_2X2(PrincipalStressVector, VolumetricStrain, DeviatoricStrain, ElasticMatrixDe, PlasticMatrixDp);
    }

    // Compute Elasto-plastic (2x2) D^ep = D^e D^p
    Matrix ElastoPlasticMatrixDep = prod(ElasticMatrixDe, PlasticMatrixDp);

    // Prepare FourthOrderIdentity and IdentityCross
    Matrix FourthOrderIdentity = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i)
        FourthOrderIdentity(i,i) = 1.0;

    for (unsigned int i = 3; i<6; ++i)
        FourthOrderIdentity(i,i) = 0.50;

    Matrix IdentityCross = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i) 
    {
        for (unsigned int j = 0; j<3; ++j) 
        {
            IdentityCross(i,j) = 1.0;
        }
    }

    // Prepare Tensor_NxN, Tensor_1xN, and Tensor_Nx1
    Matrix Tensor_NxN = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i) 
    {
        for (unsigned int j = 0; j<3; ++j) 
        {
            Tensor_NxN(i,j) = DirectionVector(i) * DirectionVector(j);
        }
    }

    Matrix Tensor_1xN = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i) 
    {
        for (unsigned int j = 0; j<3; ++j) 
        {
            Tensor_1xN(i,j) = DirectionVector(j);
        }
    }

    Matrix Tensor_Nx1 = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i) 
    {
        for (unsigned int j = 0; j<3; ++j) 
        {
            Tensor_Nx1(i,j) = DirectionVector(i);
        }
    } 

    // Perform check in case DeviatoricStrain == 0
    double DeviatoricQ_by_DeviatoricStrain; 
    if (std::abs(DeviatoricStrain) < 1.e-9 ) {DeviatoricQ_by_DeviatoricStrain = DeviatoricQ / 1.e-9;}
    else{DeviatoricQ_by_DeviatoricStrain = DeviatoricQ / DeviatoricStrain;}  

    // Compute Consistent Tangent Stiffness matrix in principal space
    Matrix DepcP = ZeroMatrix(6,6);
    DepcP  = ( ElastoPlasticMatrixDep(0,0) - 2.0 * DeviatoricQ_by_DeviatoricStrain / 9.0 ) * IdentityCross;
    DepcP += ( std::sqrt(2.0/3.0) * ElastoPlasticMatrixDep(0,1) ) * Tensor_1xN; 
    DepcP += ( std::sqrt(2.0/3.0) * ElastoPlasticMatrixDep(1,0) ) * Tensor_Nx1;
    DepcP += ( 2.0 * DeviatoricQ_by_DeviatoricStrain / 3.0 ) * (FourthOrderIdentity - Tensor_NxN);
    DepcP += ( 2.0 / 3.0 * ElastoPlasticMatrixDep(1,1) ) * Tensor_NxN;

    // Return constitutive matrix from principal space to normal space
    Matrix A = ZeroMatrix(6,6);
    Matrix ATrans = ZeroMatrix(6,6); 
    this->CalculateTransformationMatrix(rReturnMappingVariables.MainDirections, A);
    ATrans = trans(A);

    Matrix AuxMat = ZeroMatrix(6,6);
    AuxMat = prod(ATrans, DepcP);
    rConsistMatrix = prod(AuxMat, A);

}

void BorjaCamClayPlasticFlowRule::CalculateTransformationMatrix(const Matrix& rMainDirection, Matrix& rA)
{
    Matrix A1 = ZeroMatrix(3,3);
    Matrix A2 = ZeroMatrix(3,3);
    Matrix A3 = ZeroMatrix(3,3);
    Matrix A4 = ZeroMatrix(3,3);
    for (unsigned int i = 0; i<3 ; i++)
    {
        for(unsigned int j = 0; j<3 ; j++)
        {
            A1(i,j) = rMainDirection(i,j) * rMainDirection(i,j);
            rA(i,j) = A1(i,j);
        }
    }
    Vector Hj1 = ZeroVector(3);
    Hj1(0) = 0;
    Hj1(1) = 2;
    Hj1(2) = 1;

    Vector Hj2 = ZeroVector(3);
    Hj2(0) = 1;
    Hj2(1) = 0;
    Hj2(2) = 2;

    for(unsigned int k = 0; k<3; k++)
    {
        for(unsigned int l = 0; l<3; l++)
        {
            A2(k,l) = rMainDirection(k,Hj1(l)) * rMainDirection(k,Hj2(l));
            A3(k,l) = rMainDirection(Hj1(k),l) * rMainDirection(Hj2(k),l);
            A4(k,l) = rMainDirection(Hj1(k),Hj1(l)) * rMainDirection(Hj2(k),Hj2(l)) + rMainDirection(Hj2(k),Hj1(l)) * rMainDirection(Hj1(k),Hj2(l));
        }
    }

    for(unsigned int i = 0 ; i<3 ; i++)
    {
        for(unsigned int j = 3; j<6 ; j++)
        {
            unsigned int index_j = j - 3;
            unsigned int index_i = i + 3;

            rA(i,j) = A2(i, index_j);
            rA(index_i, index_j) = A3(i, index_j);
            rA(index_i, j) = A4(i, index_j);
        }
    }

    rA = trans(rA);
}


// Function that gives the elastic left cauchy green tensor B
Matrix BorjaCamClayPlasticFlowRule::GetElasticLeftCauchyGreen(RadialReturnVariables& rReturnMappingVariables)
{
    Vector Landa2 = ZeroVector(3);

    for (unsigned int i = 0; i<3; ++i)
        Landa2(i) = std::exp(2.0*mElasticPrincipalStrain(i));

    Matrix OutPut = ZeroMatrix(3,3);
    this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, Landa2, OutPut);

    return OutPut;
}

// Function that updates internal variable at every time step once the nonlinear iteration converges
bool BorjaCamClayPlasticFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{
    // Compute Delta Plastic Strain
    double NormPlasticPrincipalStrain = norm_2(mPlasticPrincipalStrain);

    // Compute Strain Components and its invariants
    double VolumetricStrain, DeviatoricStrain;
    Vector DeviatoricStrainVector;
    this->CalculateStrainInvariantsFromPrincipalStrain(mPlasticPrincipalStrain, VolumetricStrain, DeviatoricStrain, DeviatoricStrainVector);

    // Update Equivalent Plastic Strain
    mInternalVariables.DeltaPlasticStrain = NormPlasticPrincipalStrain;
    mInternalVariables.EquivalentPlasticStrain += NormPlasticPrincipalStrain;

    // Update Accumulated Plastic Volumetric Strain
    mInternalVariables.DeltaPlasticVolumetricStrain = VolumetricStrain;
    mInternalVariables.AccumulatedPlasticVolumetricStrain += VolumetricStrain;

    // Update Accumulated Plastic Deviatoric Strain
    mInternalVariables.DeltaPlasticDeviatoricStrain = DeviatoricStrain;
    mInternalVariables.AccumulatedPlasticDeviatoricStrain += DeviatoricStrain;

    // Update Preconsolidation Stress for the next time step
    double newPreconsolidationStress = mpYieldCriterion->GetHardeningLaw().CalculateHardening(newPreconsolidationStress, VolumetricStrain, mMaterialParameters.PreconsolidationPressure);
    mMaterialParameters.PreconsolidationPressure = newPreconsolidationStress;

    return true;
}

void BorjaCamClayPlasticFlowRule::UpdateStateVariables(const Vector rPrincipalStress, const double rAlpha, const double rConsistencyParameter)
{
    // Calculate final state function and derivatives
    mStateFunction = mpYieldCriterion->CalculateYieldCondition(mStateFunction, rPrincipalStress, rAlpha, mMaterialParameters.PreconsolidationPressure);
    mpYieldCriterion->CalculateYieldFunctionDerivative(rPrincipalStress, mStateFunctionFirstDerivative, rAlpha, mMaterialParameters.PreconsolidationPressure);
    mpYieldCriterion->CalculateYieldFunctionSecondDerivative(rPrincipalStress, mStateFunctionSecondDerivative);
    
    // Update plastic hardening modulus K_p
    const double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    const double OtherSlope    = mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];
    double K_p = mpYieldCriterion->GetHardeningLaw().CalculateHardening(K_p, rAlpha, mMaterialParameters.PreconsolidationPressure);
    K_p *= 1.0 / (OtherSlope-SwellingSlope);
    mMaterialParameters.PlasticHardeningModulus = K_p;

    // Update Consistency Parameter
    mMaterialParameters.ConsistencyParameter = rConsistencyParameter;
}

double BorjaCamClayPlasticFlowRule::GetPI()
{
    return std::atan(1.0)*4.0;
}

unsigned int BorjaCamClayPlasticFlowRule::GetPlasticRegion()
{
    return mRegion;
}

void BorjaCamClayPlasticFlowRule::save( Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMFlowRule )
}

void BorjaCamClayPlasticFlowRule::load( Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMFlowRule )

}

} //end namespace kratos
