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
    const double swelling_slope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    const double other_slope    = mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];

    mMaterialParameters.PreconsolidationPressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    mMaterialParameters.PlasticHardeningModulus  = mMaterialParameters.PreconsolidationPressure/ (other_slope-swelling_slope);
    mMaterialParameters.ConsistencyParameter     = 0.0;
}


bool BorjaCamClayPlasticFlowRule::CalculateReturnMapping(
            RadialReturnVariables& rReturnMappingVariables,
            const Matrix& rIncrementalDeformationGradient,
            Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{
    bool plasticity_active = false;
    rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);

    Vector principal_stress = ZeroVector(3);
    Vector main_strain      = ZeroVector(3);

    for (unsigned int i = 0; i<3; ++i)
        main_strain[i] = rNewElasticLeftCauchyGreen(i,i);

    for(unsigned int i=0; i<3; i++)
    {
        // the rStressMatrix is the precomputed principal stress or trial principal stress
        principal_stress[i] = rStressMatrix(i,i);
    }

    // Sorting Principal Stress and Strain - "0" is the largest one and "2" is the lowest one
    MPMStressPrincipalInvariantsUtility::SortPrincipalStress(principal_stress, main_strain, rReturnMappingVariables.MainDirections);

    // Assigning to local variables
    mElasticPrincipalStrain = main_strain;

    // Check for the yield Condition -- calling the yield criterion
    rReturnMappingVariables.TrialStateFunction = 0.0;
    rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, principal_stress, 0.0, mMaterialParameters.PreconsolidationPressure);

    // If yield is reached, do return mapping
    if (rReturnMappingVariables.TrialStateFunction <= 0.0)
    {
        mRegion = 0;
        mPrincipalStressUpdated = principal_stress;
        plasticity_active = false;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);

        this->UpdateStateVariables(mPrincipalStressUpdated);

    }
    else
    {
        unsigned int region = 0;
        BoundedVector<double,3> principal_stress_updated = ZeroVector(3);

        // Perform return mapping to the yield surface: Will update mElasticPrincipalStrain, Region, and principal_stress_updated
        bool converged = this->CalculateConsistencyCondition(rReturnMappingVariables, principal_stress, mElasticPrincipalStrain, region, principal_stress_updated);
        KRATOS_ERROR_IF(!converged) << "Warning:: Constitutive Law does not converge! "<<std::endl;

        mRegion = region;
        mPrincipalStressUpdated = principal_stress_updated;

        plasticity_active = true;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,true);
    }

    // rStressMatrix is the matrix of the updated stress in cartesian configuration -- this function perform back transformation
    this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, mPrincipalStressUpdated, rStressMatrix);

    // Delta plastic strain
    mPlasticPrincipalStrain = main_strain - mElasticPrincipalStrain ;

    // We're saving the updated info in terms of principal strain and stress in these matrix
    // these information will be used for the evaluation of the second contribution in the
    // consistent tangent matrix
    for (unsigned int i=0; i<3; i++)
    {
        rReturnMappingVariables.StrainMatrix(i,i) = mElasticPrincipalStrain[i];
        rReturnMappingVariables.TrialIsoStressMatrix(i,i) = mPrincipalStressUpdated[i];
    }

    rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED,true);

    return plasticity_active;

}

bool BorjaCamClayPlasticFlowRule::CalculateConsistencyCondition(RadialReturnVariables& rReturnMappingVariables, const BoundedVector<double,3>& rPrincipalStress,
    BoundedVector<double,3>& rPrincipalStrain, unsigned int& region, BoundedVector<double,3>& rPrincipalStressUpdated)
{

    // Calculate stress return in principal stress
    // The flow rule is written for associated plasticity assumption
    // Refer to paper by (Borja, 1998) for the theoretical description
    bool converged = false;

    // Initiate Main Matrices and Vectors for Newton Iteration
    BoundedVector<double,3> unknown_vector, delta_unknown_vector;
    BoundedVector<double,3> rhs_vector = ZeroVector(3);
    Matrix lhs_matrix = ZeroMatrix(3,3);
    Matrix inv_lhs_matrix = ZeroMatrix(3,3);

    // Initiate iterator variable
    unsigned int counter = 0;
    unsigned int maxcounter = 20;
    double initial_norm_residual;
    const double tolerance = 5e-03;
    const double norm_tolerance = 5e-012;

    // Initial unknown_vector
    unknown_vector = ZeroVector(3);
    delta_unknown_vector = ZeroVector(3);
    double trial_volumetric_strain, trial_deviatoric_strain;
    BoundedVector<double,3> trial_deviatoric_strain_vector;
    this->CalculateStrainInvariantsFromPrincipalStrain(rPrincipalStrain, trial_volumetric_strain, trial_deviatoric_strain, trial_deviatoric_strain_vector);
    unknown_vector[0] = trial_volumetric_strain;
    unknown_vector[1] = trial_deviatoric_strain;

    // Constant direction of deviatoric strain vector -- with check if trial_deviatoric_strain == 0
    BoundedVector<double,3> direction_strain_vector = ZeroVector(3);
    if(std::abs(trial_deviatoric_strain) > 1.e-9)
        direction_strain_vector = std::sqrt(2.0/3.0) * trial_deviatoric_strain_vector / trial_deviatoric_strain;

    // Initialize additional temporary Matrices and Vectors;
    mStateFunction = 0.0;
    mStateFunctionFirstDerivative  = ZeroVector(3);
    mStateFunctionSecondDerivative = ZeroVector(6);
    BoundedVector<double,3> principal_stress_vector = rPrincipalStress;
    BoundedVector<double,3> principal_strain_vector = ZeroVector(3);

    // Initialize Material parameters
    const double swelling_slope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    const double other_slope    = mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];

    // Begin Newton Iteration
    while(!converged)
    {
        // Calculate state function and derivatives
        double alpha  = trial_volumetric_strain - unknown_vector[0];
        mStateFunction = mpYieldCriterion->CalculateYieldCondition(mStateFunction, principal_stress_vector, alpha, mMaterialParameters.PreconsolidationPressure);
        mpYieldCriterion->CalculateYieldFunctionDerivative(principal_stress_vector, mStateFunctionFirstDerivative, alpha, mMaterialParameters.PreconsolidationPressure);
        mpYieldCriterion->CalculateYieldFunctionSecondDerivative(principal_stress_vector, mStateFunctionSecondDerivative);

        // Calculate plastic hardening modulus k_p
        double k_p = mpYieldCriterion->GetHardeningLaw().CalculateHardening(k_p, alpha, mMaterialParameters.PreconsolidationPressure);
        k_p *= 1.0 / (other_slope-swelling_slope);

        // Calculate RHS Vector
        rhs_vector[0] = unknown_vector[0] - trial_volumetric_strain + unknown_vector[2] * mStateFunctionFirstDerivative[0];
        rhs_vector[1] = unknown_vector[1] - trial_deviatoric_strain + unknown_vector[2] * mStateFunctionFirstDerivative[1];
        rhs_vector[2] = mStateFunction;

        // Calculate RHS Norm (Residual Norm)
        if (counter == 0) initial_norm_residual = norm_2(rhs_vector);
        double current_norm_residual = norm_2(rhs_vector);
        double norm_ratio = current_norm_residual/initial_norm_residual;

        // Calculate LHS Matrix
        this->CalculateLHSMatrix(lhs_matrix, principal_stress_vector, unknown_vector, k_p);

        // Compute Inverse LHS Matrix
        double det_lhs = MathUtils<double>::Det(lhs_matrix);
        MathUtils<double>::InvertMatrix( lhs_matrix, inv_lhs_matrix, det_lhs);

        // Update delta_unknown_vector
        delta_unknown_vector = prod(inv_lhs_matrix, rhs_vector);

        // Update Unknown Vector: x^(k+1) = x^(k) + dx^(k+1)
        counter += 1;
        unknown_vector += delta_unknown_vector;

        // Update principal_strain_vector from unknown_vector
        this->CalculatePrincipalStrainFromStrainInvariants(principal_strain_vector, unknown_vector[0], unknown_vector[1], direction_strain_vector);

        // Update principal_stress_vector from the new unknown_vector
        this->CalculatePrincipalStressVector(principal_strain_vector, principal_stress_vector);

        // Weighted residual Convergence criteria - to exit Newton iteration loop
        if( std::abs(mStateFunction) <= tolerance || std::abs(norm_ratio) <= norm_tolerance || counter == maxcounter)
        {
            // These updates are done since the following variables will be used during ComputeElastoPlasticTangentMatrix
            alpha  = trial_volumetric_strain - unknown_vector[0];
            this->UpdateStateVariables(principal_stress_vector, alpha, unknown_vector[2]);

            // Update rPrincipalStressUpdated and rPrincipalStrain as the final elastic principal stress and strain
            rPrincipalStressUpdated = principal_stress_vector;
            rPrincipalStrain = principal_strain_vector;

            region = 1;
            converged = true;
        }
    }

    return converged;
}

void BorjaCamClayPlasticFlowRule::CalculateLHSMatrix(Matrix& rLHSMatrix, const BoundedVector<double,3>& rPrincipalStressVector, const BoundedVector<double,3>& rUnknownVector, const double& rK_p)
{
    // Reset Zero
    rLHSMatrix = ZeroMatrix(3,3);

    // Compute ElasticMatrix D^e
    BoundedMatrix<double,2,2> elastic_matrix_D_e = ZeroMatrix(2,2);
    this->ComputeElasticMatrix_2X2(rPrincipalStressVector, rUnknownVector[0], rUnknownVector[1], elastic_matrix_D_e);

    // Compute Hessian Matrix H
    BoundedMatrix<double,2,2> hessian_matrix_H = ZeroMatrix(2,2);
    this->CalculateHessianMatrix_2x2(hessian_matrix_H);

    // Compute matrix_G
    BoundedMatrix<double,2,2> matrix_G = prod(hessian_matrix_H, elastic_matrix_D_e);

    // Arrange rLHSMatrix
    rLHSMatrix(0,0) = -1.0 * (1.0 + rUnknownVector[2] * (matrix_G(0,0) + rK_p * mStateFunctionSecondDerivative[5]));
    rLHSMatrix(0,1) = -1.0 * (rUnknownVector[2] * matrix_G(0,1));
    rLHSMatrix(0,2) = -1.0 * (mStateFunctionFirstDerivative[0]);

    rLHSMatrix(1,0) = -1.0 * (rUnknownVector[2] * (matrix_G(1,0) + rK_p * mStateFunctionSecondDerivative[4]));
    rLHSMatrix(1,1) = -1.0 * (1.0 + rUnknownVector[2] * matrix_G(1,1));
    rLHSMatrix(1,2) = -1.0 * (mStateFunctionFirstDerivative[1]);

    rLHSMatrix(2,0) = -1.0 * (elastic_matrix_D_e(0,0) * mStateFunctionFirstDerivative[0] + elastic_matrix_D_e(1,0) * mStateFunctionFirstDerivative[1] + rK_p * mStateFunctionFirstDerivative[2]);
    rLHSMatrix(2,1) = -1.0 * (elastic_matrix_D_e(0,1) * mStateFunctionFirstDerivative[0] + elastic_matrix_D_e(1,1) * mStateFunctionFirstDerivative[1]);
    rLHSMatrix(2,2) = 0.0;

}

void BorjaCamClayPlasticFlowRule::CalculateHessianMatrix_2x2(BoundedMatrix<double,2,2>& rHessianMatrix)
{
    // Material parameters
    const double shear_M = mpYieldCriterion->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];

    // Assemble matrix
    rHessianMatrix(0,0) = 2.0;
    rHessianMatrix(1,1) = 2.0 / std::pow(shear_M, 2);
    rHessianMatrix(0,1) = 0.0;
    rHessianMatrix(1,0) = 0.0;

}

void BorjaCamClayPlasticFlowRule::ComputeElasticMatrix_2X2(const BoundedVector<double,3>& rPrincipalStressVector, const double& rVolumetricStrain, const double& rDeviatoricStrain, BoundedMatrix<double,2,2>& rElasticMatrix)
{
    // Material parameters
    const double swelling_slope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    const double alpha_shear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];

    double ref_pressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    const double ocr = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ref_pressure /= ocr;

    const double constant_shear_modulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[INITIAL_SHEAR_MODULUS];
    const double shear_modulus = alpha_shear * ref_pressure * std::exp( -(rVolumetricStrain - mInitialVolumetricStrain) / swelling_slope);

    // Decompose principal_stress
    double mean_stress_p, deviatoric_q;
    MPMStressPrincipalInvariantsUtility::CalculateStressInvariants( rPrincipalStressVector, mean_stress_p, deviatoric_q);

    // Assemble matrix
    rElasticMatrix(0,0) = - mean_stress_p / swelling_slope;
    rElasticMatrix(1,1) = 3.0 * (constant_shear_modulus - shear_modulus);
    rElasticMatrix(0,1) = 3.0 * shear_modulus * rDeviatoricStrain / swelling_slope;
    rElasticMatrix(1,0) = rElasticMatrix(0,1);
}

void BorjaCamClayPlasticFlowRule::ComputePlasticMatrix_2X2(const BoundedVector<double,3>& rPrincipalStressVector, const double& rVolumetricStrain, const double& rDeviatoricStrain, const BoundedMatrix<double,2,2>& rElasticMatrix, BoundedMatrix<double,2,2>& rPlasticMatrix)
{
    // Initialize used temporary matrices and vectors
    BoundedVector<double,2> a = ZeroVector(2);
    BoundedMatrix<double,2,2> b = ZeroMatrix(2,2);
    BoundedVector<double,2> c = ZeroVector(2);
    BoundedVector<double,2> d = ZeroVector(2);
    double e;

    // Compute Hessian Matrix H
    BoundedMatrix<double,2,2> hessian_matrix_H = ZeroMatrix(2,2);
    this->CalculateHessianMatrix_2x2(hessian_matrix_H);

    // Compute matrix_G
    BoundedMatrix<double,2,2> matrix_G = prod(hessian_matrix_H, rElasticMatrix);

    // Initiate coefficient parameters
    const double k_p      = mMaterialParameters.PlasticHardeningModulus;
    const double k_p_trial = -k_p;
    const double delta_phi = mMaterialParameters.ConsistencyParameter;

    // Construct Matrix b
    b(0,0) = 1.0 + delta_phi * ( matrix_G(0,0) + k_p * mStateFunctionSecondDerivative[5] );
    b(0,1) = delta_phi * matrix_G(0,1);
    b(1,0) = delta_phi * ( matrix_G(1,0) + k_p * mStateFunctionSecondDerivative[4] );
    b(1,1) = 1.0 + delta_phi * matrix_G(1,1);
    double det_b = MathUtils<double>::Det(b);

    // Construct Vector c
    c[0] = 1.0 - delta_phi * k_p_trial * mStateFunctionSecondDerivative[5];
    c[1] = - delta_phi * k_p_trial * mStateFunctionSecondDerivative[4];

    // Construct Vector d
    d[0] = rElasticMatrix(0,0) * mStateFunctionFirstDerivative[0] + rElasticMatrix(1,0) * mStateFunctionFirstDerivative[1] + k_p * mStateFunctionFirstDerivative[2];
    d[1] = rElasticMatrix(0,1) * mStateFunctionFirstDerivative[0] + rElasticMatrix(1,1) * mStateFunctionFirstDerivative[1];

    // Construct Coefficient e
    e  = d[0] * ( b(1,1) * mStateFunctionFirstDerivative[0] - b(0,1) * mStateFunctionFirstDerivative[1] )
       + d[1] * ( b(0,0) * mStateFunctionFirstDerivative[1] - b(1,0) * mStateFunctionFirstDerivative[0] );

    // Construct Vector a
    a[0] = ( d[0] * (b(1,1) * c[0] - b(0,1) * c[1] ) + d[1] * (b(0,0) * c[1] - b(1,0) * c[0]) + det_b * k_p_trial * mStateFunctionFirstDerivative[2] );
    a[1] = std::sqrt(2.0/3.0) * ( d[1] * b(0,0) - d[0] * b(0,1) );

    // Check if e == 0
    if (std::abs(e) < 1.e-9){ a *= 1.0 / 1.e-9; }
    else{ a *= 1.0 / e; }

    // Arrange rPlasticMatrix from all the constructed variables
    rPlasticMatrix(0,0) = b(1,1) * (c[0] - a[0] * mStateFunctionFirstDerivative[0]) - b(0,1) * (c[1] - a[0] * mStateFunctionFirstDerivative[1]);
    rPlasticMatrix(0,1) = b(0,1) * (-1.0 + std::sqrt(3.0/2.0) * a[1] * mStateFunctionFirstDerivative[1]) - std::sqrt(3.0/2.0) * b(1,1) * a[1] * mStateFunctionFirstDerivative[0];
    rPlasticMatrix(1,0) = b(0,0) * (c[1] - a[0] * mStateFunctionFirstDerivative[1]) - b(1,0) * (c[0] - a[0] * mStateFunctionFirstDerivative[0]);
    rPlasticMatrix(1,1) = b(0,0) * (1.0  - std::sqrt(3.0/2.0) * a[1] * mStateFunctionFirstDerivative[1]) + std::sqrt(3.0/2.0) * b(1,0) * a[1] * mStateFunctionFirstDerivative[0];

    // Check if det_b == 0
    if (std::abs(det_b) < 1.e-9){ rPlasticMatrix *= 1.0 / 1.e-9; }
    else{ rPlasticMatrix *= 1.0 / det_b; }

}

// Compute Trial elastic principal stress matrix from Trial elastic principal strain matrix
void BorjaCamClayPlasticFlowRule::CalculatePrincipalStressTrial(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix)
{
    BoundedVector<double,3> main_strain      = ZeroVector(3);
    BoundedVector<double,3> principal_stress = ZeroVector(3);

    for (unsigned int i = 0; i<3; ++i)
    {
        main_strain[i] = rNewElasticLeftCauchyGreen(i,i);
    }

    this->CalculatePrincipalStressVector(main_strain, principal_stress);

    // Evalute the Kirchhoff principal stress
    for(unsigned int i=0; i<3; i++)
    {
        rStressMatrix(i,i) = principal_stress[i];
    }

}


// Function to compute Principal Stress Vector from Principal Strain Vector
void BorjaCamClayPlasticFlowRule::CalculatePrincipalStressVector(const BoundedVector<double,3>& rPrincipalStrain, BoundedVector<double,3>& rPrincipalStress)
{
    // Calculate volumetric and deviatoric strains from princial strain
    double volumetric_strain, deviatoric_strain;
    BoundedVector<double,3> deviatoric_strain_vector;
    this->CalculateStrainInvariantsFromPrincipalStrain(rPrincipalStrain, volumetric_strain, deviatoric_strain, deviatoric_strain_vector);

    // Calculate mean_stress_p and DeviatoricStressQ
    double mean_stress_p;
    this->CalculateMeanStress(volumetric_strain, deviatoric_strain, mean_stress_p);
    this->CalculateDeviatoricStress(volumetric_strain, deviatoric_strain_vector, rPrincipalStress);

    // Combine into principal_stress
    for (unsigned int i = 0; i<3; ++i)
        rPrincipalStress[i] += mean_stress_p;

}


// Function calculate Mean Stress from given volumetric_strain and deviatoric_strain
void BorjaCamClayPlasticFlowRule::CalculateMeanStress(const double& rVolumetricStrain, const double& rDeviatoricStrain, double& rMeanStress)
{
    const double swelling_slope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    const double alpha_shear    = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];

    double ref_pressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    const double ocr = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ref_pressure /= ocr;

    rMeanStress = ref_pressure * std::exp( -(rVolumetricStrain - mInitialVolumetricStrain) / swelling_slope) * (1.0 + 1.5 * alpha_shear * std::pow(rDeviatoricStrain, 2) / swelling_slope);

}


// Function calculate DeviatoricStressVector from given volumetric_strain and deviatoric_strain_vector
void BorjaCamClayPlasticFlowRule::CalculateDeviatoricStress(const double& rVolumetricStrain, const BoundedVector<double,3> & rDeviatoricStrainVector, BoundedVector<double,3>& rDeviatoricStress)
{
    double ref_pressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    const double ocr = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ref_pressure /= ocr;

    const double swelling_slope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    const double alpha_shear    = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];
    const double constant_shear_modulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[INITIAL_SHEAR_MODULUS];

    rDeviatoricStress = rDeviatoricStrainVector;
    const double shear_modulus = alpha_shear * ref_pressure * std::exp( -(rVolumetricStrain - mInitialVolumetricStrain) / swelling_slope);
    rDeviatoricStress *= 2.0 * ( constant_shear_modulus - shear_modulus );

}

// Function which returns principal strains from volumetric and deviatoric strain components
void BorjaCamClayPlasticFlowRule::CalculatePrincipalStrainFromStrainInvariants(BoundedVector<double,3>& rPrincipalStrain, const double& rVolumetricStrain, const double& rDeviatoricStrain, const BoundedVector<double,3>& rDirectionVector)
{
    rPrincipalStrain = ZeroVector(3);

    for (unsigned int i = 0; i<3; ++i)
    {
        rPrincipalStrain[i] += 1.0/3.0 * rVolumetricStrain;
    }
    rPrincipalStrain += std::sqrt(3.0/2.0) * rDeviatoricStrain * rDirectionVector;
}


// Function which returns volumetric and deviatoric strain components from principal strain
void BorjaCamClayPlasticFlowRule::CalculateStrainInvariantsFromPrincipalStrain(const BoundedVector<double,3>& rPrincipalStrain, double& rVolumetricStrain, double& rDeviatoricStrain, BoundedVector<double,3>& rDeviatoricStrainVector)
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
void BorjaCamClayPlasticFlowRule::ReturnStressFromPrincipalAxis(const Matrix& rEigenVectors, const BoundedVector<double,3>& rPrincipalStress, Matrix& rStressMatrix)
{
    rStressMatrix = ZeroMatrix(3,3);
    Vector aux_N  = ZeroVector(3);
    BoundedMatrix<double,3,3> aux_M  = ZeroMatrix(3,3);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
            aux_N[j] = rEigenVectors(j,i);
        aux_M = MathUtils<double>::TensorProduct3(aux_N, aux_N);
        rStressMatrix += rPrincipalStress[i]*aux_M;
    }
}

// Function that compute the consistent tangent stiffness matrix (in normal space) considering both elastic and elasto-plastic case
void BorjaCamClayPlasticFlowRule::ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& alfa, Matrix& rConsistMatrix)
{
    // Prepare principal_stress_vector and its invariants
    BoundedVector<double,3> principal_stress_vector = mPrincipalStressUpdated;
    double mean_stress_p, deviatoric_q;
    MPMStressPrincipalInvariantsUtility::CalculateStressInvariants( principal_stress_vector, mean_stress_p, deviatoric_q);

    // Compute StrainComponents and direction -- with check if deviatoric_strain == 0
    double volumetric_strain, deviatoric_strain;
    BoundedVector<double,3> deviatoric_strain_vector;
    this->CalculateStrainInvariantsFromPrincipalStrain(mElasticPrincipalStrain, volumetric_strain, deviatoric_strain, deviatoric_strain_vector);

    BoundedVector<double,3> direction_strain_vector = ZeroVector(3);
    if (std::abs(deviatoric_strain) > 1.e-9)
        direction_strain_vector = std::sqrt(2.0/3.0) * deviatoric_strain_vector / deviatoric_strain;

    // Compute ElasticMatrix (2x2) D^e
    BoundedMatrix<double,2,2> elastic_matrix_D_e = ZeroMatrix(2,2);
    this->ComputeElasticMatrix_2X2(principal_stress_vector, volumetric_strain, deviatoric_strain, elastic_matrix_D_e);

    // Compute PlasticMatrix (2x2) D^p
    BoundedMatrix<double,2,2> plastic_matrix_D_p = IdentityMatrix(2);
    if (rReturnMappingVariables.Options.Is(MPMFlowRule::PLASTIC_REGION))
    {
        this->ComputePlasticMatrix_2X2(principal_stress_vector, volumetric_strain, deviatoric_strain, elastic_matrix_D_e, plastic_matrix_D_p);
    }

    // Compute Elasto-plastic (2x2) D^ep = D^e D^p
    BoundedMatrix<double,2,2> matrix_D_ep = prod(elastic_matrix_D_e, plastic_matrix_D_p);

    // Prepare fourth_order_identity and identity_cross
    BoundedMatrix<double,6,6> fourth_order_identity = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i)
        fourth_order_identity(i,i) = 1.0;

    for (unsigned int i = 3; i<6; ++i)
        fourth_order_identity(i,i) = 0.50;

    BoundedMatrix<double,6,6> identity_cross = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            identity_cross(i,j) = 1.0;
        }
    }

    // Prepare tensor_NxN, tensor_1xN, and tensor_Nx1
    BoundedMatrix<double,6,6> tensor_NxN = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            tensor_NxN(i,j) = direction_strain_vector[i] * direction_strain_vector[j];
        }
    }

    BoundedMatrix<double,6,6> tensor_1xN = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            tensor_1xN(i,j) = direction_strain_vector[j];
        }
    }

    BoundedMatrix<double,6,6> tensor_Nx1 = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            tensor_Nx1(i,j) = direction_strain_vector[i];
        }
    }

    // Perform check in case deviatoric_strain == 0
    double deviatoric_q_by_strain;
    if (std::abs(deviatoric_strain) < 1.e-9 ) {deviatoric_q_by_strain = deviatoric_q / 1.e-9;}
    else{deviatoric_q_by_strain = deviatoric_q / deviatoric_strain;}

    // Compute Consistent Tangent Stiffness matrix in principal space
    BoundedMatrix<double,6,6> D_elasto_plastic = ZeroMatrix(6,6);
    D_elasto_plastic  = ( matrix_D_ep(0,0) - 2.0 * deviatoric_q_by_strain / 9.0 ) * identity_cross;
    D_elasto_plastic += ( std::sqrt(2.0/3.0) * matrix_D_ep(0,1) ) * tensor_1xN;
    D_elasto_plastic += ( std::sqrt(2.0/3.0) * matrix_D_ep(1,0) ) * tensor_Nx1;
    D_elasto_plastic += ( 2.0 * deviatoric_q_by_strain / 3.0 ) * (fourth_order_identity - tensor_NxN);
    D_elasto_plastic += ( 2.0 / 3.0 * matrix_D_ep(1,1) ) * tensor_NxN;

    // Return constitutive matrix from principal space to normal space
    BoundedMatrix<double,6,6> A = ZeroMatrix(6,6);
    BoundedMatrix<double,6,6> A_trans = ZeroMatrix(6,6);
    this->CalculateTransformationMatrix(rReturnMappingVariables.MainDirections, A);
    A_trans = trans(A);

    BoundedMatrix<double,6,6> aux_mat = ZeroMatrix(6,6);
    aux_mat = prod(A_trans, D_elasto_plastic);
    rConsistMatrix = prod(aux_mat, A);

}

void BorjaCamClayPlasticFlowRule::CalculateTransformationMatrix(const BoundedMatrix<double,3,3>& rMainDirection, BoundedMatrix<double,6,6>& rA)
{
    BoundedMatrix<double,3,3> A1 = ZeroMatrix(3);
    BoundedMatrix<double,3,3> A2 = ZeroMatrix(3);
    BoundedMatrix<double,3,3> A3 = ZeroMatrix(3);
    BoundedMatrix<double,3,3> A4 = ZeroMatrix(3);
    for (unsigned int i = 0; i<3 ; i++)
    {
        for(unsigned int j = 0; j<3 ; j++)
        {
            A1(i,j) = rMainDirection(i,j) * rMainDirection(i,j);
            rA(i,j) = A1(i,j);
        }
    }
    BoundedVector<double,3> Hj1 = ZeroVector(3);
    Hj1[0] = 0;
    Hj1[1] = 2;
    Hj1[2] = 1;

    BoundedVector<double,3> Hj2 = ZeroVector(3);
    Hj2[0] = 1;
    Hj2[1] = 0;
    Hj2[2] = 2;

    for(unsigned int k = 0; k<3; k++)
    {
        for(unsigned int l = 0; l<3; l++)
        {
            A2(k,l) = rMainDirection(k,Hj1[l]) * rMainDirection(k,Hj2[l]);
            A3(k,l) = rMainDirection(Hj1[k],l) * rMainDirection(Hj2[k],l);
            A4(k,l) = rMainDirection(Hj1[k],Hj1[l]) * rMainDirection(Hj2[k],Hj2[l]) + rMainDirection(Hj2[k],Hj1[l]) * rMainDirection(Hj1[k],Hj2[l]);
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
    BoundedVector<double,3> landa_2 = ZeroVector(3);

    for (unsigned int i = 0; i<3; ++i)
        landa_2[i] = std::exp(2.0*mElasticPrincipalStrain[i]);

    Matrix output = ZeroMatrix(3,3);
    this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, landa_2, output);

    return output;
}

// Function that updates internal variable at every time step once the nonlinear iteration converges
bool BorjaCamClayPlasticFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{
    // Compute Delta Plastic Strain
    double norm_plastic_principal_strain = norm_2(mPlasticPrincipalStrain);

    // Compute Strain Components and its invariants
    double volumetric_strain, deviatoric_strain;
    BoundedVector<double,3> deviatoric_strain_vector;
    this->CalculateStrainInvariantsFromPrincipalStrain(mPlasticPrincipalStrain, volumetric_strain, deviatoric_strain, deviatoric_strain_vector);

    // Update Equivalent Plastic Strain
    mInternalVariables.DeltaPlasticStrain = norm_plastic_principal_strain;
    mInternalVariables.EquivalentPlasticStrain += norm_plastic_principal_strain;

    // Update Accumulated Plastic Volumetric Strain
    mInternalVariables.DeltaPlasticVolumetricStrain = volumetric_strain;
    mInternalVariables.AccumulatedPlasticVolumetricStrain += volumetric_strain;

    // Update Accumulated Plastic Deviatoric Strain
    mInternalVariables.DeltaPlasticDeviatoricStrain = deviatoric_strain;
    mInternalVariables.AccumulatedPlasticDeviatoricStrain += deviatoric_strain;

    // Update Preconsolidation Stress for the next time step
    double new_preconsolidation_stress = mpYieldCriterion->GetHardeningLaw().CalculateHardening(new_preconsolidation_stress, volumetric_strain, mMaterialParameters.PreconsolidationPressure);
    mMaterialParameters.PreconsolidationPressure = new_preconsolidation_stress;

    return true;
}

void BorjaCamClayPlasticFlowRule::UpdateStateVariables(const BoundedVector<double,3> rPrincipalStress, const double rAlpha, const double rConsistencyParameter)
{
    // Calculate final state function and derivatives
    mStateFunction = mpYieldCriterion->CalculateYieldCondition(mStateFunction, rPrincipalStress, rAlpha, mMaterialParameters.PreconsolidationPressure);
    mpYieldCriterion->CalculateYieldFunctionDerivative(rPrincipalStress, mStateFunctionFirstDerivative, rAlpha, mMaterialParameters.PreconsolidationPressure);
    mpYieldCriterion->CalculateYieldFunctionSecondDerivative(rPrincipalStress, mStateFunctionSecondDerivative);

    // Update plastic hardening modulus k_p
    const double swelling_slope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    const double other_slope    = mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];
    double k_p = mpYieldCriterion->GetHardeningLaw().CalculateHardening(k_p, rAlpha, mMaterialParameters.PreconsolidationPressure);
    k_p *= 1.0 / (other_slope-swelling_slope);
    mMaterialParameters.PlasticHardeningModulus = k_p;

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
