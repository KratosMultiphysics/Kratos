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
#include "includes/ublas_interface.h"

// Project includes
#include "custom_constitutive/flow_rules/mc_plastic_flow_rule.hpp"

#include "particle_mechanics_application.h"
#include "custom_utilities/mpm_stress_principal_invariants_utility.h"

namespace Kratos
{



//************ CONSTRUCTOR ***********
MCPlasticFlowRule::MCPlasticFlowRule()
    :MPMFlowRule()
{
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

MCPlasticFlowRule::MCPlasticFlowRule(YieldCriterionPointer pYieldCriterion)
    :MPMFlowRule(pYieldCriterion)
{

}

//********* ASSIGMENT OPERATOR
MCPlasticFlowRule& MCPlasticFlowRule::operator=(MCPlasticFlowRule const& rOther)
{
    MPMFlowRule::operator=(rOther);
    return *this;

}



//********** COPY CONSTRUCTOR *********
MCPlasticFlowRule::MCPlasticFlowRule(MCPlasticFlowRule const& rOther)
    :MPMFlowRule(rOther)
{
}

//*******   CLONE ********
MPMFlowRule::Pointer MCPlasticFlowRule::Clone() const
{
    MPMFlowRule::Pointer p_clone(new MCPlasticFlowRule(*this));
    return p_clone;
}



// ********** DESTRUCTOR **************
MCPlasticFlowRule::~MCPlasticFlowRule()
{
}

void MCPlasticFlowRule::InitializeMaterial(YieldCriterionPointer& pYieldCriterionPointer, HardeningLawPointer& pHardeningPointer, const Properties& rProp)
{
    MPMFlowRule::InitializeMaterial(pYieldCriterionPointer, pHardeningPointer, rProp);

    mElasticPrincipalStrain = ZeroVector(3);
    mPlasticPrincipalStrain = ZeroVector(3);
    mElasticPreviousPrincipalStrain = ZeroVector(3);
    mPrincipalStressTrial = ZeroVector(3);
    mPrincipalStressUpdated = ZeroVector(3);
    mLargeStrainBool = true;
    mRegion = 0;

    mEquivalentPlasticStrain = 0.0;

    this->InitializeMaterialParameters();
}

// Initialize material parameters which are able to change
void MCPlasticFlowRule::InitializeMaterialParameters(){
    mMaterialParameters.Cohesion      = mpYieldCriterion->GetHardeningLaw().GetProperties()[COHESION];
    mMaterialParameters.FrictionAngle = mpYieldCriterion->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
    mMaterialParameters.DilatancyAngle= mpYieldCriterion->GetHardeningLaw().GetProperties()[INTERNAL_DILATANCY_ANGLE];
}

bool MCPlasticFlowRule::CalculateReturnMapping(
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
        // the rStressMatrix is the precomputed principal stress
        principal_stress[i] = rStressMatrix(i,i);
    }

    // Sorting Principal Stress and Strain - "0" is the largest one and "2" is the lowest one
    MPMStressPrincipalInvariantsUtility::SortPrincipalStress(principal_stress, main_strain, rReturnMappingVariables.MainDirections);

    mPrincipalStressTrial   = principal_stress;
    mElasticPrincipalStrain = main_strain;
    mElasticPreviousPrincipalStrain = main_strain;

    // Check for the yield Condition -- calling the yield criterion
    rReturnMappingVariables.TrialStateFunction = 0.0;
    rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, principal_stress, mMaterialParameters.Cohesion, mMaterialParameters.FrictionAngle);

    // If yield is reached, do return mapping
    if (rReturnMappingVariables.TrialStateFunction <= 0.0)
    {
        mRegion = 0;
        mPrincipalStressUpdated = principal_stress;
        plasticity_active = false;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);
    }
    else
    {
        unsigned int region = 0;
        BoundedVector<double,3> principal_stress_updated = ZeroVector(3);

        bool converged = this->CalculateConsistencyCondition(rReturnMappingVariables, principal_stress, mElasticPrincipalStrain, region, principal_stress_updated);
        KRATOS_ERROR_IF(!converged) << "Warning:: Constitutive Law does not converge! "<<std::endl;

        mRegion = region;
        mPrincipalStressUpdated = principal_stress_updated;

        plasticity_active = true;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,true);
    }

    // rStressMatrix is the matrix of the updated stress in cartesian configuration -- this function perform back transformation
    this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, mPrincipalStressUpdated, rStressMatrix);

    BoundedVector<double,3> DeltaPrincipalStress = principal_stress - mPrincipalStressUpdated;

    // Updated the PrincipalStrain vector
    BoundedMatrix<double,3,3> inv_elastic_matrix = ZeroMatrix(3,3);
    this->CalculateInverseElasticMatrix(rReturnMappingVariables, inv_elastic_matrix);

    // Delta plastic strain
    BoundedVector<double,3> plastic_strain = prod(inv_elastic_matrix, DeltaPrincipalStress);

    // Now the component of mElasticPrincipalStrain are sorted in the same way as plastic_strain!
    mElasticPrincipalStrain -= plastic_strain;
    mPlasticPrincipalStrain = plastic_strain;

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

bool MCPlasticFlowRule::CalculateConsistencyCondition(RadialReturnVariables& rReturnMappingVariables, const BoundedVector<double,3>& rPrincipalStress, const BoundedVector<double,3>& rPrincipalStrain, unsigned int& rRegion, BoundedVector<double,3>& rPrincipalStressUpdated)
{

    // Calculate stress return in principal stress
    // The flow rule is written for associated and non-associated plasticity
    // if the internal friction angle is = to the dilatancy angle the plasticity is associated.
    // Refer to paper by Clausen for MC flow rule for the theoretical implementation

    // Material parameters
    const double cohesion       = mMaterialParameters.Cohesion;
    const double friction_angle  = mMaterialParameters.FrictionAngle;
    const double dilatancy_angle = mMaterialParameters.DilatancyAngle;

    // Necessary coefficients
    const double friction_coefficient  = (1 + std::sin(friction_angle))/(1 - std::sin(friction_angle));
    const double dilatancy_coefficient = (1 + std::sin(dilatancy_angle))/(1 - std::sin(dilatancy_angle));
    const double cohesion_coefficient  = 2 * cohesion * std::sqrt(friction_coefficient);

    // Stress coordinate of the criterions apex
    const double apex = cohesion_coefficient/(friction_coefficient-1);

    // Compute elastic matrix which takes account only for normal stresses
    BoundedMatrix<double,3,3> D = ZeroMatrix(3,3);
    this->ComputeElasticMatrix_3X3(rReturnMappingVariables, D);

    // Compute the direction of the plastic return stress R_p
    BoundedVector<double,3> R_p = ZeroVector(3);
    double den_p = friction_coefficient *(D(0,0)*dilatancy_coefficient - D(0,2)) - D(2,0) * dilatancy_coefficient + D(2,2);
    if (std::abs(den_p) < 1.e-9) den_p = 1.e-9;
    R_p[0] = (D(0,0)*dilatancy_coefficient - D(0,2) )/den_p;
    R_p[1] = (D(1,0)*dilatancy_coefficient - D(1,2) )/den_p;
    R_p[2] = (D(2,0)*dilatancy_coefficient - D(2,2) )/den_p;

    // Vector from predictor stress to the apex
    BoundedVector<double,3> sigma_P_apex = ZeroVector(3);
    sigma_P_apex[0] = rPrincipalStress[0] - apex;
    sigma_P_apex[1] = rPrincipalStress[1] - apex;
    sigma_P_apex[2] = rPrincipalStress[2] - apex;

    // Boundary plane between region I and II: evaluated as the cross product between R_p and R1, the direction of line 1
    BoundedVector<double,3> NI_II = ZeroVector(3);
    NI_II[0] = R_p[1] * friction_coefficient - R_p[2];
    NI_II[1] = R_p[2] - R_p[0] * friction_coefficient;
    NI_II[2] = R_p[0] - R_p[1];
    const double pI_II = NI_II[0] * sigma_P_apex[0] + NI_II[1] * sigma_P_apex[1] + NI_II[2] * sigma_P_apex[2];

    // Boundary plane between region I and III: evaluated as the cross product between R_p and R2, the direction of line 2
    BoundedVector<double,3> NI_III = ZeroVector(3);
    NI_III[0] = R_p[1] * friction_coefficient - R_p[2] * friction_coefficient;
    NI_III[1] = R_p[2] - R_p[0] * friction_coefficient;
    NI_III[2] = R_p[0] * friction_coefficient - R_p[1];
    const double pI_III = NI_III[0] * sigma_P_apex[0] + NI_III[1] * sigma_P_apex[1] + NI_III[2] * sigma_P_apex[2];

    // t-parameters for region determination
    // Secondary surface in region II a = [0 k -1], b  = [0 m -1] -- needed to calculate t1
    double den_p2 = friction_coefficient * (D(1,1) * dilatancy_coefficient - D(1,2)) - D(1,2) * dilatancy_coefficient + D(2,2);
    if (std::abs(den_p2) < 1.e-9) den_p2 = 1.e-9;
    BoundedVector<double,3> R_p2 = ZeroVector(3);
    R_p2[0] = (D(0,1) * dilatancy_coefficient - D(0,2))/den_p2;
    R_p2[1] = (D(1,1) * dilatancy_coefficient - D(1,2))/den_p2;
    R_p2[2] = (D(2,1) * dilatancy_coefficient - D(2,2))/den_p2;

    BoundedVector<double,3> N2 = ZeroVector(3);
    N2[0] = R_p[1]*R_p2[2] - R_p[2]*R_p2[1];
    N2[1] = R_p[2]*R_p2[0] - R_p[0]*R_p2[2];
    N2[2] = R_p[0]*R_p2[1] - R_p[1]*R_p2[0];

    const double num1 = N2[0] * sigma_P_apex[0] + N2[1] * sigma_P_apex[1] + N2[2] * sigma_P_apex[2];
    double den_1 = N2[0] + N2[1] + friction_coefficient * N2[2];
    if (std::abs(den_1) < 1.e-9) den_1 = 1.e-9;
    const double t1 = num1 / den_1 ;

    // Secondary surface in region III a = [k -1 0], b  = [m -1 0] -- needed to calculate t2
    double den = friction_coefficient * (D(0,0) * dilatancy_coefficient - D(0,1)) - D(1,0) * dilatancy_coefficient + D(1,1);
    if (std::abs(den) < 1.e-9) den = 1.e-9;
    BoundedVector<double,3> R_p3 = ZeroVector(3);
    R_p3[0] = (D(0,0) * dilatancy_coefficient - D(0,1))/den;
    R_p3[1] = (D(1,0) * dilatancy_coefficient - D(1,1))/den;
    R_p3[2] = (D(2,0) * dilatancy_coefficient - D(2,1))/den;

    BoundedVector<double,3> N3 = ZeroVector(3);
    N3[0] = R_p[1]*R_p3[2] - R_p[2]*R_p3[1];
    N3[1] = R_p[2]*R_p3[0] - R_p[0]*R_p3[2];
    N3[2] = R_p[0]*R_p3[1] - R_p[1]*R_p3[0];

    const double num2 = N3[0] * sigma_P_apex[0] + N3[1] * sigma_P_apex[1] + N3[2] * sigma_P_apex[2];
    double den_2 = N3[0] + friction_coefficient * N3[1] + friction_coefficient * N3[2];
    if (std::abs(den_2) < 1.e-9) den_2 = 1.e-9;
    const double t2 = num2 / den_2 ;

    // rRegion detection and return determination
    // Return mapping to the apex
    if(t1 > 0 || t2 > 0) //check: both the conditions have to be satisfied
    {
        rRegion = 4;
        rPrincipalStressUpdated[0] = apex;
        rPrincipalStressUpdated[1] = apex;
        rPrincipalStressUpdated[2] = apex;
    }

    // Return mapping to line 1
    else if(pI_II < 0)
    {
        rRegion = 2;
        rPrincipalStressUpdated[0] = t1 + apex;
        rPrincipalStressUpdated[1] = t1 + apex;
        rPrincipalStressUpdated[2] = t1 * friction_coefficient + apex;
    }

    // Return mapping to the yield surface
    else if(pI_III <= 0)
    {
        rRegion = 1;
        const double state_function     = rReturnMappingVariables.TrialStateFunction;
        rPrincipalStressUpdated[0] = rPrincipalStress[0] -  state_function * R_p[0];
        rPrincipalStressUpdated[1] = rPrincipalStress[1] -  state_function * R_p[1];
        rPrincipalStressUpdated[2] = rPrincipalStress[2] -  state_function * R_p[2];
    }

    // Return mapping to line 2
    else
    {
        rRegion = 3;
        rPrincipalStressUpdated[0] = t2 + apex;
        rPrincipalStressUpdated[1] = t2 * friction_coefficient + apex;
        rPrincipalStressUpdated[2] = t2 * friction_coefficient + apex;
    }

    return true;
}

void MCPlasticFlowRule::ComputeElasticMatrix_3X3(const RadialReturnVariables& rReturnMappingVariables, BoundedMatrix<double,3,3>& rElasticMatrix)
{

    const double young_modulus  = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    const double poisson_ratio  = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];

    const double diagonal    = young_modulus/(1.0+poisson_ratio)/(1.0-2.0*poisson_ratio) * (1.0-poisson_ratio);
    const double nondiagonal = young_modulus/(1.0+poisson_ratio)/(1.0-2.0*poisson_ratio) * ( poisson_ratio);

    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            if (i == j)
            {
                rElasticMatrix(i,i) = diagonal;
            }
            else
            {
                rElasticMatrix(i,j) = nondiagonal;
            }
        }
    }
}

void MCPlasticFlowRule::CalculateInverseElasticMatrix(const RadialReturnVariables& rReturnMappingVariables, BoundedMatrix<double,3,3>& rInverseElasticMatrix)
{
    const double young_modulus  = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    const double poisson_ratio  = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];

    const double lame_lambda  = (young_modulus*poisson_ratio)/((1+poisson_ratio)*(1-2*poisson_ratio));
    const double lame_mu      =  young_modulus/(2*(1+poisson_ratio));

    const double diagonal    = (lame_lambda + lame_mu)/(lame_mu*(3.0*lame_lambda+2.0*lame_mu));
    const double nondiagonal = (-lame_lambda)/( 2.0*lame_mu*(3.0*lame_lambda + 2.0*lame_mu));

    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            if (i == j)
            {
                rInverseElasticMatrix(i,i) = diagonal;
            }
            else
            {
                rInverseElasticMatrix(i,j) = nondiagonal;
            }
        }
    }
}

void MCPlasticFlowRule::CalculateElasticMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rElasticMatrix)
{
    const double young_modulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    const double poisson_ratio = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];

    const double diagonal    = young_modulus/(1.0+poisson_ratio)/(1.0-2.0*poisson_ratio) * (1.0-poisson_ratio);
    const double nondiagonal = young_modulus/(1.0+poisson_ratio)/(1.0-2.0*poisson_ratio) * ( poisson_ratio);
    const double corte       = young_modulus/(1.0+poisson_ratio)/2.0;

    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            if (i == j)
            {
                rElasticMatrix(i,i) = diagonal;
            }
            else
            {
                rElasticMatrix(i,j) = nondiagonal;
            }
        }
    }

    for (unsigned int j = 3; j<6; ++j)
        rElasticMatrix(j,j) = corte;

}

void MCPlasticFlowRule::CalculatePrincipalStressTrial(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix)
{
    BoundedVector<double,3> main_strain = ZeroVector(3);

    for (unsigned int i = 0; i<3; ++i)
    {
        main_strain[i] = rNewElasticLeftCauchyGreen(i,i);
    }

    // Calculate the elastic matrix
    BoundedMatrix<double,3,3> ElasticMatrix = ZeroMatrix(3,3);
    const double young_modulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    const double poisson_ratio = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];
    const double diagonal    = young_modulus/(1.0+poisson_ratio)/(1.0-2.0*poisson_ratio) * (1.0-poisson_ratio);
    const double nondiagonal = young_modulus/(1.0+poisson_ratio)/(1.0-2.0*poisson_ratio) * ( poisson_ratio);

    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            if (i == j)
            {
                ElasticMatrix(i,i) = diagonal;
            }
            else
            {
                ElasticMatrix(i,j) = nondiagonal;
            }
        }
    }

    BoundedVector<double,3> principal_stress = ZeroVector(3);

    // Evalute the Kirchhoff principal stress
    principal_stress = prod(ElasticMatrix, main_strain);

    for(unsigned int i=0; i<3; i++)
    {
        rStressMatrix(i,i) = principal_stress[i];
    }
}


void MCPlasticFlowRule::ReturnStressFromPrincipalAxis(const BoundedMatrix<double,3,3>& rEigenVectors, const BoundedVector<double,3>& rPrincipalStress, Matrix& rStressMatrix)
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

// this matrix is evaluated to make consistent the elastoplastic tangent matrix
// the principal stresses are in descending order
void MCPlasticFlowRule::CalculateModificationMatrix(const RadialReturnVariables& rReturnMappingVariables, BoundedMatrix<double,3,3>& rAuxT, BoundedMatrix<double,3,3>& rInvAuxT)
{
    if(mPrincipalStressUpdated[0] - mPrincipalStressUpdated[1] > 0)
    {
        rAuxT(0,0) = (mPrincipalStressUpdated[0] - mPrincipalStressUpdated[1])/(mPrincipalStressTrial[0] - mPrincipalStressTrial[1]);
        rInvAuxT(0,0) = 1/rAuxT(0,0);
    }
    if(mPrincipalStressUpdated[0] - mPrincipalStressUpdated[2] > 0)
    {
        rAuxT(1,1) = (mPrincipalStressUpdated[0] - mPrincipalStressUpdated[2])/(mPrincipalStressTrial[0] - mPrincipalStressTrial[2]);
        rInvAuxT(1,1) = 1/rAuxT(1,1);
    }
    if(mPrincipalStressUpdated[1] - mPrincipalStressUpdated[2] > 0)
    {
        rAuxT(2,2) = (mPrincipalStressUpdated[1] - mPrincipalStressUpdated[2])/(mPrincipalStressTrial[1] - mPrincipalStressTrial[2]);
        rInvAuxT(2,2) = 1/rAuxT(2,2);
    }
}

void MCPlasticFlowRule::CalculateDepSurface(BoundedMatrix<double,3,3>& rElasticMatrix, BoundedVector<double,3>& rFNorm, BoundedVector<double,3>& rGNorm, BoundedMatrix<double,3,3>& rAuxDep)
{
    BoundedVector<double,3> aux_F = prod(trans(rFNorm), rElasticMatrix);

    BoundedVector<double,3> a = prod(rElasticMatrix, rGNorm);
    BoundedVector<double,3> b = prod(trans(rFNorm),rElasticMatrix);

    BoundedMatrix<double,3,3> num = ZeroMatrix(3,3);

    for (unsigned int i = 0; i<3; i++)
    {
        for(unsigned int j = 0; j<3; j++)
        {
            num(i,j) = a[i] * b[j];
        }
    }

    double den = MathUtils<double>::Dot(aux_F,rGNorm);
    rAuxDep = rElasticMatrix - num / den;

}

void MCPlasticFlowRule::CalculateDepLine(BoundedMatrix<double,3,3>& rInvD, BoundedVector<double,3>& rFNorm, BoundedVector<double,3>& rGNorm, BoundedMatrix<double,3,3>& rAuxDep)
{
    BoundedMatrix<double,3,3> num = ZeroMatrix(3,3);

    for (unsigned int i = 0; i<3; i++)
    {
        for(unsigned int j = 0; j<3; j++)
        {
            num(i,j) = rFNorm[i] * rGNorm[j];
        }
    }

    BoundedVector<double,3> den_1 = prod(rInvD, rGNorm);
    double den = MathUtils<double>::Dot(trans(rFNorm),den_1);

    rAuxDep = num / den;

}

void MCPlasticFlowRule::CalculateTransformationMatrix(const BoundedMatrix<double,3,3>& rMainDirection, BoundedMatrix<double,6,6>& rA)
{
    BoundedMatrix<double,3,3> A1 = ZeroMatrix(3,3);
    BoundedMatrix<double,3,3> A2 = ZeroMatrix(3,3);
    BoundedMatrix<double,3,3> A3 = ZeroMatrix(3,3);
    BoundedMatrix<double,3,3> A4 = ZeroMatrix(3,3);
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

void MCPlasticFlowRule::ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& alfa, Matrix& rConsistMatrix)
{
    // Elastoplastic constitutive matrix
    if (rReturnMappingVariables.Options.Is(MPMFlowRule::PLASTIC_REGION))
    {
        BoundedMatrix<double,6,6> t = identity_matrix<double> (6);
        BoundedMatrix<double,3,3> aux_T = ZeroMatrix(3,3);
        BoundedMatrix<double,3,3> aux_T_inv = ZeroMatrix(3,3);

        //1. Calculate  the modification matrix t
        this->CalculateModificationMatrix( rReturnMappingVariables, aux_T, aux_T_inv);
        for(unsigned int i = 3; i<6; i++)
        {
            int index_i = i-3;

            t(i,i) = t(i,i) * aux_T(index_i, index_i);
        }

        BoundedVector<double,3> d_principal_stress = mPrincipalStressTrial - mPrincipalStressUpdated;
        BoundedMatrix<double,6,6> D_ep = ZeroMatrix(6,6);

        //2. Calculate the ElastoPlastic Matrix depending on the region of return mapping
        this->CalculateElastoPlasticMatrix(rReturnMappingVariables, mRegion, d_principal_stress, D_ep);

        //3. Consistent Constitutive matrix
        BoundedMatrix<double,6,6> D_elasto_plastic = ZeroMatrix(6,6);
        D_elasto_plastic = prod(t, D_ep);

        //4. Return constitutive matrix from principal axis
        BoundedMatrix<double,6,6> A = ZeroMatrix(6,6);
        BoundedMatrix<double,6,6> A_trans = ZeroMatrix(6,6);

        //4.a Calculate transformation matrix
        this->CalculateTransformationMatrix(rReturnMappingVariables.MainDirections, A);
        A_trans = trans(A);

        BoundedMatrix<double,6,6> aux_mat = ZeroMatrix(6,6);
        aux_mat = prod(A_trans, D_elasto_plastic);
        rConsistMatrix = prod(aux_mat, A);
    }
    //Elastic matrix
    else
    {
        this->CalculateElasticMatrix(rReturnMappingVariables, rConsistMatrix);
    }

}

void MCPlasticFlowRule::CalculateElastoPlasticMatrix(const RadialReturnVariables& rReturnMappingVariables, unsigned int& rRegion, BoundedVector<double,3>& DiffPrincipalStress, BoundedMatrix<double,6,6>& rDep)
{

    const double young_modulus      = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    const double poisson_ratio      = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];
    const double shear_contribution = young_modulus/(1.0+poisson_ratio)/2.0;

    const double friction_angle  = mMaterialParameters.FrictionAngle;
    const double dilatancy_angle = mMaterialParameters.DilatancyAngle;

    const double friction_coefficient = (1 + std::sin(friction_angle))/(1 - std::sin(friction_angle));
    const double dilatancy_coefficient = (1 + std::sin(dilatancy_angle))/(1 - std::sin(dilatancy_angle));

    switch(rRegion)
    {
        // Return mapping on yield surface
        case 1:
        {
            // Yield plane normal
            BoundedVector<double,3> F_norm = ZeroVector(3);
            F_norm[0] = friction_coefficient;
            F_norm[1] = 0;
            F_norm[2] = -1;

            // Potential plane normal
            BoundedVector<double,3> G_norm = ZeroVector(3);
            G_norm[0] = dilatancy_coefficient;
            G_norm[1] = 0;
            G_norm[2] = -1;
            BoundedMatrix<double,3,3> aux_D_ep = ZeroMatrix(3,3);

            // Compute elastic matrix which takes account only for normal stresses
            BoundedMatrix<double,3,3> D = ZeroMatrix(3,3);
            this->ComputeElasticMatrix_3X3(rReturnMappingVariables, D);

            this->CalculateDepSurface(D, F_norm, G_norm, aux_D_ep);

            // Shear components of consistent constitutive matrix in principal stress space
            for (unsigned int j = 3; j<6; ++j)
                rDep(j,j) = shear_contribution;

            for (unsigned int i = 0; i<3 ; ++i)
            {
                for (unsigned int k = 0; k<3 ; ++k)
                {
                    rDep(i,k) = aux_D_ep(i,k);
                }
            }

            break;
        }

        // Return to a line 1, triaxial compression, sigp1 = sigp2
        case 2:
        {
            // Edge line direction
            BoundedVector<double,3> L_F_dir = ZeroVector(3);
            L_F_dir[0] = 1;
            L_F_dir[1] = 1;
            L_F_dir[2] = friction_coefficient;

            // Potential edge line direction
            BoundedVector<double,3> L_G_dir = ZeroVector(3);
            L_G_dir[0] = 1;
            L_G_dir[1] = 1;
            L_G_dir[2] = dilatancy_coefficient;

            BoundedMatrix<double,3,3> inv_elastic_matrix = ZeroMatrix(3,3);
            this->CalculateInverseElasticMatrix(rReturnMappingVariables, inv_elastic_matrix);
            BoundedMatrix<double,3,3> aux_D_ep = ZeroMatrix(3,3);

            this->CalculateDepLine(inv_elastic_matrix, L_F_dir, L_G_dir, aux_D_ep);

            for (unsigned int j = 3; j<6; ++j)
                rDep(j,j) = shear_contribution;

            for (unsigned int i = 0; i<3 ; ++i)
            {
                for (unsigned int k = 0; k<3 ; ++k)
                {
                    rDep(i,k) = aux_D_ep(i,k);
                }
            }

            break;
        }

        // Return to a line 2, triaxial extension, sigp2 = sigp3
        case 3:
        {
            //Edge line direction
            BoundedVector<double,3> L_F_dir = ZeroVector(3);
            L_F_dir[0] = 1;
            L_F_dir[1] = friction_coefficient;
            L_F_dir[2] = friction_coefficient;

            //Potential edge line direction
            BoundedVector<double,3> L_G_dir = ZeroVector(3);
            L_G_dir[0] = 1;
            L_G_dir[1] = dilatancy_coefficient;
            L_G_dir[2] = dilatancy_coefficient;

            BoundedMatrix<double,3,3> inv_elastic_matrix = ZeroMatrix(3,3);
            this->CalculateInverseElasticMatrix(rReturnMappingVariables, inv_elastic_matrix);
            BoundedMatrix<double,3,3> aux_D_ep = ZeroMatrix(3,3);

            this->CalculateDepLine(inv_elastic_matrix, L_F_dir, L_G_dir, aux_D_ep);

            for (unsigned int j = 3; j<6; ++j)
                rDep(j,j) = shear_contribution;

            for (unsigned int i = 0; i<3 ; ++i)
            {
                for (unsigned int k = 0; k<3 ; ++k)
                {
                    rDep(i,k) = aux_D_ep(i,k);
                }
            }

            break;
        }
    }
}


Matrix MCPlasticFlowRule::GetElasticLeftCauchyGreen(RadialReturnVariables& rReturnMappingVariables)
{

    BoundedVector<double,3> lambda_2 = ZeroVector(3);

    for (unsigned int i = 0; i<3; ++i)
        lambda_2[i] = std::exp(2.0*mElasticPrincipalStrain[i]);

    Matrix output = ZeroMatrix(3,3);
    this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, lambda_2, output);

    return output;
}

bool MCPlasticFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{
    // Compute Delta Plastic Strain
    double norm_plastic_principal_strain = norm_2(mPlasticPrincipalStrain);
    mInternalVariables.DeltaPlasticStrain = norm_plastic_principal_strain;

    // Compute Strain Components and its invariants
    double volumetric_plastic_principal_strain = sum(mPlasticPrincipalStrain);
    BoundedVector<double,3> deviatoric_plastic_principal_strain = mPlasticPrincipalStrain;
    for (unsigned int i = 0; i<3; ++i)
        deviatoric_plastic_principal_strain[i] -= 1.0/3.0 * volumetric_plastic_principal_strain;
    double plastic_deviatoric_strain = std::sqrt(2.0/3.0) * norm_2(deviatoric_plastic_principal_strain);

    const double friction_angle  = mMaterialParameters.FrictionAngle;
    const double dilatancy_angle = mMaterialParameters.DilatancyAngle;

    const double friction_coefficient = (1 + std::sin(friction_angle))/(1 - std::sin(friction_angle));
    const double dilatancy_coefficient = (1 + std::sin(dilatancy_angle))/(1 - std::sin(dilatancy_angle));

    double norm_state_function_derivative = 0.0;

    // Calculate the norm of state function (or potential) derivative
    if(friction_angle == dilatancy_angle) // I am using an associative flow rule
    {
        norm_state_function_derivative = std::sqrt(1 + friction_coefficient * friction_coefficient);
    }
    else //I am using a non-associative flow rule
    {
        norm_state_function_derivative = std::sqrt(1 + dilatancy_coefficient * dilatancy_coefficient);
    }

    // Update Equivalent Plastic Strain
    double delta_equivalent_plastic_strain = mInternalVariables.DeltaPlasticStrain / norm_state_function_derivative;
    mInternalVariables.EquivalentPlasticStrain    += delta_equivalent_plastic_strain;

    // Update Accumulated Plastic Deviatoric Strain
    mInternalVariables.DeltaPlasticDeviatoricStrain = plastic_deviatoric_strain;
    mInternalVariables.AccumulatedPlasticDeviatoricStrain += plastic_deviatoric_strain;

    return true;
}


double MCPlasticFlowRule::GetPI()
{
    return std::atan(1.0)*4.0;
}

unsigned int MCPlasticFlowRule::GetPlasticRegion()
{
    return mRegion;
}


void MCPlasticFlowRule::ComputePlasticHardeningParameter(const BoundedVector<double,3>& rHenckyStrainVector, const double& rAlpha, double& rH)
{
    rH = 0.0;
}

void MCPlasticFlowRule::save( Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMFlowRule )
}

void MCPlasticFlowRule::load( Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMFlowRule )

}

} //end namespace kratos
