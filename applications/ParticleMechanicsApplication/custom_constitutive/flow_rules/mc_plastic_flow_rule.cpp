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
#include<cmath>

// External includes
#include "includes/ublas_interface.h"

// Project includes
#include "custom_constitutive/flow_rules/mc_plastic_flow_rule.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

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

bool MCPlasticFlowRule::CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{
    bool PlasticityActive = false;
    rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);
    
    Vector PrincipalStress = ZeroVector(3);
    Vector MainStrain      = ZeroVector(3);
    
    for (unsigned int i = 0; i<3; ++i)
        MainStrain[i] = rNewElasticLeftCauchyGreen(i,i);

    for(unsigned int i=0; i<3; i++)
    {
        // the rStressMatrix is the precomputed principal stress
        PrincipalStress(i) = rStressMatrix(i,i);
    }

    // Sorting Principal Stress and Strain - "0" is the largest one and "2" is the lowest one
    MPMStressPrincipalInvariantsUtility::SortPrincipalStress(PrincipalStress, MainStrain, rReturnMappingVariables.MainDirections);

    mPrincipalStressTrial = PrincipalStress;
    mElasticPrincipalStrain = MainStrain;
    mElasticPreviousPrincipalStrain = MainStrain;

    // Check for the yield Condition -- calling the yield criterion
    rReturnMappingVariables.TrialStateFunction = 0.0;
    rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PrincipalStress, mMaterialParameters.Cohesion, mMaterialParameters.FrictionAngle);
    
    // If yield is reached, do return mapping
    if (rReturnMappingVariables.TrialStateFunction <= 0.0)
    {
        mRegion = 0;
        mPrincipalStressUpdated = PrincipalStress;
        PlasticityActive = false;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);
    }
    else
    {
        unsigned int Region = 0;
        Vector PrincipalStressUpdated = ZeroVector(3);

        bool converged = this->CalculateConsistencyCondition(rReturnMappingVariables, PrincipalStress, mElasticPrincipalStrain, Region, PrincipalStressUpdated);
        KRATOS_ERROR_IF(!converged) << "Warning:: Constitutive Law does not converge! "<<std::endl;

        mRegion = Region;
        mPrincipalStressUpdated = PrincipalStressUpdated;

        PlasticityActive = true;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,true);
    }

    // rStressMatrix is the matrix of the updated stress in cartesian configuration -- this function perform back transformation
    this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, mPrincipalStressUpdated, rStressMatrix);

    Vector DeltaPrincipalStress = PrincipalStress - mPrincipalStressUpdated; 

    // Updated the PrincipalStrain vector
    Matrix InverseElasticMatrix = ZeroMatrix(3,3);
    this->CalculateInverseElasticMatrix(rReturnMappingVariables, InverseElasticMatrix);
 
    // Delta plastic strain
    Vector PlasticStrain = prod(InverseElasticMatrix, DeltaPrincipalStress);

    // Now the component of mElasticPrincipalStrain are sorted in the same way as PlasticStrain!!!!!!!!!
    mElasticPrincipalStrain -= PlasticStrain;
    mPlasticPrincipalStrain = PlasticStrain;

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

bool MCPlasticFlowRule::CalculateConsistencyCondition(RadialReturnVariables& rReturnMappingVariables, Vector& rPrincipalStress, Vector& rPrincipalStrain, unsigned int& region, Vector& rPrincipalStressUpdated)
{

    // Calculate stress return in principal stress
    // The flow rule is written for associated and non-associated plasticity
    // if the internal friction angle is = to the dilatancy angle the plasticity is associated.
    // Refer to paper by Clausen for MC flow rule for the theoretical implementation

    // Material parameters
    const double Cohesion       = mMaterialParameters.Cohesion;
    const double FrictionAngle  = mMaterialParameters.FrictionAngle;
    const double DilatancyAngle = mMaterialParameters.DilatancyAngle;

    // Necessary coefficients
    const double FrictionCoefficient  = (1 + std::sin(FrictionAngle))/(1 - std::sin(FrictionAngle));
    const double DilatancyCoefficient = (1 + std::sin(DilatancyAngle))/(1 - std::sin(DilatancyAngle));
    const double CohesionCoefficient  = 2 * Cohesion * sqrt(FrictionCoefficient);

    // Stress coordinate of the criterions apex
    const double apex = CohesionCoefficient/(FrictionCoefficient-1);

    // Compute elastic matrix which takes account only for normal stresses
    Matrix D = ZeroMatrix(3,3);
    this->ComputeElasticMatrix_3X3(rReturnMappingVariables, D);

    // Compute the direction of the plastic return stress Rp
    Vector Rp = ZeroVector(3);
    double denp = FrictionCoefficient *(D(0,0)*DilatancyCoefficient - D(0,2)) - D(2,0) * DilatancyCoefficient + D(2,2);
    Rp(0) = (D(0,0)*DilatancyCoefficient - D(0,2) )/denp;
    Rp(1) = (D(1,0)*DilatancyCoefficient - D(1,2) )/denp;
    Rp(2) = (D(2,0)*DilatancyCoefficient - D(2,2) )/denp;

    // Vector from predictor stress to the apex
    Vector SigmaPApex = ZeroVector(3);
    SigmaPApex(0) = rPrincipalStress(0) - apex;
    SigmaPApex(1) = rPrincipalStress(1) - apex;
    SigmaPApex(2) = rPrincipalStress(2) - apex;

    // Boundary plane between region I and II: evaluated as the cross product between Rp and R1, the direction of line 1
    Vector NI_II = ZeroVector(3);
    NI_II(0) = Rp(1) * FrictionCoefficient - Rp(2);
    NI_II(1) = Rp(2) - Rp(0) * FrictionCoefficient;
    NI_II(2) = Rp(0) - Rp(1);
    const double pI_II = NI_II(0) * SigmaPApex(0) + NI_II(1) * SigmaPApex(1) + NI_II(2) * SigmaPApex(2);

    // Boundary plane between region I and III: evaluated as the cross product between Rp and R2, the direction of line 2
    Vector NI_III = ZeroVector(3);
    NI_III(0) = Rp(1) * FrictionCoefficient - Rp(2) * FrictionCoefficient;
    NI_III(1) = Rp(2) - Rp(0) * FrictionCoefficient;
    NI_III(2) = Rp(0) * FrictionCoefficient - Rp(1);
    const double pI_III = NI_III(0) * SigmaPApex(0) + NI_III(1) * SigmaPApex(1) + NI_III(2) * SigmaPApex(2);

    // t-parameters for region determination
    // Secondary surface in region II a = [0 k -1], b  = [0 m -1] -- needed to calculate t1
    double denp2 = FrictionCoefficient * (D(1,1) * DilatancyCoefficient - D(1,2)) - D(1,2) * DilatancyCoefficient + D(2,2);
    Vector Rp2 = ZeroVector(3);
    Rp2(0) = (D(0,1) * DilatancyCoefficient - D(0,2))/denp2;
    Rp2(1) = (D(1,1) * DilatancyCoefficient - D(1,2))/denp2;
    Rp2(2) = (D(2,1) * DilatancyCoefficient - D(2,2))/denp2;

    Vector N2 = ZeroVector(3);
    N2(0) = Rp(1)*Rp2(2) - Rp(2)*Rp2(1);
    N2(1) = Rp(2)*Rp2(0) - Rp(0)*Rp2(2);
    N2(2) = Rp(0)*Rp2(1) - Rp(1)*Rp2(0);

    double num1 = N2(0) * SigmaPApex(0) + N2(1) * SigmaPApex(1) + N2(2) * SigmaPApex(2);
    double den1 = N2(0) + N2(1) + FrictionCoefficient * N2(2);
    const double t1 = num1 / den1 ;

    // Secondary surface in region III a = [k -1 0], b  = [m -1 0] -- needed to calculate t2
    double den = FrictionCoefficient * (D(0,0) * DilatancyCoefficient - D(0,1)) - D(1,0) * DilatancyCoefficient + D(1,1);
    Vector Rp3 = ZeroVector(3);
    Rp3(0) = (D(0,0) * DilatancyCoefficient - D(0,1))/den;
    Rp3(1) = (D(1,0) * DilatancyCoefficient - D(1,1))/den;
    Rp3(2) = (D(2,0) * DilatancyCoefficient - D(2,1))/den;

    Vector N3 = ZeroVector(3);
    N3(0) = Rp(1)*Rp3(2) - Rp(2)*Rp3(1);
    N3(1) = Rp(2)*Rp3(0) - Rp(0)*Rp3(2);
    N3(2) = Rp(0)*Rp3(1) - Rp(1)*Rp3(0);

    double num2 = N3(0) * SigmaPApex(0) + N3(1) * SigmaPApex(1) + N3(2) * SigmaPApex(2);
    double den2 = N3(0) + FrictionCoefficient * N3(1) + FrictionCoefficient * N3(2);
    const double t2 = num2 / den2 ;

    // Region detection and return determination
    // Return mapping to the apex
    if(t1 > 0 || t2 > 0) //check: both the conditions have to be satisfied
    {
        region = 4;
        rPrincipalStressUpdated(0) = apex;
        rPrincipalStressUpdated(1) = apex;
        rPrincipalStressUpdated(2) = apex;
    }
    
    // Return mapping to line 1     
    else if(pI_II < 0)
    {
        region = 2;
        rPrincipalStressUpdated(0) = t1 + apex;
        rPrincipalStressUpdated(1) = t1 + apex;
        rPrincipalStressUpdated(2) = t1 * FrictionCoefficient + apex;
    }

    // Return mapping to the yield surface
    else if(pI_III <= 0)
    {
        region = 1;
        double StateFunction     = rReturnMappingVariables.TrialStateFunction;
        rPrincipalStressUpdated(0) = rPrincipalStress(0) -  StateFunction * Rp(0);
        rPrincipalStressUpdated(1) = rPrincipalStress(1) -  StateFunction * Rp(1);
        rPrincipalStressUpdated(2) = rPrincipalStress(2) -  StateFunction * Rp(2);
    }

    // Return mapping to line 2
    else
    {
        region = 3;
        rPrincipalStressUpdated(0) = t2 + apex;
        rPrincipalStressUpdated(1) = t2 * FrictionCoefficient + apex;
        rPrincipalStressUpdated(2) = t2 * FrictionCoefficient + apex;
    }

    return true;
}

void MCPlasticFlowRule::ComputeElasticMatrix_3X3(const RadialReturnVariables& rReturnMappingVariables, Matrix& rElasticMatrix)
{

    const double YoungModulus        = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    const double PoissonCoefficient  = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];

    const double diagonal   = YoungModulus/(1.0+PoissonCoefficient)/(1.0-2.0*PoissonCoefficient) * (1.0-PoissonCoefficient);
    const double nodiagonal = YoungModulus/(1.0+PoissonCoefficient)/(1.0-2.0*PoissonCoefficient) * ( PoissonCoefficient);

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
                rElasticMatrix(i,j) = nodiagonal;
            }
        }
    }
}

void MCPlasticFlowRule::CalculateInverseElasticMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rInverseElasticMatrix)
{
    const double YoungModulus        = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    const double PoissonCoefficient  = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];

    const double LameLambda         = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    const double LameMu             =  YoungModulus/(2*(1+PoissonCoefficient));

    const double Diagonal    = (LameLambda + LameMu)/(LameMu*(3.0*LameLambda+2.0*LameMu));
    const double NonDiagonal = (-LameLambda)/( 2.0*LameMu*(3.0*LameLambda + 2.0*LameMu));

    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            if (i == j)
            {
                rInverseElasticMatrix(i,i) = Diagonal;
            }
            else
            {
                rInverseElasticMatrix(i,j) = NonDiagonal;
            }
        }
    }
}

void MCPlasticFlowRule::CalculateElasticMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rElasticMatrix)
{
    const double Young      = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    const double Nu         = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];
    
    const double diagonal   = Young/(1.0+Nu)/(1.0-2.0*Nu) * (1.0-Nu);
    const double nodiagonal = Young/(1.0+Nu)/(1.0-2.0*Nu) * ( Nu);
    const double corte      = Young/(1.0+Nu)/2.0;

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
                rElasticMatrix(i,j) = nodiagonal;
            }
        }
    }

    for (unsigned int j = 3; j<6; ++j)
        rElasticMatrix(j,j) = corte;

}

void MCPlasticFlowRule::CalculatePrincipalStressTrial(const RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix)
{
    Vector MainStrain      = ZeroVector(3);

    for (unsigned int i = 0; i<3; ++i)
    {
        MainStrain[i] = rNewElasticLeftCauchyGreen(i,i);
    }

    // Calculate the elastic matrix
    Matrix ElasticMatrix = ZeroMatrix(3,3);
    const double& Young     = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    const double& Nu        = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];
    const double diagonal   = Young/(1.0+Nu)/(1.0-2.0*Nu) * (1.0-Nu);
    const double nodiagonal = Young/(1.0+Nu)/(1.0-2.0*Nu) * ( Nu);

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
                ElasticMatrix(i,j) = nodiagonal;
            }
        }
    }

    Vector PrincipalStress = ZeroVector(3);

    // Evalute the Kirchhoff principal stress
    PrincipalStress = prod(ElasticMatrix, MainStrain);

    for(unsigned int i=0; i<3; i++)
    {
        rStressMatrix(i,i) = PrincipalStress(i);
    }
}
    

void MCPlasticFlowRule::ReturnStressFromPrincipalAxis(const Matrix& rEigenVectors, const Vector& rPrincipalStress, Matrix& rStressMatrix)
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

// this matrix is evaluated to make consistent the elastoplastic tangent matrix
// the principal stresses are in descending order
void MCPlasticFlowRule::CalculateModificationMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rAuxT, Matrix& rInvAuxT)
{
    if(mPrincipalStressUpdated(0) - mPrincipalStressUpdated(1) > 0)
    {
        rAuxT(0,0) = (mPrincipalStressUpdated(0) - mPrincipalStressUpdated(1))/(mPrincipalStressTrial(0) - mPrincipalStressTrial(1));
        rInvAuxT(0,0) = 1/rAuxT(0,0);
    }
    if(mPrincipalStressUpdated(0) - mPrincipalStressUpdated(2) > 0)
    {
        rAuxT(1,1) = (mPrincipalStressUpdated(0) - mPrincipalStressUpdated(2))/(mPrincipalStressTrial(0) - mPrincipalStressTrial(2));
        rInvAuxT(1,1) = 1/rAuxT(1,1);
    }
    if(mPrincipalStressUpdated(1) - mPrincipalStressUpdated(2) > 0)
    {
        rAuxT(2,2) = (mPrincipalStressUpdated(1) - mPrincipalStressUpdated(2))/(mPrincipalStressTrial(1) - mPrincipalStressTrial(2));
        rInvAuxT(2,2) = 1/rAuxT(2,2);
    }
}

void MCPlasticFlowRule::CalculateDepSurface(Matrix& rElasticMatrix, Vector& rFNorm, Vector& rGNorm, Matrix& rAuxDep)
{
    Vector AuxF = prod(trans(rFNorm), rElasticMatrix);

    Vector A = prod(rElasticMatrix, rGNorm);
    Vector B = prod(trans(rFNorm),rElasticMatrix);

    Matrix Num = ZeroMatrix(3,3);

    for (unsigned int i = 0; i<3; i++)
    {
        for(unsigned int j = 0; j<3; j++)
        {
            Num(i,j) = A(i) * B(j);
        }
    }

    double Den = MathUtils<double>::Dot(AuxF,rGNorm);
    rAuxDep = rElasticMatrix - Num / Den;

}

void MCPlasticFlowRule::CalculateDepLine(Matrix& rInvD, Vector& rFNorm, Vector& rGNorm, Matrix& rAuxDep)
{
    Matrix Num = ZeroMatrix(3,3);

    for (unsigned int i = 0; i<3; i++)
    {
        for(unsigned int j = 0; j<3; j++)
        {
            Num(i,j) = rFNorm(i) * rGNorm(j);
        }
    }

    Vector Den1 = prod(rInvD, rGNorm);
    double Den = MathUtils<double>::Dot(trans(rFNorm),Den1);

    rAuxDep = Num / Den;

}

void MCPlasticFlowRule::CalculateTransformationMatrix(const Matrix& rMainDirection, Matrix& rA)
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

void MCPlasticFlowRule::ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& alfa, Matrix& rConsistMatrix)
{
    //Elastoplastic constitutive matrix
    if (rReturnMappingVariables.Options.Is(MPMFlowRule::PLASTIC_REGION))
    {
        Matrix T = identity_matrix<double> (6);
        Matrix AuxT = ZeroMatrix(3,3);
        Matrix InvAuxT = ZeroMatrix(3,3);

        //1. Calculate  the modification matrix T
        this->CalculateModificationMatrix( rReturnMappingVariables, AuxT, InvAuxT);
        for(unsigned int i = 3; i<6; i++)
        {
            int index_i = i-3;

            T(i,i) = T(i,i) * AuxT(index_i, index_i);
        }

        Vector DeltaPrincipalStresses = mPrincipalStressTrial - mPrincipalStressUpdated;
        Matrix Dep = ZeroMatrix(6,6);

        //2. Calculate the ElastoPlastic Matrix depending on the region of return mapping
        this->CalculateElastoPlasticMatrix(rReturnMappingVariables, mRegion, DeltaPrincipalStresses, Dep);

        //3. Consistent Constitutive matrix
        Matrix DepcP = ZeroMatrix(6,6);
        DepcP = prod(T, Dep);

        //4. Return constitutive matrix from principal axis
        Matrix A = ZeroMatrix(6,6);
        Matrix ATrans = ZeroMatrix(6,6);
        
        //4.a Calculate transformation matrix
        this->CalculateTransformationMatrix(rReturnMappingVariables.MainDirections, A);
        ATrans = trans(A);

        Matrix AuxMat = ZeroMatrix(6,6);
        AuxMat = prod(ATrans, DepcP);
        rConsistMatrix = prod(AuxMat, A);
    }
    //Elastic matrix
    else
    {
        this->CalculateElasticMatrix(rReturnMappingVariables, rConsistMatrix);
    }

}

void MCPlasticFlowRule::CalculateElastoPlasticMatrix(const RadialReturnVariables& rReturnMappingVariables, unsigned int& rRegion, Vector& DiffPrincipalStress, Matrix& rDep)
{
    
    const double Young      = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    const double Nu         = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];
    const double shear_contribution = Young/(1.0+Nu)/2.0;

    const double FrictionAngle  = mMaterialParameters.FrictionAngle; 
    const double DilatancyAngle = mMaterialParameters.DilatancyAngle;

    const double FrictionCoefficient = (1 + std::sin(FrictionAngle))/(1 - std::sin(FrictionAngle));
    const double DilatancyCoefficient = (1 + std::sin(DilatancyAngle))/(1 - std::sin(DilatancyAngle));

    // Return mapping on yield surface
    if(rRegion == 1) 
    {
        // Yield plane normal
        Vector FNorm = ZeroVector(3);
        FNorm(0) = FrictionCoefficient;
        FNorm(1) = 0;
        FNorm(2) = -1;

        // Potential plane normal
        Vector GNorm = ZeroVector(3);
        GNorm(0) = DilatancyCoefficient;
        GNorm(1) = 0;
        GNorm(2) = -1;
        Matrix AuxDep = ZeroMatrix(3,3);

        // Compute elastic matrix which takes account only for normal stresses
        Matrix D = ZeroMatrix(3,3);
        this->ComputeElasticMatrix_3X3(rReturnMappingVariables, D); 

        this->CalculateDepSurface(D, FNorm, GNorm, AuxDep);

        // Shear components of consistent constitutive matrix in principal stress space
        for (unsigned int j = 3; j<6; ++j)
            rDep(j,j) = shear_contribution;

        for (unsigned int i = 0; i<3 ; ++i)
        {
            for (unsigned int k = 0; k<3 ; ++k)
            {
                rDep(i,k) = AuxDep(i,k);
            }
        }
    }
    // Return to a line 1, triaxial compression, sigp1 = sigp2
    else if(rRegion == 2) 
    {
        // Edge line direction
        Vector LFDir = ZeroVector(3);
        LFDir(0) = 1;
        LFDir(1) = 1;
        LFDir(2) = FrictionCoefficient;

        // Potential edge line direction
        Vector LGDir = ZeroVector(3);
        LGDir(0) = 1;
        LGDir(1) = 1;
        LGDir(2) = DilatancyCoefficient;

        Matrix InverseElasticMatrix = ZeroMatrix(3,3);
        this->CalculateInverseElasticMatrix(rReturnMappingVariables, InverseElasticMatrix);
        Matrix AuxDep = ZeroMatrix(3,3);

        this->CalculateDepLine(InverseElasticMatrix, LFDir, LGDir, AuxDep);

        for (unsigned int j = 3; j<6; ++j)
            rDep(j,j) = shear_contribution;

        for (unsigned int i = 0; i<3 ; ++i)
        {
            for (unsigned int k = 0; k<3 ; ++k)
            {
                rDep(i,k) = AuxDep(i,k);
            }
        }
    }
    // Return to a line 2, triaxial extension, sigp2 = sigp3
    else if(rRegion == 3) 
    {
        //Edge line direction
        Vector LFDir = ZeroVector(3);
        LFDir(0) = 1;
        LFDir(1) = FrictionCoefficient;
        LFDir(2) = FrictionCoefficient;

        //Potential edge line direction
        Vector LGDir = ZeroVector(3);
        LGDir(0) = 1;
        LGDir(1) = DilatancyCoefficient;
        LGDir(2) = DilatancyCoefficient;

        Matrix InverseElasticMatrix = ZeroMatrix(3,3);
        this->CalculateInverseElasticMatrix(rReturnMappingVariables, InverseElasticMatrix);
        Matrix AuxDep = ZeroMatrix(3,3);

        this->CalculateDepLine(InverseElasticMatrix, LFDir, LGDir, AuxDep);

        for (unsigned int j = 3; j<6; ++j)
            rDep(j,j) = shear_contribution;

        for (unsigned int i = 0; i<3 ; ++i)
        {
            for (unsigned int k = 0; k<3 ; ++k)
            {
                rDep(i,k) = AuxDep(i,k);
            }
        }
    }
}


Matrix MCPlasticFlowRule::GetElasticLeftCauchyGreen(RadialReturnVariables& rReturnMappingVariables)
{

    Vector Landa2 = ZeroVector(3);

    for (unsigned int i = 0; i<3; ++i)
        Landa2(i) = std::exp(2.0*mElasticPrincipalStrain(i));

    Matrix OutPut = ZeroMatrix(3,3);
    this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, Landa2, OutPut);

    return OutPut;
}

bool MCPlasticFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{
    // Compute Delta Plastic Strain
    double NormPlasticPrincipalStrain = norm_2(mPlasticPrincipalStrain);
    mInternalVariables.DeltaPlasticStrain = NormPlasticPrincipalStrain;

    // Compute Strain Components and its invariants
    double VolumetricPlasticPrincipalStrain = sum(mPlasticPrincipalStrain);
    Vector DeviatoricPlasticPrincipalStrain = mPlasticPrincipalStrain;
    for (unsigned int i = 0; i<3; ++i)
        DeviatoricPlasticPrincipalStrain(i) -= 1.0/3.0 * VolumetricPlasticPrincipalStrain;
    double DeltaAccumulatedPlasticDeviatoricStrain = sqrt(2.0/3.0) * norm_2(DeviatoricPlasticPrincipalStrain);

    const double FrictionAngle  = mMaterialParameters.FrictionAngle; 
    const double DilatancyAngle = mMaterialParameters.DilatancyAngle;

    const double FrictionCoefficient = (1 + std::sin(FrictionAngle))/(1 - std::sin(FrictionAngle));
    const double DilatancyCoefficient = (1 + std::sin(DilatancyAngle))/(1 - std::sin(DilatancyAngle));
    
    double NormStateFunctionDerivative = 0.0;

    // Calculate the norm of state function (or potential) derivative
    if(FrictionAngle == DilatancyAngle) // I am using an associative flow rule
    {
        NormStateFunctionDerivative = sqrt(1 + FrictionCoefficient * FrictionCoefficient);
    }
    else //I am using a non-associative flow rule
    {
        NormStateFunctionDerivative = sqrt(1 + DilatancyCoefficient * DilatancyCoefficient);
    }

    // Update Equivalent Plastic Strain
    double DeltaEquivalentPlasticStrain = mInternalVariables.DeltaPlasticStrain / NormStateFunctionDerivative;
    mInternalVariables.EquivalentPlasticStrain    += DeltaEquivalentPlasticStrain;

    // Update Accumulated Plastic Deviatoric Strain
    mInternalVariables.DeltaPlasticDeviatoricStrain = DeltaAccumulatedPlasticDeviatoricStrain;
    mInternalVariables.AccumulatedPlasticDeviatoricStrain += DeltaAccumulatedPlasticDeviatoricStrain;

    return true;
}

double MCPlasticFlowRule::GetSmoothingLodeAngle()
{
    return 29.0*GetPI()/180.0;
}


double MCPlasticFlowRule::GetPI()
{
    return std::atan(1.0)*4.0;
}

unsigned int MCPlasticFlowRule::GetPlasticRegion()
{
    return mRegion;
}


void MCPlasticFlowRule::ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH)
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
