//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
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
//NonAssociativePlasticFlowRule::NonAssociativePlasticFlowRule( const NonAssociativePlasticFlowRule & rOther)
//:FlowRule(rOther)
//{

//}
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
    //mTrialEigenValues = ZeroVector(3);
    //mEigenVectors = ZeroMatrix(3,3);
    mElasticPrincipalStrain = ZeroVector(3);
    mPlasticPrincipalStrain = ZeroVector(3);
    mElasticPreviousPrincipalStrain = ZeroVector(3);
    mPrincipalStressTrial = ZeroVector(3);
    mPrincipalStressUpdated = ZeroVector(3);
    mLargeStrainBool = true;
    mRegion = 0;

    mEquivalentPlasticStrain = 0.0;
}
bool MCPlasticFlowRule::CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{

    //std::cout<<" in CalculateReturnMapping "<<std::endl;

    bool PlasticityActive = false;
    rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);
    //1.-Compute Principal Axis

    //Computes the polar decomposition and Gives the PrincipalStrain
    // mTrialEigenValues and mEigenValues are computed and stored (crec que és una guarrada, pero ...)
    Vector MainStrain      = ZeroVector(3);
    // Matrix EigenVectors    = ZeroMatrix(3,3);
    // this-> ComputePrincipalAxisStrain(rReturnMappingVariables, rNewElasticLeftCauchyGreen, MainStrain, EigenVectors);
    //std::cout<<"rNewElasticLeftCauchyGreen"<<rNewElasticLeftCauchyGreen<<std::endl;
    for (unsigned int i = 0; i<3; ++i)
        MainStrain[i] = rNewElasticLeftCauchyGreen(i,i);
    ////std::cout<<"MainStrain before 1"<<MainStrain<<std::endl;

//***************************************************************************************************************************

    //double Young        = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    //double Nu           = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];

    ////double LameLambda         = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    ////double LameMu             =  YoungModulus/(2*(1+PoissonCoefficient));

    ////2.-Compute ElasticMatrix &  Trial State
    //Matrix ElasticMatrix = ZeroMatrix(3,3);

    //double diagonal   = Young/(1.0+Nu)/(1.0-2.0*Nu) * (1.0-Nu);
    //double nodiagonal = Young/(1.0+Nu)/(1.0-2.0*Nu) * ( Nu);

    //for (unsigned int i = 0; i<3; ++i) {
    //for (unsigned int j = 0; j<3; ++j) {
    //if (i == j) {
    //ElasticMatrix(i,i) = diagonal;
    //}
    //else {
    //ElasticMatrix(i,j) = nodiagonal;
    //}
    //}
    //}

    //Vector PrincipalStress = ZeroVector(3);
    //PrincipalStress = prod(ElasticMatrix, MainStrain);
//***************************************************************************************************************************
    //Compute elastic matrix which takes account only for normal stresses
    //Matrix D = ZeroMatrix(3,3);
    //this->ComputeElasticMatrix_3X3(rReturnMappingVariables, D); //check to input

    Vector PrincipalStress = ZeroVector(3);

    for(unsigned int i=0; i<3; i++)
    {

        PrincipalStress(i) = rStressMatrix(i,i);
    }


    //std::cout<<" MainStrain"<<MainStrain<<std::endl;
    //std::cout<<" ElasticMatrix"<<ElasticMatrix<<std::endl;
    //std::cout<<" PrincipalStress before sorting"<<PrincipalStress<<std::endl;
    //mPrincipalStressTrial = PrincipalStress;
    //Sorting principal stress (to check is in ascending or descending order)
    //std::cout<<" PrincipalStress before"<<PrincipalStress<<std::endl;
//**********************************************************************************
//**********SORTING EIGENVALUES AND EIGENVECTORS************************************
    Matrix PrincipalDirection1 = ZeroMatrix(3,1);
    Matrix PrincipalDirection2 = ZeroMatrix(3,1);
    Matrix PrincipalDirection3 = ZeroMatrix(3,1);

    for(unsigned int i=0; i<3; i++)
    {
        PrincipalDirection1(i,0) = rReturnMappingVariables.MainDirections(0,i);
    }
    for(unsigned int i=0; i<3; i++)
    {
        PrincipalDirection2(i,0) = rReturnMappingVariables.MainDirections(1,i);
    }
    for(unsigned int i=0; i<3; i++)
    {
        PrincipalDirection3(i,0) = rReturnMappingVariables.MainDirections(2,i);
    }

    if(PrincipalStress(0)<PrincipalStress(1))
    {
        std::swap(PrincipalStress(0),PrincipalStress(1));

        std::swap(MainStrain(0),MainStrain(1));

        Matrix TempMatrix = PrincipalDirection1;

        PrincipalDirection1 = PrincipalDirection2;

        PrincipalDirection2 = TempMatrix;
    }

    if(PrincipalStress(1)<PrincipalStress(2))
    {
        std::swap(PrincipalStress(1),PrincipalStress(2));

        std::swap(MainStrain(1),MainStrain(2));

        Matrix TempMatrix = PrincipalDirection2;

        PrincipalDirection2 = PrincipalDirection3;

        PrincipalDirection3 = TempMatrix;
    }

    if(PrincipalStress(0)<PrincipalStress(1))
    {
        std::swap(PrincipalStress(0),PrincipalStress(1));

        std::swap(MainStrain(0),MainStrain(1));

        Matrix TempMatrix = PrincipalDirection1;

        PrincipalDirection1 = PrincipalDirection2;

        PrincipalDirection2 = TempMatrix;
    }
    for(unsigned int i=0; i<3; i++)
    {
        rReturnMappingVariables.MainDirections(i,0) = PrincipalDirection1(i,0);
    }
    for(unsigned int i=0; i<3; i++)
    {
        rReturnMappingVariables.MainDirections(i,1) = PrincipalDirection2(i,0);
    }
    for(unsigned int i=0; i<3; i++)
    {
        rReturnMappingVariables.MainDirections(i,2) = PrincipalDirection3(i,0);
    }
    mPrincipalStressTrial = PrincipalStress;
    //Now the component of mElasticPrincipalStrain are sorted in the same way as MainStrain!!!!!!!!!
    mElasticPrincipalStrain = MainStrain;
    mElasticPreviousPrincipalStrain = MainStrain;
    //std::cout<<" PrincipalStress after sorting"<<PrincipalStress<<std::endl;
    //std::cout<<" PrincipalDirection after sorting"<<rReturnMappingVariables.MainDirections<<std::endl;
//**********************************************************************************
//**********************************************************************************
    //std::sort(PrincipalStress.begin(), PrincipalStress.end());
    //double MinStress = PrincipalStress(0);
    //double MaxStress = PrincipalStress(2);

    //PrincipalStress(0) = MaxStress;
    //PrincipalStress(2) = MinStress;

    //std::cout<<" PrincipalStress after sorting"<<PrincipalStress<<std::endl;
    //std::cout<<" MainStrain after sorting"<<MainStrain<<std::endl;
    //3.- Check for the yield Condition
    InternalVariables PlasticVariables = mInternalVariables;
    rReturnMappingVariables.TrialStateFunction = 0.0;
    rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PrincipalStress, PlasticVariables.EquivalentPlasticStrain );
    //std::cout<<" TrialStateFunctionipalStress "<<rReturnMappingVariables.TrialStateFunction<<std::endl;
    if (rReturnMappingVariables.TrialStateFunction <= 0.0)
    {
        PlasticVariables.DeltaPlasticStrain = 0;
        mRegion = 0;
        mPrincipalStressUpdated = PrincipalStress;
        PlasticityActive = false;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);
        //std::cout<<"ELASTIC REGIME "<<std::endl;
    }
    else
    {
        //std::cout<<"PLASTIC REGIME "<<std::endl;
        //std::cout<<"TrialStateFunction "<<rReturnMappingVariables.TrialStateFunction<<std::endl;
        //Matrix InverseElasticMatrix = ZeroMatrix(3,3);
        //CalculateInverseElasticMatrix(rReturnMappingVariables, InverseElasticMatrix);
        //std::cout<<" 	InverseElasticMatrix "<<InverseElasticMatrix<<std::endl;

        //mElasticPrincipalStrain = MainStrain;

        int Region = 0;
        Vector PrincipalStressUpdated = ZeroVector(3);

        bool converged = this->CalculateConsistencyCondition(rReturnMappingVariables, PrincipalStress, mElasticPrincipalStrain, Region, PrincipalStressUpdated);
        if (!converged)
            std::cout<<" Constit Law Did Not Converge "<<std::endl;


        mRegion = Region;
        mPrincipalStressUpdated = PrincipalStressUpdated;
        //std::cout<<" mPrincipalStressTrial in return mapping "<< mPrincipalStressTrial<<std::endl;
        //std::cout<<" mPrincipalStressUpdated in return mapping "<< mPrincipalStressUpdated<<std::endl;
        //std::cout<<" mPrincipalStress in return mapping "<< PrincipalStress<<std::endl;
        //std::cout<<" Region in return mapping "<< Region<<std::endl;
        //std::cout<<" mRegion in return mapping "<< mRegion<<std::endl;
        PlasticityActive = true;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,true);

    }

//*********************************************************************************************
    //4. Return stress matrix from principal axis
    //Matrix A = ZeroMatrix(6,6);
    //Matrix ATrans = ZeroMatrix(6,6);
    ////4.a Calculate transformation matrix
    //this->CalculateTransformationMatrix(rReturnMappingVariables.MainDirections, A);
    //ATrans = trans(A);

    ////Vector StressUpdated = ZeroVector(6);
    //Matrix AuxMat = ZeroMatrix(6,3);
    //for (unsigned int i = 0; i<6; i++)
    //{
    //for(unsigned int j = 0; j<3;j++)
    //{
    //AuxMat(i,j) = ATrans(i,j);
    //}
    //}

    //Vector AuxStress = ZeroVector(3);
    //AuxStress = prod(AuxMat, mPrincipalStressUpdated);

    //rStressMatrix = MathUtils<double>::StressVectorToTensor(AuxStress);
//*************************************************************************************************
    //rStressMatrix is the matrix of the updated stress in cartesian configuration
    this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, mPrincipalStressUpdated, rStressMatrix);
    //std::cout<<" rStressMatrix "<<rStressMatrix<<std::endl;
    //std::cout<<" rStressMatrix from principal axis"<<rStressMatrix<<std::endl;
    Vector DeltaPrincipalStress = PrincipalStress - mPrincipalStressUpdated; //check if the order of subtraction is correct
    //std::cout<<" PrincipalStress in return mapping "<< PrincipalStress<<std::endl;
    //std::cout<<" mPrincipalStressUpdated in return mapping "<< mPrincipalStressUpdated<<std::endl;
    //std::cout<<" DeltaPrincipalStress in return mapping "<< DeltaPrincipalStress<<std::endl;
    //Updated the PrincipalStrain vector
    Matrix InverseElasticMatrix = ZeroMatrix(3,3);
    this->CalculateInverseElasticMatrix(rReturnMappingVariables, InverseElasticMatrix);
    //std::cout<<"mElasticPrincipalStrain before "<<mElasticPrincipalStrain<<std::endl;



    //Delta plastic strain
    Vector PlasticStrain = prod( InverseElasticMatrix, DeltaPrincipalStress);

    //MainStrain -= PlasticStrain;
    //mElasticPrincipalStrain = MainStrain;

    //Now the component of mElasticPrincipalStrain are sorted in the same way as PlasticStrain!!!!!!!!!
    mElasticPrincipalStrain -= PlasticStrain;
    mPlasticPrincipalStrain = PlasticStrain;

    //std::cout<<"DeltaPrincipalStress "<<DeltaPrincipalStress<<std::endl;
    //std::cout<<"delta PlasticStrain "<<PlasticStrain<<std::endl;
    //std::cout<<"mElasticPrincipalStrain after"<<mElasticPrincipalStrain<<std::endl;

    //i'm saving the updated info in terms of principal strain and stress in these matrix
    //these information will be used for the evaluation of the second contribution in the
    //consistent tangent matrix
    for (unsigned int i=0; i<3; i++)
    {
        rReturnMappingVariables.StrainMatrix(i,i) = mElasticPrincipalStrain(i);
        rReturnMappingVariables.TrialIsoStressMatrix(i,i) = mPrincipalStressUpdated(i);
    }

    //std::cout<<"mPrincipalStressUpdated "<<mPrincipalStressUpdated<<std::endl;
    //std::cout<<"mElasticPrincipalStrain "<<mElasticPrincipalStrain<<std::endl;
    //std::cout<<"rReturnMappingVariables.StrainMatrix"<<rReturnMappingVariables.StrainMatrix<<std::endl;
    //std::cout<<"rReturnMappingVariables.TrialIsoStressMatrix"<<rReturnMappingVariables.TrialIsoStressMatrix<<std::endl;

    //PlasticVariables.DeltaPlasticStrain = PlasticStrain;

    //std::cout<<" mElasticPrincipalStrain "<<mElasticPrincipalStrain<<std::endl;
    //rNewElasticLeftCauchyGreen = this->GetElasticLeftCauchyGreen(rReturnMappingVariables);

    rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED,true);

    return PlasticityActive;


}

bool MCPlasticFlowRule::CalculateConsistencyCondition(RadialReturnVariables& rReturnMappingVariables, Vector& rPrincipalStress, Vector& rPrincipalStrain, int& region, Vector& rPrincipalStressUpdated)
{

    //std::cout<<"in the consistency condition"<<std::endl;
    //Calculate stress return in principal stress
    //The flow rule is written for associated and non-associated plasticity
    //if the internal friction angle is = to the dilatancy angle the plasticity is associated.

    double Cohesion = mpYieldCriterion->GetHardeningLaw().GetProperties()[COHESION];
    double FrictionAngle = mpYieldCriterion->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
    double DilatancyAngle = mpYieldCriterion->GetHardeningLaw().GetProperties()[INTERNAL_DILATANCY_ANGLE];

    //FrictionAngle *= GetPI() / 180.0;
    //DilatancyAngle *= GetPI() / 180.0;

    //k
    double FrictionCoefficient = (1 + std::sin(FrictionAngle))/(1 - std::sin(FrictionAngle));

    //m
    double DilatancyCoefficient = (1 + std::sin(DilatancyAngle))/(1 - std::sin(DilatancyAngle));

    //comp
    double CohesionCoefficient = 2 * Cohesion * sqrt(FrictionCoefficient);

    //Stress coordinate of the criterions apex
    double apex = CohesionCoefficient/(FrictionCoefficient-1);

    //Compute elastic matrix which takes account only for normal stresses
    Matrix D = ZeroMatrix(3,3);
    this->ComputeElasticMatrix_3X3(rReturnMappingVariables, D); //check to input

    //Compute the direction of the plastic return stress Rp
    Vector Rp = ZeroVector(3);
    double denp = FrictionCoefficient *(D(0,0)*DilatancyCoefficient - D(0,2)) - D(2,0) * DilatancyCoefficient + D(2,2);

    Rp(0) = (D(0,0)*DilatancyCoefficient - D(0,2) )/denp;
    Rp(1) = (D(1,0)*DilatancyCoefficient - D(1,2) )/denp;
    Rp(2) = (D(2,0)*DilatancyCoefficient - D(2,2) )/denp;

    //Vector from predictor stress to the apex
    Vector SigmaPApex = ZeroVector(3);

    SigmaPApex(0) = rPrincipalStress(0) - apex;
    SigmaPApex(1) = rPrincipalStress(1) - apex;
    SigmaPApex(2) = rPrincipalStress(2) - apex;

    //Boundary plane between region I and II: evaluated as the cross product between Rp and R1, the direction of line 1
    Vector NI_II = ZeroVector(3);
    NI_II(0) = Rp(1) * FrictionCoefficient - Rp(2);
    NI_II(1) = Rp(2) - Rp(0) * FrictionCoefficient;
    NI_II(2) = Rp(0) - Rp(1);

    double pI_II = NI_II(0) * SigmaPApex(0) + NI_II(1) * SigmaPApex(1) + NI_II(2) * SigmaPApex(2);

    //Boundary plane between region I and III: evaluated as the cross product between Rp and R2, the direction of line 2
    Vector NI_III = ZeroVector(3);
    NI_III(0) = Rp(1) * FrictionCoefficient - Rp(2) * FrictionCoefficient;
    NI_III(1) = Rp(2) - Rp(0) * FrictionCoefficient;
    NI_III(2) = Rp(0) * FrictionCoefficient - Rp(1);

    double pI_III = NI_III(0) * SigmaPApex(0) + NI_III(1) * SigmaPApex(1) + NI_III(2) * SigmaPApex(2);

    //t-paramters for region determination
    //secondary surface in region II a = [0 k -1], b  = [0 m -1]
    double denp2 = FrictionCoefficient * (D(1,1) * DilatancyCoefficient - D(1,2)) - D(1,2) * DilatancyCoefficient + D(2,2);

    Vector Rp2 = ZeroVector(3);
    Rp2(0) = (D(0,1) * DilatancyCoefficient - D(0,2))/denp2;
    Rp2(1) = (D(1,1) * DilatancyCoefficient - D(1,2))/denp2;
    Rp2(2) = (D(2,1) * DilatancyCoefficient - D(2,2))/denp2;

    //N2 = cross(Rp,Rp2)
    Vector N2 = ZeroVector(3);
    N2(0) = Rp(1)*Rp2(2) - Rp(2)*Rp2(1);
    N2(1) = Rp(2)*Rp2(0) - Rp(0)*Rp2(2);
    N2(2) = Rp(0)*Rp2(1) - Rp(1)*Rp2(0);

    double num1 = N2(0) * SigmaPApex(0) + N2(1) * SigmaPApex(1) + N2(2) * SigmaPApex(2);
    double den1 = N2(0) + N2(1) + FrictionCoefficient * N2(2);

    double t1 = num1 / den1 ;

    //secondary surface in region III a = [k -1 0], b  = [m -1 0]
    double den = FrictionCoefficient * (D(0,0) * DilatancyCoefficient - D(0,1)) - D(1,0) * DilatancyCoefficient + D(1,1);

    Vector Rp3 = ZeroVector(3);
    Rp3(0) = (D(0,0) * DilatancyCoefficient - D(0,1))/den;
    Rp3(1) = (D(1,0) * DilatancyCoefficient - D(1,1))/den;
    Rp3(2) = (D(2,0) * DilatancyCoefficient - D(2,1))/den;

    //N3 = cross(Rp,Rp3)
    Vector N3 = ZeroVector(3);
    N3(0) = Rp(1)*Rp3(2) - Rp(2)*Rp3(1);
    N3(1) = Rp(2)*Rp3(0) - Rp(0)*Rp3(2);
    N3(2) = Rp(0)*Rp3(1) - Rp(1)*Rp3(0);

    double num2 = N3(0) * SigmaPApex(0) + N3(1) * SigmaPApex(1) + N3(2) * SigmaPApex(2);
    double den2 = N3(0) + FrictionCoefficient * N3(1) + FrictionCoefficient * N3(2);

    double t2 = num2 / den2 ;

    //Region determination and update
    //std::cout<<"apex "<<apex<<std::endl;
    //std::cout<<"t1 "<<t1<<std::endl;
    //std::cout<<"t2 "<<t2<<std::endl;
    //std::cout<<"pI_II "<<pI_II<<std::endl;
    //std::cout<<"pI_III "<<pI_III<<std::endl;
    //return mapping to the apex
    if(t1 > 0 || t2 > 0) //check: both the conditions have to be satisfied
    {
        region = 4;
        //double StateFunction     = rReturnMappingVariables.TrialStateFunction;
        rPrincipalStressUpdated(0) = apex;
        rPrincipalStressUpdated(1) = apex;
        rPrincipalStressUpdated(2) = apex;
        //std::cout<<"region "<<region<<std::endl;
        //std::cout<<" rPrincipalStress "<<rPrincipalStress<<std::endl;
        //std::cout<<" StateFunction "<<StateFunction<<std::endl;
        //std::cout<<" rPrincipalStressUpdated "<<rPrincipalStressUpdated<<std::endl;
    }
    //return mapping to line 1       R1 = [1 1 k]
    else if(pI_II < 0)
    {
        region = 2;
        //double StateFunction     = rReturnMappingVariables.TrialStateFunction;
        rPrincipalStressUpdated(0) = t1 + apex;
        rPrincipalStressUpdated(1) = t1 + apex;
        rPrincipalStressUpdated(2) = t1 * FrictionCoefficient + apex;
        //std::cout<<"region "<<region<<std::endl;
        //std::cout<<" rPrincipalStress "<<rPrincipalStress<<std::endl;
        //std::cout<<" StateFunction "<<StateFunction<<std::endl;
        //std::cout<<" rPrincipalStressUpdated "<<rPrincipalStressUpdated<<std::endl;
    }

    //return mapping to the yield surface
    else if(pI_III <= 0)
    {
        region = 1;
        double StateFunction     = rReturnMappingVariables.TrialStateFunction;
        //std::cout<<"StateFunction "<<StateFunction<<std::endl;
        //std::cout<<"Rp "<<Rp<<std::endl;
        rPrincipalStressUpdated(0) = rPrincipalStress(0) -  StateFunction * Rp(0);
        rPrincipalStressUpdated(1) = rPrincipalStress(1) -  StateFunction * Rp(1);
        rPrincipalStressUpdated(2) = rPrincipalStress(2) -  StateFunction * Rp(2);
        //std::cout<<"region "<<region<<std::endl;
        //std::cout<<" rPrincipalStress "<<rPrincipalStress<<std::endl;
        //std::cout<<" StateFunction "<<StateFunction<<std::endl;
        //std::cout<<" rPrincipalStressUpdated "<<rPrincipalStressUpdated<<std::endl;
    }

    //return mapping to line 2
    else
    {
        region = 3;
        //double StateFunction     = rReturnMappingVariables.TrialStateFunction;
        rPrincipalStressUpdated(0) = t2 + apex;
        rPrincipalStressUpdated(1) = t2 * FrictionCoefficient + apex;
        rPrincipalStressUpdated(2) = t2 * FrictionCoefficient + apex;
        //std::cout<<"region "<<region<<std::endl;
        //std::cout<<" rPrincipalStress "<<rPrincipalStress<<std::endl;
        //std::cout<<" StateFunction "<<StateFunction<<std::endl;
        //std::cout<<" rPrincipalStressUpdated "<<rPrincipalStressUpdated<<std::endl;
        //std::cout<<"t1 "<<t1<<std::endl;
        //std::cout<<"pI_II "<<pI_II<<std::endl;
        //std::cout<<"pI_III "<<pI_III<<std::endl;
        //std::cout<<"rPrincipalStressUpdated "<<rPrincipalStressUpdated<<std::endl;
        //std::cout<<"apex "<<apex<<std::endl;
        //std::cout<<"t2 "<<t2<<std::endl;
        //std::cout<<"FrictionCoefficient "<<FrictionCoefficient<<std::endl;

    }
    //std::cout<<"rPrincipalStress in consistency condition "<<rPrincipalStress<<std::endl;
    //std::cout<<"rPrincipalStressUpdated in consistency condition"<<rPrincipalStressUpdated<<std::endl;
    //if (region != 0)
    //{
    //return true;
    //}
    return true;
}



void MCPlasticFlowRule::ComputeElasticMatrix_3X3(const RadialReturnVariables& rReturnMappingVariables, Matrix& rElasticMatrix)
{

    double YoungModulus        = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    double PoissonCoefficient  = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];



    double diagonal   = YoungModulus/(1.0+PoissonCoefficient)/(1.0-2.0*PoissonCoefficient) * (1.0-PoissonCoefficient);
    double nodiagonal = YoungModulus/(1.0+PoissonCoefficient)/(1.0-2.0*PoissonCoefficient) * ( PoissonCoefficient);

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

    //double LameLambda         = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    //double LameMu             =  YoungModulus/(2*(1+PoissonCoefficient));

    ////2.-Compute ElasticMatrix &  Trial State
    ////Matrix rElasticMatrix = ZeroMatrix(3,3);
    //for (unsigned int i = 0; i<3; ++i)
    //{
    //for (unsigned int j = 0; j<3; ++j)
    //rElasticMatrix(i,j) = LameLambda + 2.0/3.0*LameMu;

    //rElasticMatrix(i,i) += 2.0*LameMu;
    //}

}


void MCPlasticFlowRule::CalculateInverseElasticMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rInverseElasticMatrix)
{

    double YoungModulus        = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    double PoissonCoefficient  = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];

    double LameLambda         = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    double LameMu             =  YoungModulus/(2*(1+PoissonCoefficient));

    double Diagonal    = (LameLambda + LameMu)/(LameMu*(3.0*LameLambda+2.0*LameMu));
    double NonDiagonal = (-LameLambda)/( 2.0*LameMu*(3.0*LameLambda + 2.0*LameMu));

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

    //Matrix Aux = ZeroMatrix(6,6);
    //rElasticMatrix = Aux;

    double Young      = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    double Nu         = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];
    double diagonal   = Young/(1.0+Nu)/(1.0-2.0*Nu) * (1.0-Nu);
    double nodiagonal = Young/(1.0+Nu)/(1.0-2.0*Nu) * ( Nu);
    double corte      = Young/(1.0+Nu)/2.0;



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

    //std::cout<<"rElasticMatrix "<<rElasticMatrix<<std::endl;
}

void MCPlasticFlowRule::ReturnStressFromPrincipalAxis(const Matrix& rEigenVectors, const Vector& rPrincipalStress, Matrix& rStressMatrix)
{
    rStressMatrix = ZeroMatrix(3,3); //Esto seguro que está prohibido
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

//this matrix is evaluated to make consistent the elastoplastic tangent matrix
//the principal stresses are in descending order
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
    //Vector DeltaPrincipalStresses = mPrincipalStressTrial - mPrincipalStressUpdated;
    //rAuxT(0,0) = (DeltaPrincipalStresses(0) - DeltaPrincipalStresses(1))/(mPrincipalStressTrial(0) - mPrincipalStressTrial(1));
    //rInvAuxT(0,0) = 1/rAuxT(0,0);

    //rAuxT(1,1) = (DeltaPrincipalStresses(0) - DeltaPrincipalStresses(2))/(mPrincipalStressTrial(0) - mPrincipalStressTrial(2));
    //rInvAuxT(1,1) = 1/rAuxT(1,1);

    //rAuxT(2,2) = (DeltaPrincipalStresses(1) - DeltaPrincipalStresses(2))/(mPrincipalStressTrial(1) - mPrincipalStressTrial(2));
    //rInvAuxT(2,2) = 1/rAuxT(2,2);

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
    //std::cout<<" Num "<<Num<<std::endl;

    Vector Den1 = prod(rInvD, rGNorm);
    //std::cout<<" Den1 "<<Den1<<std::endl;
    //double den = prod(trans(rFNorm), Den1);
    double Den = MathUtils<double>::Dot(trans(rFNorm),Den1);
    //std::cout<<" Den "<<Den<<std::endl;
    rAuxDep = Num / Den;
    //std::cout<<" rAuxDep "<<rAuxDep<<std::endl;

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


//void MCPlasticFlowRule::GetPrincipalStressAndStrain(Vector& rPrincipalStresses, Vector& rPrincipalStrains)
//{
//rPrincipalStresses = mPrincipalStressUpdated;
//rPrincipalStrains = mElasticPreviousPrincipalStrain;
//}

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
        //std::cout<<"AuxT "<<AuxT<<std::endl;
        for(unsigned int i = 3; i<6; i++)
        {
            int index_i = i-3;

            T(i,i) = T(i,i) * AuxT(index_i, index_i);
        }
        //std::cout<<"T "<<T<<std::endl;
        Vector DeltaPrincipalStresses = mPrincipalStressTrial - mPrincipalStressUpdated;
        //std::cout<<"mPrincipalStressTrial in ept matrix "<<mPrincipalStressTrial<<std::endl;
        //std::cout<<"mPrincipalStressUpdated in ept matrix "<<mPrincipalStressUpdated<<std::endl;
        //std::cout<<"DeltaPrincipalStresses in ept matrix "<<DeltaPrincipalStresses<<std::endl;
        Matrix Dep = ZeroMatrix(6,6);

        //2. Calculate the ElastoPlastic Matrix depending on the region of return mapping
        this->CalculateElastoPlasticMatrix(rReturnMappingVariables, mRegion, DeltaPrincipalStresses, Dep);
        //std::cout<<"Dep "<<Dep<<std::endl;
        //3. Consistent Constitutive matrix
        Matrix DepcP = ZeroMatrix(6,6);
        DepcP = prod(T, Dep);
        //std::cout<<"DepcP "<<DepcP<<std::endl;
        //4. Return constitutive matrix from principal axis
        Matrix A = ZeroMatrix(6,6);
        Matrix ATrans = ZeroMatrix(6,6);
        //4.a Calculate transformation matrix
        this->CalculateTransformationMatrix(rReturnMappingVariables.MainDirections, A);
        //std::cout<<"A "<<A<<std::endl;
        ATrans = trans(A);
        //std::cout<<"ATrans "<<ATrans<<std::endl;

        //Vector StressUpdated = ZeroVector(6);
        Matrix AuxMat = ZeroMatrix(6,6);
        AuxMat = prod(ATrans, DepcP);
        rConsistMatrix = prod(AuxMat, A);
        //std::cout<<"rConsistMatrix "<<rConsistMatrix<<std::endl;
    }
//Elastic matrix
    else
    {

        //Matrix rConsistMatrix = ZeroMatrix(6,6);
        this->CalculateElasticMatrix(rReturnMappingVariables, rConsistMatrix);
        //std::cout<<"rConsistMatrix "<<rConsistMatrix<<std::endl;

    }
//if (rReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION))
    //{
    //Vector DeltaPrincipalStresses = mPrincipalStressTrial - mPrincipalStressUpdated;
    //Matrix Dep = ZeroMatrix(3,3);

    ////2. Calculate the ElastoPlastic Matrix depending on the region of return mapping
    //this->CalculateElastoPlasticMatrix(rReturnMappingVariables, mRegion, DeltaPrincipalStresses, Dep);
    //}
    //else
//{

    ////Matrix rConsistMatrix = ZeroMatrix(6,6);
    //this->ComputeElasticMatrix_3X3(rReturnMappingVariables, rConsistMatrix);
    ////std::cout<<"rConsistMatrix "<<rConsistMatrix<<std::endl;

//}

}

void MCPlasticFlowRule::CalculateElastoPlasticMatrix(const RadialReturnVariables& rReturnMappingVariables, int& rRegion, Vector& DiffPrincipalStress, Matrix& rDep)
{

    double Young      = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    double Nu         = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];
    double shear_contribution = Young/(1.0+Nu)/2.0;

    //double Cohesion = mpYieldCriterion->GetHardeningLaw().GetProperties()[COHESION];
    double FrictionAngle = mpYieldCriterion->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
    double DilatancyAngle = mpYieldCriterion->GetHardeningLaw().GetProperties()[INTERNAL_DILATANCY_ANGLE];

    //FrictionAngle *= GetPI() / 180.0;
    //DilatancyAngle *= GetPI() / 180.0;

    //k
    double FrictionCoefficient = (1 + std::sin(FrictionAngle))/(1 - std::sin(FrictionAngle));

    //m
    double DilatancyCoefficient = (1 + std::sin(DilatancyAngle))/(1 - std::sin(DilatancyAngle));

    //comp
    //double CohesionCoefficient = 2 * Cohesion * sqrt(FrictionCoefficient);

    if(rRegion == 1) //return mapping on yield surface
    {
        //std::cout<<"rRegion "<<rRegion<<std::endl;
        //Yield plane normal
        Vector FNorm = ZeroVector(3);
        FNorm(0) = FrictionCoefficient;
        FNorm(1) = 0;
        FNorm(2) = -1;
        //Potential plane normal
        Vector GNorm = ZeroVector(3);
        GNorm(0) = DilatancyCoefficient;
        GNorm(1) = 0;
        GNorm(2) = -1;
        Matrix AuxDep = ZeroMatrix(3,3);
        //Compute elastic matrix which takes account only for normal stresses
        Matrix D = ZeroMatrix(3,3);
        this->ComputeElasticMatrix_3X3(rReturnMappingVariables, D); //check to input

        this->CalculateDepSurface(D, FNorm, GNorm, AuxDep);

        //Shear components of consistent constitutive matrix in principal stress space
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
    else if(rRegion == 2) //return to a line 1, triaxial compression, sigp1 = sigp2
    {
        //std::cout<<"rRegion "<<rRegion<<std::endl;
        //Edge line direction
        Vector LFDir = ZeroVector(3);
        LFDir(0) = 1;
        LFDir(1) = 1;
        LFDir(2) = FrictionCoefficient;
        //std::cout<<" LFDir "<<LFDir<<std::endl;

        //Potential edge line direction
        Vector LGDir = ZeroVector(3);
        LGDir(0) = 1;
        LGDir(1) = 1;
        LGDir(2) = DilatancyCoefficient;
        //std::cout<<" LGDir "<<LGDir<<std::endl;
        Matrix InverseElasticMatrix = ZeroMatrix(3,3);
        this->CalculateInverseElasticMatrix(rReturnMappingVariables, InverseElasticMatrix);
        //std::cout<<" InverseElasticMatrix "<<InverseElasticMatrix<<std::endl;
        Matrix AuxDep = ZeroMatrix(3,3);

        this->CalculateDepLine(InverseElasticMatrix, LFDir, LGDir, AuxDep);
        //std::cout<<" AuxDep "<<AuxDep<<std::endl;
        //Shear components of consistent constitutive matrix in principal stress space
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
    else if(rRegion == 3) //return to a line 1, triaxial extension, sigp2 = sigp3
    {
        //std::cout<<"rRegion "<<rRegion<<std::endl;
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

        //Shear components of consistent constitutive matrix in principal stress space
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
    else //apex return
    {
        //std::cout<<"rRegion "<<rRegion<<std::endl;
        //Dep could be the zero matrix or DNorm
        //int alfa = 10000;

        //Matrix InverseElasticMatrix = ZeroMatrix(3,3);
        //this->CalculateInverseElasticMatrix(rReturnMappingVariables, InverseElasticMatrix);
        //Vector KoitDir = prod(InverseElasticMatrix, DiffPrincipalStress);

        //Matrix AuxDep = ZeroMatrix(3,3);
        ////Compute elastic matrix which takes account only for normal stresses
        //Matrix D = ZeroMatrix(3,3);
        //this->ComputeElasticMatrix_3X3(rReturnMappingVariables, D); //check to input

        //this->CalculateDepSurface(D, KoitDir, KoitDir, AuxDep);
        //for (unsigned int i = 0; i<3 ; ++i)
        //{
        //for (unsigned int k = 0; k<3 ; ++k)
        //{
        //rDep(i,k) = AuxDep(i,k) / alfa;
        //}
        //}
        //std::cout<<" Cep apex "<<rDep<<std::endl;
    }
}


Matrix MCPlasticFlowRule::GetElasticLeftCauchyGreen(RadialReturnVariables& rReturnMappingVariables)
{

    //std::cout<<" in GetElasticLeftCauchyGreen "<<std::endl;
    Vector Landa2 = ZeroVector(3);

    //std::cout<<" mElasticPrincipalStrain "<<mElasticPrincipalStrain<<std::endl;
    //std::cout<<" mLargeStrainBool "<<mLargeStrainBool<<std::endl;
    //if (mLargeStrainBool)
    //{
    for (unsigned int i = 0; i<3; ++i)
        Landa2(i) = std::exp(2.0*mElasticPrincipalStrain(i));
    //}

    //else {
    //Landa2 = mElasticPrincipalStrain;
    //}
    //std::cout<<" mElasticPrincipalStrain "<<mElasticPrincipalStrain<<std::endl;
    //std::cout<<" Landa2 in principal direction"<<Landa2<<std::endl;
    Matrix OutPut = ZeroMatrix(3,3);
    this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, Landa2, OutPut);
    //std::cout<<" be in cartesian direction "<<OutPut<<std::endl;
    return OutPut;

}

bool MCPlasticFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{

    //mInternalVariables.EquivalentPlasticStrainOld  = PlasticVariables.EquivalentPlasticStrain;
    //mInternalVariables.EquivalentPlasticStrainOld = 0.0;

    double NormPlasticPrincipalStrain = sqrt((mPlasticPrincipalStrain(0) * mPlasticPrincipalStrain(0) + mPlasticPrincipalStrain(1) * mPlasticPrincipalStrain(1) +mPlasticPrincipalStrain(2) * mPlasticPrincipalStrain(2) ));


    mInternalVariables.DeltaPlasticStrain          = NormPlasticPrincipalStrain;//sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

//Calculate the norm of state function (or potential) derivative

    double FrictionAngle = mpYieldCriterion->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
    double DilatancyAngle = mpYieldCriterion->GetHardeningLaw().GetProperties()[INTERNAL_DILATANCY_ANGLE];

    //FrictionAngle *= GetPI() / 180.0;
    //DilatancyAngle *= GetPI() / 180.0;

    //k
    double FrictionCoefficient = (1 + std::sin(FrictionAngle))/(1 - std::sin(FrictionAngle));
    //m
    double DilatancyCoefficient = (1 + std::sin(DilatancyAngle))/(1 - std::sin(DilatancyAngle));
    double NormStateFunctionDerivative = 0.0;
    //double CohesionCoefficient = 2 * Cohesion * sqrt(FrictionCoefficient);

    if(FrictionAngle == DilatancyAngle) // I am using an associative flow rule
    {
        NormStateFunctionDerivative = sqrt(1 + FrictionCoefficient * FrictionCoefficient);
    }
    else //I am using a non-associative flow rule
    {
        NormStateFunctionDerivative = sqrt(1 + DilatancyCoefficient * DilatancyCoefficient);
    }

    //double NormStateFunctionDerivative = mpYieldCriterion->CalculateNormYieldFunctionDerivative(NormStateFunctionDerivative);

    double DeltaEquivalentPlasticStrain = mInternalVariables.DeltaPlasticStrain / NormStateFunctionDerivative;

    mInternalVariables.EquivalentPlasticStrain    += DeltaEquivalentPlasticStrain;
    //if(mPlasticPrincipalStrain(0) != 0.0)
    //{
    //std::cout<<" mInternalVariables.EquivalentPlasticStrain "<<mInternalVariables.EquivalentPlasticStrain<<std::endl;
    //std::cout<<" NormPlasticPrincipalStrain "<<NormPlasticPrincipalStrain<<std::endl;
    //}

    //mInternalVariables.DeltaPlasticStrain         *= ( 1.0/rReturnMappingVariables.DeltaTime );

    return true;
}




double MCPlasticFlowRule::GetSmoothingLodeAngle()
{
    return 29.0*GetPI()/180.0;
    //return 22.8*GetPI()/180.0;
    //return 10.0*GetPI()/180.0;
}


double MCPlasticFlowRule::GetPI()
{
    return 3.14159265359;
}

double MCPlasticFlowRule::GetSmoothingHiperbolic()
{
    return 10.0;
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
