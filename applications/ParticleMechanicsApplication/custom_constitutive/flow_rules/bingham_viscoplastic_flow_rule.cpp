//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Duan Wenjie
//


// System includes
#include <string>
#include <iostream>

// External includes
#include "solid_mechanics_application_variables.h"

// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "particle_mechanics_application.h"
#include "bingham_viscoplastic_flow_rule.hpp"

namespace Kratos
{
//*******************************CONSTRUCTOR******************************************
//************************************************************************************

BinghamViscoplasticFlowRule::BinghamViscoplasticFlowRule()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

BinghamViscoplasticFlowRule::BinghamViscoplasticFlowRule(YieldCriterionPointer pYieldCriterion)
    :NonLinearAssociativePlasticFlowRule(pYieldCriterion)
{

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

BinghamViscoplasticFlowRule& BinghamViscoplasticFlowRule::operator=(BinghamViscoplasticFlowRule const& rOther)
{
    FlowRule::operator=(rOther);
    return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

BinghamViscoplasticFlowRule::BinghamViscoplasticFlowRule(BinghamViscoplasticFlowRule const& rOther)
    :NonLinearAssociativePlasticFlowRule(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

FlowRule::Pointer BinghamViscoplasticFlowRule::Clone() const
{
    FlowRule::Pointer p_clone(new BinghamViscoplasticFlowRule(*this));
    return p_clone;
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

BinghamViscoplasticFlowRule::~BinghamViscoplasticFlowRule()
{
}

//***************************CALCULATE DELTA GAMMA*************************
//************************************************************************************


bool BinghamViscoplasticFlowRule::CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters )
{
    //Set aux variables
    double eta = GetProperties().GetValue(VISCOSITY);
    //KRATOS_WATCH(eta);
    double dt = rReturnMappingVariables.DeltaTime;

    //Set convergence parameters
    unsigned int iter    = 0;
    double Tolerance     = 1e-5;
    double MaxIterations = 50;

    //start
    double DeltaDeltaGamma    = 0;
    double DeltaStateFunction = 0;
    rReturnMappingVariables.DeltaGamma    = 0;

    double Qtrial = rReturnMappingVariables.NormIsochoricStress;
    double G = rCriterionParameters.GetLameMu_bar();

    HardeningLaw::Parameters NewHardeningParameters(rCriterionParameters.GetHardeningParameters());
    NewHardeningParameters.SetEquivalentPlasticStrain(rPlasticVariables.EquivalentPlasticStrain);
    NewHardeningParameters.SetDeltaGamma(rReturnMappingVariables.DeltaGamma);

    double Hardening = mpYieldCriterion->GetHardeningLaw().CalculateHardening(Hardening,NewHardeningParameters);
    //double StateFunction = (Qtrail-2.0*G*rReturnMappingVariables.DeltaGamma)*pow((dt/(mu*rReturnMappingVariables.DeltaGamma+dt)),RateSensitivity)-sqrt(2.0/3.0) * Hardening;
    double StateFunction = Qtrial-std::sqrt(2.0/3.0) * Hardening;

//    std::cout << "Qtrail:" << Qtrail << std::endl;
//    std::cout << "OldTrialStateFunction:" << rReturnMappingVariables.TrialStateFunction << std::endl;
//    std::cout << "OldHardening:" << Qtrail - rReturnMappingVariables.TrialStateFunction << std::endl;
//    std::cout << "NewHardening:" << sqrt(2.0/3.0) * Hardening << std::endl;
//    std::cout << "StateFunction:" << StateFunction << std::endl;
//    std::cin.get();

    while ( std::abs(StateFunction)>=Tolerance && iter<=MaxIterations)
    {
        //Calculate Delta State Function:
        //DeltaStateFunction = mpYieldCriterion->CalculateDeltaStateFunction( DeltaStateFunction, rCriterionParameters );
        double DeltaHardening = mpYieldCriterion->GetHardeningLaw().CalculateDeltaHardening(Hardening,NewHardeningParameters);
        //DeltaStateFunction = -(2.0*G+RateSensitivity*mu*(Qtrail-2.0*G*rReturnMappingVariables.DeltaGamma)/(mu*rReturnMappingVariables.DeltaGamma+dt))*pow((dt/(mu*rReturnMappingVariables.DeltaGamma+dt)),RateSensitivity)-2.0/3.0 * DeltaHardening;
        DeltaStateFunction = -(eta/dt+2.0*G+2.0/3.0 * DeltaHardening);
        //Calculate DeltaGamma:
        DeltaDeltaGamma  = -StateFunction/DeltaStateFunction;
        rReturnMappingVariables.DeltaGamma += DeltaDeltaGamma;

        //Update Equivalent Plastic Strain:
        rPlasticVariables.DeltaPlasticStrain       = std::sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;
        rPlasticVariables.EquivalentPlasticStrain  = rPlasticVariables.EquivalentPlasticStrainOld + rPlasticVariables.DeltaPlasticStrain;

        //Calculate State Function:
        NewHardeningParameters.SetEquivalentPlasticStrain(rPlasticVariables.EquivalentPlasticStrain);
        NewHardeningParameters.SetDeltaGamma(rReturnMappingVariables.DeltaGamma);
        Hardening = mpYieldCriterion->GetHardeningLaw().CalculateHardening(Hardening,NewHardeningParameters);
        StateFunction = Qtrial-2.0*G*rReturnMappingVariables.DeltaGamma-std::sqrt(2.0/3.0) * Hardening-eta/dt*rReturnMappingVariables.DeltaGamma;

        iter++;
//        std::cout << "Iter:" << iter << std::endl;
//        std::cout << "Qtrail:" << Qtrail << std::endl;
//        std::cout << "StateFunction:" << StateFunction << std::endl;
//        std::cout << "DeltaStateFunction:" << DeltaStateFunction << std::endl;
//        std::cout << "DeltaDeltaGamma:" << DeltaDeltaGamma << std::endl;
//        std::cout << "DeltaGamma:" << rReturnMappingVariables.DeltaGamma << std::endl;
//        std::cin.get();
    }


    if(iter>MaxIterations)
        return false;


    return true;

}
void BinghamViscoplasticFlowRule::CalculateScalingFactors(const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors )
{

    //1.-Identity build
    Matrix IdentityMatrix       = identity_matrix<double> (3);

    //2.-Auxiliar matrices
    rScalingFactors.Normal      = rReturnMappingVariables.TrialIsoStressMatrix * ( 1.0 / rReturnMappingVariables.NormIsochoricStress );

    Matrix Norm_Normal          = prod( rScalingFactors.Normal, trans(rScalingFactors.Normal) );

    double Trace_Norm_Normal    = Norm_Normal( 0, 0 ) + Norm_Normal( 1, 1 )	+ Norm_Normal( 2, 2 );

    rScalingFactors.Dev_Normal  = Norm_Normal;
    rScalingFactors.Dev_Normal -= (1.0/3.0) * Trace_Norm_Normal * IdentityMatrix;


    //3.-Auxiliar constants
    double EquivalentPlasticStrain = mInternalVariables.EquivalentPlasticStrain + std::sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;
    double DeltaHardening = 0;

    double eta = GetProperties().GetValue(VISCOSITY);
    //double dt = rReturnMappingVariables.DeltaTime;

    HardeningLaw::Parameters HardeningParameters;
    HardeningParameters.SetTemperature(rReturnMappingVariables.Temperature);
    HardeningParameters.SetEquivalentPlasticStrain(EquivalentPlasticStrain);

    DeltaHardening = mpYieldCriterion->GetHardeningLaw().CalculateDeltaHardening( DeltaHardening, HardeningParameters );

    rScalingFactors.Beta0 = eta / (2 * rReturnMappingVariables.LameMu_bar) + 1.0 + DeltaHardening/(3.0 * rReturnMappingVariables.LameMu_bar);

    rScalingFactors.Beta1 = 2.0 * rReturnMappingVariables.LameMu_bar * rReturnMappingVariables.DeltaGamma / rReturnMappingVariables.NormIsochoricStress;

    rScalingFactors.Beta2 = ( ( 1.0 - ( 1.0 / rScalingFactors.Beta0 ) ) * (2.0/3.0) * rReturnMappingVariables.NormIsochoricStress * rReturnMappingVariables.DeltaGamma )/(rReturnMappingVariables.LameMu_bar) ;

    rScalingFactors.Beta3 = ( ( 1.0 / rScalingFactors.Beta0 ) - rScalingFactors.Beta1 + rScalingFactors.Beta2 );

    rScalingFactors.Beta4 = ( ( 1.0 / rScalingFactors.Beta0 ) - rScalingFactors.Beta1 ) * rReturnMappingVariables.NormIsochoricStress / ( rReturnMappingVariables.LameMu_bar ) ;

}



void BinghamViscoplasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
}

void BinghamViscoplasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )
}

}// namespace Kratos.
