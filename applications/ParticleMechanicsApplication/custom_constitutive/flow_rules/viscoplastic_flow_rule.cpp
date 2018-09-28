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

// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "particle_mechanics_application.h"
#include "viscoplastic_flow_rule.hpp"

namespace Kratos
{
//*******************************CONSTRUCTOR******************************************
//************************************************************************************

ViscoplasticFlowRule::ViscoplasticFlowRule()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

ViscoplasticFlowRule::ViscoplasticFlowRule(YieldCriterionPointer pYieldCriterion)
    :NonLinearAssociativePlasticFlowRule(pYieldCriterion)
{

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ViscoplasticFlowRule& ViscoplasticFlowRule::operator=(ViscoplasticFlowRule const& rOther)
{
    FlowRule::operator=(rOther);
    return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

ViscoplasticFlowRule::ViscoplasticFlowRule(ViscoplasticFlowRule const& rOther)
    :NonLinearAssociativePlasticFlowRule(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

FlowRule::Pointer ViscoplasticFlowRule::Clone() const
{
    FlowRule::Pointer p_clone(new ViscoplasticFlowRule(*this));
    return p_clone;
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

ViscoplasticFlowRule::~ViscoplasticFlowRule()
{
}

/// Operations.
///
//***************************UPDATE INTERNAL VARIABLES********************************
//************************************************************************************

bool ViscoplasticFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{

    mInternalVariables.EquivalentPlasticStrainOld  = mInternalVariables.EquivalentPlasticStrain;

    mInternalVariables.DeltaPlasticStrain          = std::sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

    mInternalVariables.EquivalentPlasticStrain    += mInternalVariables.DeltaPlasticStrain;

    mInternalVariables.DeltaPlasticStrain         *= ( 1.0/rReturnMappingVariables.DeltaTime );

    //update thermal variables
    mThermalVariables = rReturnMappingVariables.Thermal;

    // mInternalVariables.print();

    // mThermalVariables.print();

    return true;
}

//***************************UPDATE STRESS CONFIGURATION *****************************
//************************************************************************************

void ViscoplasticFlowRule::UpdateConfiguration( RadialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix )
{
    //Back Stress update

    //std::cout<< " ElasticIsoStress "<<rIsoStressMatrix<<std::endl;

    //Plastic Strain Update
    if( rReturnMappingVariables.NormIsochoricStress > 0 )
    {

        //Stress Update:
        double Auxiliar   = 2.0 * rReturnMappingVariables.LameMu_bar * rReturnMappingVariables.DeltaGamma;

        Matrix Normal     = rIsoStressMatrix * ( 1.0 / rReturnMappingVariables.NormIsochoricStress );

        rIsoStressMatrix -= ( Normal * Auxiliar );

    }

    //std::cout<< " PlasticIsoStress "<<rIsoStressMatrix<<std::endl;
}

//***************************CALCULATE LOCAL NEWTON PROCEDURE*************************
//************************************************************************************


bool ViscoplasticFlowRule::CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters )
{
    //Set aux variables
    //std::cout << GetProperties()[1][DYNAMIC_VISCOSITY] << std::endl;
    //std::cout << GetProperties()[VISCOSITY] << std::endl;
    //std::cout << GetProperties().Has(YOUNG_MODULUS) << std::endl;
    double mu = GetProperties().GetValue(VISCOSITY);
    //mu = 1.0;
    //double G = GetProperties().GetValue(YOUNG_MODULUS)/ (1.0 + mu) / 2.0;
    double dt = rReturnMappingVariables.DeltaTime;
    double RateSensitivity = 1;

//    std::cout << "VISCOSITY:" << mu << std::endl;
//    std::cout << "G:" << G << std::endl;
//    std::cout << "TimeStep:" << dt << std::endl;

    //Set convergence parameters
    unsigned int iter    = 0;
    double Tolerance     = 1e-5;
    double MaxIterations = 50;

    //start
    double DeltaDeltaGamma    = 0;
    double DeltaStateFunction = 0;
    rReturnMappingVariables.DeltaGamma    = 0;

    double Qtrail = rReturnMappingVariables.NormIsochoricStress;
    double G = rCriterionParameters.GetLameMu_bar();
    //double StateFunction                  = rReturnMappingVariables.TrialStateFunction;
    //double OldTrialStateFunction = rReturnMappingVariables.TrialStateFunction;

    HardeningLaw::Parameters NewHardeningParameters(rCriterionParameters.GetHardeningParameters());
    NewHardeningParameters.SetEquivalentPlasticStrain(rPlasticVariables.EquivalentPlasticStrain);
    NewHardeningParameters.SetDeltaGamma(rReturnMappingVariables.DeltaGamma);


    double Hardening = mpYieldCriterion->GetHardeningLaw().CalculateHardening(Hardening,NewHardeningParameters);
    double StateFunction = (Qtrail-2.0*G*rReturnMappingVariables.DeltaGamma)*std::pow((dt/(mu*rReturnMappingVariables.DeltaGamma+dt)),RateSensitivity)-std::sqrt(2.0/3.0) * Hardening;

//    std::cout << "Qtrail:" << Qtrail << std::endl;
//    std::cout << "OldTrialStateFunction:" << rReturnMappingVariables.TrialStateFunction << std::endl;
//    std::cout << "StateFunction:" << StateFunction << std::endl;
//    std::cin.get();

    while ( std::abs(StateFunction)>=Tolerance && iter<=MaxIterations)
    {
        //Calculate Delta State Function:
        //DeltaStateFunction = mpYieldCriterion->CalculateDeltaStateFunction( DeltaStateFunction, rCriterionParameters );
        double DeltaHardening = mpYieldCriterion->GetHardeningLaw().CalculateDeltaHardening(Hardening,NewHardeningParameters);
        DeltaStateFunction = -(2.0*G+RateSensitivity*mu*(Qtrail-2.0*G*rReturnMappingVariables.DeltaGamma)/(mu*rReturnMappingVariables.DeltaGamma+dt))*std::pow((dt/(mu*rReturnMappingVariables.DeltaGamma+dt)),RateSensitivity)-2.0/3.0 * DeltaHardening;
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
        StateFunction = (Qtrail-2.0*G*rReturnMappingVariables.DeltaGamma)*std::pow((dt/(mu*rReturnMappingVariables.DeltaGamma+dt)),RateSensitivity)-std::sqrt(2.0/3.0) * Hardening;

        iter++;
//        std::cout << "Iter:" << iter << std::endl;
//        std::cout << "Qtrail:" << Qtrail << std::endl;
//        std::cout << "OldTrialStateFunction:" << rReturnMappingVariables.TrialStateFunction << std::endl;
//        std::cout << "StateFunction:" << StateFunction << std::endl;
//        std::cout << "DeltaDeltaGamma:" << DeltaDeltaGamma << std::endl;
//        std::cout << "DeltaGamma:" << rReturnMappingVariables.DeltaGamma << std::endl;
//        std::cin.get();
    }


    if(iter>MaxIterations)
        return false;


    return true;
}

void ViscoplasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
}

void ViscoplasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )
}



}// namespace Kratos.

