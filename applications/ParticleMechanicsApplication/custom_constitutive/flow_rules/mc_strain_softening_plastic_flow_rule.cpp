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
#include<cmath>

// External includes
#include "includes/ublas_interface.h"

// Project includes
#include "custom_constitutive/flow_rules/mc_strain_softening_plastic_flow_rule.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

#include "particle_mechanics_application.h"

namespace Kratos
{



//************ CONSTRUCTOR ***********
MCStrainSofteningPlasticFlowRule::MCStrainSofteningPlasticFlowRule()
    :MCPlasticFlowRule()
{
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

MCStrainSofteningPlasticFlowRule::MCStrainSofteningPlasticFlowRule(YieldCriterionPointer pYieldCriterion)
    :MCPlasticFlowRule(pYieldCriterion)
{

}

//********* ASSIGMENT OPERATOR
MCStrainSofteningPlasticFlowRule& MCStrainSofteningPlasticFlowRule::operator=(MCStrainSofteningPlasticFlowRule const& rOther)
{
    MCPlasticFlowRule::operator=(rOther);
    return *this;

}



//********** COPY CONSTRUCTOR *********
MCStrainSofteningPlasticFlowRule::MCStrainSofteningPlasticFlowRule(MCStrainSofteningPlasticFlowRule const& rOther)
    :MCPlasticFlowRule(rOther)
{
}

//*******   CLONE ********
MPMFlowRule::Pointer MCStrainSofteningPlasticFlowRule::Clone() const
{
    MPMFlowRule::Pointer p_clone(new MCStrainSofteningPlasticFlowRule(*this));
    return p_clone;
}



// ********** DESTRUCTOR **************
MCStrainSofteningPlasticFlowRule::~MCStrainSofteningPlasticFlowRule()
{
}

bool MCStrainSofteningPlasticFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{
    MCPlasticFlowRule::UpdateInternalVariables( rReturnMappingVariables );

    this->UpdateMaterialParameters();

    return true;
}

void MCStrainSofteningPlasticFlowRule::UpdateMaterialParameters()
{
    // Calculate hardening for each parameters: cohesion, frictionangle, and dilatancyangle
    double Hardening;

    Hardening = mpYieldCriterion->GetHardeningLaw().CalculateHardening(Hardening, mInternalVariables.AccumulatedPlasticDeviatoricStrain, COHESION);
    Hardening = Hardening * mInternalVariables.DeltaPlasticDeviatoricStrain;
    mMaterialParameters.Cohesion += Hardening;

    Hardening = mpYieldCriterion->GetHardeningLaw().CalculateHardening(Hardening, mInternalVariables.AccumulatedPlasticDeviatoricStrain, INTERNAL_FRICTION_ANGLE);
    Hardening = Hardening * mInternalVariables.DeltaPlasticDeviatoricStrain;
    mMaterialParameters.FrictionAngle += Hardening;

    Hardening = mpYieldCriterion->GetHardeningLaw().CalculateHardening(Hardening, mInternalVariables.AccumulatedPlasticDeviatoricStrain, INTERNAL_DILATANCY_ANGLE);
    Hardening = Hardening * mInternalVariables.DeltaPlasticDeviatoricStrain;
    mMaterialParameters.DilatancyAngle += Hardening;
}



void MCStrainSofteningPlasticFlowRule::save( Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMFlowRule )
}

void MCStrainSofteningPlasticFlowRule::load( Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMFlowRule )

}

} //end namespace kratos
